from collections import defaultdict
import traceback, json
from itertools import product
from pathlib import Path
from subprocess import TimeoutExpired
from typing import Optional, List, Any, Callable
from zipfile import ZipFile

import boto3
import pandas as pd
import uvicorn
import pydna.all as pyd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fastapi import FastAPI, Response, Request, HTTPException
from fastapi.routing import APIRoute
from fastapi.responses import FileResponse, JSONResponse
from gfapy import Gfa
from mangum import Mangum
from mangum.types import LambdaEvent
from sqlmodel import Session, select
from starlette.middleware.exceptions import ExceptionMiddleware
from aws_lambda_powertools import Logger

if 'AppContainer' in __file__:
    from .denovo_from_paper import *
    from .helpers import add_files, plasmids_from_gfa
    from .annotation_propagation import align_sequences
    from .basespace import basespace_get
    from .db_model import make_engine, PlasmidSeqRun, ExperimentStartModel, PlasmidSeqRunList
else:
    from denovo_from_paper import *
    from helpers import add_files, plasmids_from_gfa
    from annotation_propagation import align_sequences
    from basespace import basespace_get
    from db_model import make_engine, PlasmidSeqRun, ExperimentStartModel, PlasmidSeqRunList

from Bio import SeqIO


class LoggerRouteHandler(APIRoute):
    def get_route_handler(self) -> Callable:
        original_route_handler = super().get_route_handler()

        async def route_handler(request: Request) -> Response:
            # Add fastapi context to logs
            ctx = {
                "path": request.url.path,
                "route": self.path,
                "method": request.method,
                "body": request.scope['aws.event']['body']
            }
            logger.append_keys(fastapi=ctx)
            logger.info("Received request")

            response = await original_route_handler(request)
            ctx["body"] = response.body
            logger.append_keys(fastapi=ctx)
            logger.info("Successfully completed request")

            return response

        return route_handler


app = FastAPI()
app.router.route_class = LoggerRouteHandler
app.add_middleware(ExceptionMiddleware, handlers=app.exception_handlers)

logger: Logger = Logger(log_uncaught_exceptions=True)

engine = make_engine()


def set_last_step(session: Session, job_data: PlasmidSeqRun, step_name: str) -> PlasmidSeqRun:
    existing_data = session.get(PlasmidSeqRun, job_data.id)

    if existing_data:
        job_data = existing_data
    else:
        session.add(job_data)

    job_data.last_step = step_name
    job_data.error_type = None
    job_data.error_message = None
    job_data.error_path = None
    session.commit()
    session.refresh(job_data)
    return job_data


@app.exception_handler(Exception)
async def unhandled_exception_handler(request: Request, err: Exception):
    logger.exception("Unhandled exception")

    event_body = request.scope['aws.event']['body']
    psr = PlasmidSeqRun.parse_raw(event_body)
    error_message = '; '.join(map(str, err.args))
    error_type = type(err).__name__
    if psr.id and psr.data_id and psr.experiment_id:
        with Session(engine) as session:
            db_psr = session.get(PlasmidSeqRun, psr.id)
            db_psr.last_step = 'Error'
            db_psr.error_type = error_type
            db_psr.error_message = error_message
            db_psr.error_path = request.url.path
            session.commit()

    return JSONResponse(status_code=500, content={"detail": "Internal Server Error",
                                                  "ExceptionType": error_type,
                                                  "ExceptionMessage": error_message})


@app.post('/start_experiment', response_model=PlasmidSeqRunList)
async def start_expt(exp_data: ExperimentStartModel):
    logger.info(exp_data.json())
    with Session(engine) as session:
        jobs = session.scalars(select(PlasmidSeqRun).where(PlasmidSeqRun.experiment_id == exp_data.name)).all()
        logger.info(jobs)
    run_list = PlasmidSeqRunList(runs=jobs)
    return run_list


@app.post('/setup', response_model=PlasmidSeqRun)
async def setup(job: PlasmidSeqRun):
    if job.last_step not in ['Queued', 'Error']:
        return job

    logger.info(f'{job.data_id} - Importing LabGuruAPI')
    from LabGuruAPI import Plasmid, SESSION
    SESSION.set_config_value('AUTH', 'TOKEN', os.environ['LG_AUTH_TOKEN'])
    tmp_folder = Path('/tmp') / job.data_path('setup_step')
    logger.info(f'{job.data_id} - Making Temp Folder')
    tmp_folder.mkdir(exist_ok=True, parents=True)

    # Get the fastq files
    logger.info(f'{job.data_id} - Querying Basespace')
    response = basespace_get(job.basespace_href)
    logger.info(response.content)
    logger.info(f'{job.data_id} - Downloading FastQ Files')
    fastq_folder = tmp_folder / 'fastq'
    fastq_folder.mkdir(exist_ok=True, parents=True)
    for cur_item in response.json()['Items']:
        logger.info(f'{job.data_id} - Downloading FastQ File {cur_item["Name"]}')
        file_response = basespace_get(cur_item['HrefContent'])
        out_file = fastq_folder / cur_item['Name']
        logger.info(f'Writing to {out_file.as_posix()}')
        out_file.write_bytes(file_response.content)

    logger.info(f'{job.data_id} - Writing FastQ to S3')
    add_files(job.data_path('raw_fastq'), *fastq_folder.iterdir())

    # get the template file
    logger.info(f'{job.data_id} - Downloading template gbk from LG')
    template_folder = tmp_folder / 'template_gbk'
    template_folder.mkdir(exist_ok=True, parents=True)
    gbk_file = template_folder / f"{job.template_name}.gb"
    plasmid = Plasmid.from_name(job.template_name)
    if not plasmid:
        logger.debug(f'{job.template_name} does not exist as a plasmid.')
        plasmid = await Plasmid.async_find_one(Plasmid.clone_no == job.template_name)
        if plasmid:
            with Session(engine) as session:
                j = session.get(PlasmidSeqRun, job.id)
                j.template_name = plasmid.name
                session.commit()
                session.refresh(j)
                job = j
        else:
            detail = f'{job.template_name} does not exist as a plasmid or clone ID.'
            raise ValueError(detail)

    if plasmid.sequence is not None:
        plasmid.sequence.write(str(gbk_file))
    else:
        raise ValueError(f'{job.template_name} does not have a sequence in LabGuru')

    logger.info(f'{job.data_id} - Writing Template to S3')
    add_files(job.data_path('template'), gbk_file)

    logger.info(f'{job.data_id} - Updating Job in RDS')
    with Session(engine) as session:
        job = set_last_step(session, job, 'Setup')
    return job


@app.post('/trim', response_model=PlasmidSeqRun)
async def trim_fastqs(job: PlasmidSeqRun):
    if job.last_step != 'Setup':
        return job
    # copy files to the temp folder
    tmp_folder = Path('/tmp') / job.data_path('trim_step')
    tmp_folder.mkdir(exist_ok=True, parents=True)

    s3 = boto3.resource('s3', region_name='us-east-1')
    bucket_resource = s3.Bucket('foundry-plasmid-seq')

    r1_file: Optional[Path] = None
    r2_file: Optional[Path] = None

    for obj in bucket_resource.objects.filter(Prefix=job.data_path('raw_fastq')):
        # Define the local path for the file
        target = tmp_folder / Path(obj.key).name

        # Download the file from S3 to the target path
        bucket_resource.download_file(obj.key, str(target))

        # Set the r1 or r2 file
        if '_R1_' in target.name:
            r1_file = target
        elif '_R2_' in target.name:
            r2_file = target

    assert r1_file is not None
    assert r2_file is not None

    dest_dir = tmp_folder / 'filtered'
    dest_dir.mkdir(exist_ok=True, parents=True)

    fwd_filtered, rev_filtered = trimmomatic(r1_file, r2_file, TRIMMOMATIC_MIN_LENGTH, DEFAULT_WINDOW,
                                             TRIMMOMATIC_DEFAULT_QUALITY, dest_dir, TRIMMOMATIC_PATH, True)

    add_files(job.data_path('trimmed_fastq'), Path(fwd_filtered), Path(rev_filtered))

    with Session(engine) as session:
        job = set_last_step(session, job, 'Trim')
    return job


@app.post('/check_transition/{from_state}/{to_state}', response_model=PlasmidSeqRun)
async def assembly_status(job: PlasmidSeqRun, from_state: str, to_state: str):
    with Session(engine) as session:
        job = session.get(PlasmidSeqRun, job.id)
    if job.last_step.lower() == from_state.lower():
        return JSONResponse(status_code=527,
                            content={"detail": f'Job {job.data_path()} is still in state {from_state}.'})
    elif job.last_step.lower() == to_state.lower():
        return job
    elif job.last_step == 'Error':
        return JSONResponse(status_code=500, content={"detail": "Internal Server Error",
                                                      "ExceptionType": job.error_type,
                                                      "ExceptionMessage": job.error_message})
    else:
        return job


@app.post('/assemble', response_model=PlasmidSeqRun)
async def assemble(job: PlasmidSeqRun):
    if job.last_step != 'Trim':
        return job

    # copy files to the temp folder
    tmp_folder = Path('/tmp') / job.data_path('assemble_step')
    result_folder = tmp_folder / 'assembly'
    tmp_folder.mkdir(exist_ok=True, parents=True)

    s3 = boto3.resource('s3', region_name='us-east-1')
    bucket_resource = s3.Bucket('foundry-plasmid-seq')

    r1_file: Optional[Path] = None
    r2_file: Optional[Path] = None

    for obj in bucket_resource.objects.filter(Prefix=job.data_path('trimmed_fastq')):
        # Define the local path for the file
        target = tmp_folder / Path(obj.key).name

        # Download the file from S3 to the target path
        bucket_resource.download_file(obj.key, str(target))

        # Set the r1 or r2 file
        if '_R1_' in target.name:
            r1_file = target
        elif '_R2_' in target.name:
            r2_file = target

    assert r1_file is not None
    assert r2_file is not None

    # queries, max_readlen = load_queries([r1_file, r2_file])

    try:
        assembly_fasta = denovo_assembly(r1_file, r2_file, result_folder, verbose=True)
    except OSError as ose:
        raise ose
    except TimeoutExpired:
        raise TimeoutError(f'{job.data_path()} could not be assembled in 14 min.')

    validate_file(assembly_fasta)

    assembly_gfa = Gfa.from_file(assembly_fasta.replace('.fasta', '.gfa'))
    assembly_seqs = plasmids_from_gfa(assembly_gfa)
    records = [pyd.Dseqrecord(s, id=f"Component_{i + 1}", circular=True) for i, s in
               enumerate(assembly_seqs)]
    SeqIO.write(records, assembly_fasta, 'fasta')

    if len(assembly_seqs) > 1:
        job.error_message = f"{len(assembly_seqs)} Assembly Options"

    add_files(job.data_path('assembly'), *result_folder.iterdir())
    with Session(engine) as session:
        job = set_last_step(session, job, 'Assemble')
    return job


@app.post('/transfer_annotations', response_model=PlasmidSeqRun)
async def transfer_annotations(job: PlasmidSeqRun):
    if job.last_step != 'Assemble':
        return job

    temp_dir = Path('/tmp') / job.data_path('analysis_step')
    out_dir = temp_dir / 'outputs'
    out_dir.mkdir(exist_ok=True, parents=True)

    s3 = boto3.resource('s3', region_name='us-east-1')
    bucket_resource = s3.Bucket('foundry-plasmid-seq')

    # Download template genbank
    for obj in bucket_resource.objects.filter(Prefix=job.data_path('template')):
        # Define the local path for the file
        target = temp_dir / Path(obj.key).name

        # Download the file from S3 to the target path
        if '.gb' in target.suffix:
            initial_file = target
            logger.info(f"Writing template to {initial_file.as_posix()}")
            bucket_resource.download_file(obj.key, str(target))
            break
    else:
        raise FileNotFoundError(f"Could not find a genbank file for template {job.template_name}.")

    # Download Assembly Files
    for obj in bucket_resource.objects.filter(Prefix=job.data_path('assembly')):
        # Define the local path for the file
        target = temp_dir / Path(obj.key).name

        # Download the file from S3 to the target path
        bucket_resource.download_file(obj.key, str(target))

    template_record: SeqRecord = SeqIO.read(str(initial_file), "genbank")
    annotated_assembly, position_map, mut_df, cds_df = align_sequences(str(initial_file),
                                                                       str(temp_dir / 'assembly.fasta'))

    assembly_file = out_dir / f'{job.data_id}.gbk'
    SeqIO.write([annotated_assembly], str(assembly_file), 'genbank')

    cds_list = []
    for cur_cds in (f for f in annotated_assembly.features if f.type == 'CDS'):
        cds_list.append({
            'Template': template_record.name,
            'Template CDS': cur_cds.qualifiers.get('template_label',
                                                   cur_cds.qualifiers.get('label', ['Unlabeled CDS']))[0],
            'Assembly CDS': cur_cds.qualifiers.get('label', ['Unlabeled CDS'])[0]
        })
    cds_summary = pd.DataFrame(cds_list).pivot(index='Template', columns='Template CDS', values='Assembly CDS')
    mut_df['Template'] = template_record.name
    cds_df['Template'] = template_record.name

    results_xlsx = out_dir / f'{template_record.name}_results.xlsx'
    with pd.ExcelWriter(results_xlsx) as xlsx:
        cds_summary.to_excel(xlsx, 'CDS Summary')
        mut_df.set_index('Template').to_excel(xlsx, 'Polymorphisms')
        cds_df.set_index('Template').to_excel(xlsx, 'CDS Changes')

    add_files(job.data_path('results/data'), results_xlsx)

    t_seq_list = list(str(template_record.seq))
    a_seq_list = list(str(annotated_assembly.seq))

    for t, a in position_map.items():
        a_seq_list[a] = t_seq_list[t]

    wt_like_sequence = ''.join(a_seq_list)
    annotated_assembly.seq = Seq(wt_like_sequence)

    wt_like_assembly_gbk = out_dir / f'{job.template_name}_with_{job.data_id}_indels.gbk'
    SeqIO.write([annotated_assembly], str(wt_like_assembly_gbk), 'genbank')
    add_files(job.data_path('results/sequences'), wt_like_assembly_gbk, assembly_file)

    with Session(engine) as session:
        job = set_last_step(session, job, 'Analysis')
    return job


@app.post('/combine', response_model=PlasmidSeqRunList)
async def combine_results(jobs: PlasmidSeqRunList):
    runs_by_type = defaultdict(list)
    for r in jobs.runs:
        r_type = type(r)
        runs_by_type[r_type].append(r)

    jobs.runs = runs_by_type[PlasmidSeqRun]
    for cur_list in runs_by_type[PlasmidSeqRunList]:
        jobs.runs.extend(cur_list.runs)

    if len(jobs.runs) == 0:
        raise ValueError("There are no runs to combine data for.")

    # Download Result Files
    experiment_id = jobs.runs[0].experiment_id
    result_dir = Path('/tmp') / experiment_id / 'results'
    temp_dir = result_dir / 'full'
    temp_dir.mkdir(parents=True, exist_ok=True)

    s3 = boto3.resource('s3', region_name='us-east-1')
    bucket_resource = s3.Bucket('foundry-plasmid-seq')

    for job in jobs.runs:
        for obj in bucket_resource.objects.filter(Prefix=job.data_path('results')):
            # Define the local path for the file
            target = temp_dir / Path(obj.key).name

            # Download the file from S3 to the target path
            logger.info(f'Downloading {target.name}')
            bucket_resource.download_file(obj.key, str(target))

    # Divide them by type
    zipfiles = {'genbanks': [], 'wt_like': [], 'data': [], 'other': []}
    for cur_file in temp_dir.iterdir():
        if '_indels.gbk' in cur_file.name:
            zipfiles['wt_like'].append(cur_file)
        elif cur_file.suffix == '.gbk':
            zipfiles['genbanks'].append(cur_file)
        elif cur_file.suffix == '.xlsx':
            zipfiles['data'].append(cur_file)
        else:
            zipfiles['other'].append(cur_file)

    # combine the data files
    logger.info('Combining Data Files')
    sheet_names = ['CDS Summary', 'Polymorphisms', 'CDS Changes']
    dataframes = defaultdict(list)
    for cur_sheet, cur_file in product(sheet_names, zipfiles['data']):
        cur_data = pd.read_excel(cur_file, cur_sheet)
        dataframes[cur_sheet].append(cur_data)

    results_xlsx = result_dir / f'{experiment_id}_full_results.xlsx'
    logger.info(f'Writing combined data to {results_xlsx.as_posix()}')
    with pd.ExcelWriter(results_xlsx) as xlsx:
        for cur_sheet, cur_data in dataframes.items():
            cur_df = pd.concat(cur_data)
            cur_df.to_excel(xlsx, cur_sheet)

    zip_path = result_dir / f'{experiment_id}_NGS_analysis.zip'
    with ZipFile(zip_path, 'w') as zf:
        zf.write(results_xlsx, results_xlsx.name)
        for z_folder, z_file_list in zipfiles.items():
            for z_file in z_file_list:
                zf.write(z_file, f'{z_folder}/{z_file.name}')

    add_files(f'{experiment_id}/results', results_xlsx, zip_path)

    with Session(engine) as session:
        for j in jobs.runs:
            db_j = session.get(PlasmidSeqRun, j.id)
            db_j.last_step = 'Complete'

        session.commit()

    return jobs


@app.get('/zipfile/{expt_id}')
async def download_results(expt_id: str):
    tmp_path = Path('/tmp') / expt_id / 'result_zip'
    tmp_path.mkdir(parents=True, exist_ok=True)
    tmp_file = tmp_path / f'{expt_id}_NGS_analysis.zip'

    s3 = boto3.resource('s3', region_name='us-east-1')
    bucket_resource = s3.Bucket('foundry-plasmid-seq')

    bucket_resource.download_file(f'{expt_id}/results/{expt_id}_NGS_analysis.zip', tmp_file)

    return FileResponse(tmp_file)


if __name__ == "__main__":
    uvicorn.run("app:handler")
else:
    handler = Mangum(app)
    handler = logger.inject_lambda_context(handler, clear_state=True)
