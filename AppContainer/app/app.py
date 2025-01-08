import os
import re
from collections import defaultdict
from pathlib import Path
from subprocess import TimeoutExpired
from typing import Optional, List, Callable, TypeAlias, Annotated
from zipfile import ZipFile

import boto3
import pandas as pd
import uvicorn
import pydna.all as pyd
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from botocore.exceptions import ClientError
from fastapi import FastAPI, Response, Request, Depends
from fastapi.routing import APIRoute
from fastapi.responses import FileResponse, JSONResponse
from gfapy import Gfa
from mangum import Mangum
from mypy_boto3_s3.service_resource import Bucket
from sqlmodel import Session, select
from starlette.middleware.exceptions import ExceptionMiddleware
from aws_lambda_powertools import Logger

if 'AppContainer' in __file__:
    from .denovo_from_paper import *
    from .helpers import add_files, plasmids_from_gfa, clean_genbank_file
    from .assembly_analysis import assembly_analysis_pipeline, AssemblyTooShortError
    from .basespace import basespace_get
    from .db_model import make_engine, PlasmidSeqRun, ExperimentStartModel, PlasmidSeqRunList, AssemblyList, \
        PlasmidSeqAssembly

    PRODUCTION = False
else:
    from denovo_from_paper import *
    from helpers import add_files, plasmids_from_gfa, clean_genbank_file
    from assembly_analysis import assembly_analysis_pipeline, AssemblyTooShortError
    from basespace import basespace_get
    from db_model import (make_engine, PlasmidSeqRun, ExperimentStartModel, PlasmidSeqRunList, AssemblyList,
                          PlasmidSeqAssembly)
    PRODUCTION = True

from Bio import SeqIO


# Folder constants
class S3Folders:
    RAW_FASTQ: str = 'raw_fastq'
    TEMPLATE: str = 'template'
    TRIMMED_FASTQ: str = 'trimmed_fastq'
    ASSEMBLY: str = 'assembly'
    PLASMIDS: str = 'assembly_plasmids'
    RESULTS: str = 'results'
    RESULT_SEQS: str = f'{RESULTS}/sequences'


class LoggerRouteHandler(APIRoute):
    def get_route_handler(self) -> Callable:
        original_route_handler = super().get_route_handler()

        async def route_handler(request: Request) -> Response:
            # Add fastapi context to logs
            ctx = {
                "path": request.url.path,
                "route": self.path,
                "method": request.method,
                "body": request.scope.get('aws.event', {}).get('body', 'DIRECT_CALL')
            }
            logger.append_keys(fastapi=ctx)
            logger.info("Received request")

            response = await original_route_handler(request)
            ctx["body"] = getattr(response, 'body', '')
            logger.append_keys(fastapi=ctx)
            logger.info("Successfully completed request")

            return response

        return route_handler


app = FastAPI()
app.router.route_class = LoggerRouteHandler
app.add_middleware(ExceptionMiddleware, handlers=app.exception_handlers)

logger: Logger = Logger(log_uncaught_exceptions=True)

engine = make_engine()


async def get_session() -> Session:
    new_session = Session(engine)
    try:
        yield new_session
    finally:
        new_session.close()


DSession: TypeAlias = Annotated[Session, Depends(get_session)]


async def pseq_bucket() -> Bucket:
    s3 = boto3.resource('s3', region_name='us-east-1')
    bucket_resource: Bucket = s3.Bucket('foundry-plasmid-seq')
    return bucket_resource


S3Bucket: TypeAlias = Annotated[Bucket, Depends(pseq_bucket)]


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


class WrongStateException(ValueError):
    def __init__(self, expected_status: str, job: PlasmidSeqRun, *args):
        self.expected_status = expected_status
        self.job = job
        super().__init__(*args)

    @staticmethod
    def check_job(job: PlasmidSeqRun, *allowed_statuses: str):
        if job.last_step not in allowed_statuses:
            raise WrongStateException(allowed_statuses[0] if len(allowed_statuses) == 1 else str(allowed_statuses), job)


@app.exception_handler(WrongStateException)
async def handle_wrong_state_exception(request: Request, err: WrongStateException):
    info_string = f"Job {err.job.data_path()} was in state {err.job.last_step}, {err.expected_status} was expected"
    logger.info(info_string)
    return JSONResponse(status_code=521, content=dict(detail=info_string))


@app.exception_handler(AssemblyTooShortError)
async def handle_too_short_error(request: Request, err: AssemblyTooShortError):
    info_string = ', '.join(err.args)
    logger.info(info_string)
    return JSONResponse(status_code=538, content=dict(detail=info_string))


@app.exception_handler(Exception)
async def unhandled_exception_handler(request: Request, err: Exception):
    logger.exception("Unhandled exception")

    error_message = '; '.join(map(str, err.args))
    error_type = type(err).__name__
    try:
        event_body = request.scope['aws.event']['body']
        psr = PlasmidSeqRun.parse_raw(event_body)
        if psr.id and psr.data_id and psr.experiment_id:
            with Session(engine) as session:
                db_psr = session.get(PlasmidSeqRun, psr.id)
                db_psr.last_step = 'Error'
                db_psr.error_type = error_type
                db_psr.error_message = error_message
                db_psr.error_path = request.url.path
                session.commit()
    except:
        pass

    return JSONResponse(status_code=500, content={"detail": "Internal Server Error",
                                                  "ExceptionType": error_type,
                                                  "ExceptionMessage": error_message})


@app.post('/start_experiment', response_model=PlasmidSeqRunList)
async def start_expt(exp_data: ExperimentStartModel):
    logger.info(exp_data.json())
    with Session(engine) as session:
        jobs: List[PlasmidSeqRun] = session.scalars(select(PlasmidSeqRun).where(PlasmidSeqRun.experiment_id == exp_data.name)).all()
        stripped_jobs = []
        for job in jobs:
            job.template_name = job.template_name.strip() if job.template_name else job.template_name
            job.data_id = job.data_id.strip() if job.data_id else job.data_id
            session.add(job)
            session.commit()
            session.refresh(job)
            logger.info(job)
            stripped_jobs.append(job.copy())
        logger.info(stripped_jobs)
    run_list = PlasmidSeqRunList(runs=stripped_jobs)
    return run_list


@app.post('/setup', response_model=PlasmidSeqRun)
async def setup(job: PlasmidSeqRun, session: DSession):
    if job.last_step not in ['Queued', 'Error']:
        return job

    logger.info(f'{job.data_id} - Importing LabGuruAPI')
    from LabGuruAPI import Plasmid, SESSION
    SESSION.set_config_value('AUTH', 'TOKEN', os.environ['LG_AUTH_TOKEN'])
    tmp_folder: Path = Path('/tmp') / job.data_path('setup_step')
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
    template_folder: Path = tmp_folder / 'template_gbk'
    template_folder.mkdir(exist_ok=True, parents=True)
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

    # Check for a foundry workflow template
    if plasmid.sequence is not None and (m := re.search(r'Clonal plasmid isolate of (.+)', plasmid.description)):
        parent_plasmid = Plasmid.from_name(m.group(1))
        if parent_plasmid:
            plasmid = parent_plasmid
            j = session.get(PlasmidSeqRun, job.id)
            j.template_name = plasmid.name
            session.commit()
            session.refresh(j)
            job = j
        else:
            linked_plasmids: List[Plasmid] = plasmid.get_linked_items(Plasmid)
            for lp in linked_plasmids:
                if 'GBFT' in lp.name:
                    plasmid = lp
                    j = session.get(PlasmidSeqRun, job.id)
                    j.template_name = plasmid.name
                    session.commit()
                    session.refresh(j)
                    job = j
                    break


    gbk_file: Path = template_folder / job.template_gb
    if plasmid.sequence is not None:
        out_seq = plasmid.sequence.copy()
        out_seq.name = plasmid.name[:16]
        out_seq.id = plasmid.name[:16]
        out_seq.locus = plasmid.name[:16]
        out_seq.write(str(gbk_file))
    else:
        raise ValueError(f'{job.template_name} does not have a sequence in LabGuru')

    logger.info(f'{job.data_id} - Writing Template to S3')
    add_files(job.data_path(S3Folders.TEMPLATE, 'bak'), gbk_file)
    await clean_genbank_file(gbk_file)  # clean up the genbank file from the terrible Geneious annotations
    add_files(job.data_path(S3Folders.TEMPLATE), gbk_file)

    logger.info(f'{job.data_id} - Updating Job in RDS')
    job = set_last_step(session, job, 'Setup')
    return job


@app.post('/trim', response_model=PlasmidSeqRun)
async def trim_fastqs(job: PlasmidSeqRun, session: DSession, bucket: S3Bucket):
    if job.last_step != 'Setup':
        return job
    # copy files to the temp folder
    tmp_folder = Path('/tmp') / job.data_path('trim_step')
    tmp_folder.mkdir(exist_ok=True, parents=True)

    r1_file: Optional[Path] = None
    r2_file: Optional[Path] = None

    for obj in bucket.objects.filter(Prefix=job.data_path(S3Folders.RAW_FASTQ)):
        # Define the local path for the file
        target = tmp_folder / Path(obj.key).name

        # Download the file from S3 to the target path
        bucket.download_file(obj.key, str(target))

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

    add_files(job.data_path(S3Folders.TRIMMED_FASTQ), Path(fwd_filtered), Path(rev_filtered))

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
async def assemble(job: PlasmidSeqRun, session: DSession, bucket: S3Bucket):
    if job.last_step != 'Trim':
        return job

    # copy fastq files to the temp folder
    tmp_folder = Path('/tmp') / job.data_path('assemble_step')
    result_folder = tmp_folder / 'assembly'
    tmp_folder.mkdir(exist_ok=True, parents=True)

    r1_file: Optional[Path] = None
    r2_file: Optional[Path] = None

    for obj in bucket.objects.filter(Prefix=job.data_path(S3Folders.TRIMMED_FASTQ)):
        # Define the local path for the file
        target = tmp_folder / Path(obj.key).name

        # Download the file from S3 to the target path
        bucket.download_file(obj.key, str(target))

        # Set the r1 or r2 file
        if '_R1_' in target.name:
            r1_file = target
        elif '_R2_' in target.name:
            r2_file = target

    assert r1_file is not None
    assert r2_file is not None

    # read in the template to find the rep origin
    template_file = tmp_folder / job.template_gb
    ori_file = tmp_folder / 'ori.fasta'

    bucket.download_file(job.data_path(S3Folders.TEMPLATE, job.template_gb), template_file.as_posix())
    template_record: SeqRecord = SeqIO.read(template_file.as_posix(), 'gb')

    cur_feature: SeqFeature
    for cur_feature in template_record.features:
        if cur_feature.type == 'rep_origin':
            f_name = cur_feature.qualifiers.get('label', ['origin'])[0]
            cur_feature.location.strand = 1
            f_seq = str(cur_feature.extract(template_record.seq)).upper()
            ori_file.write_text(f">{f_name}\n{f_seq}")
            unicycler_rotation_opts = ['--start_genes', ori_file.as_posix()]
            break
    else:
        unicycler_rotation_opts = ['--no_rotate']

    # run unicycler
    unicycler_debug = True
    unicycler_command = ['unicycler',
                         '-1', str(r1_file),
                         '-2', str(r2_file),
                         '-o', str(result_folder),
                         '--keep', str(3 if unicycler_debug else 1)]
    unicycler_command.extend(unicycler_rotation_opts)

    try:
        run_command(' '.join(unicycler_command))
        assembly_fasta = result_folder / 'assembly.fasta'
        assert assembly_fasta.is_file()
        add_files(job.data_path(S3Folders.ASSEMBLY), *result_folder.iterdir())
    except OSError as ose:
        raise ose
    except AssertionError:
        raise FileNotFoundError(f'{job.data_path()} did not result in an assembly file.')
    except TimeoutExpired:
        raise TimeoutError(f'{job.data_path()} could not be assembled in 14 min.')

    # Get assembly depth information
    spades_logfile = result_folder / 'spades_assembly' / 'spades.log'
    coverage_dict = {}
    cur_k = None
    for c_line in spades_logfile.read_text().splitlines(keepends=False):
        if cur_k is None and c_line.startswith('===== K') and 'started' in c_line:
            cur_k = c_line.split()[1]
        elif cur_k is not None and 'kmer_coverage_model.cpp   : 309' in c_line:
            coverage_data = c_line.split('   ')[-1]
            coverage_split = coverage_data.replace('std.', 'std').replace('. ', ': ').split(': ')
            coverage_dict[cur_k] = {'mean': float(coverage_split[1]), 'std': float(coverage_split[3])}
            cur_k = None

    unicycler_logfile = result_folder / 'unicycler.log'
    best_key = ''
    for c_line in unicycler_logfile.read_text().splitlines(keepends=False):
        if ' ← best' in c_line:
            best_key = 'K' + c_line.split()[0]
            break

    job = session.get(PlasmidSeqRun, job.id)
    job.assembly_coverage_mean = coverage_dict[best_key]['mean']
    job.assembly_coverage_std = coverage_dict[best_key]['std']
    session.commit()

    return set_last_step(session, job, 'Assemble')


@app.post('/make_plasmids', response_model=AssemblyList)
async def construct_plasmids(job: PlasmidSeqRun, bucket_resource: S3Bucket, session: DSession):
    job = session.get(PlasmidSeqRun, job.id)
    if job.last_step != 'Assemble':
        return AssemblyList(assemblies=job.assemblies)

    # copy files to the temp folder
    tmp_folder = Path('/tmp') / job.data_path('make_plasmid_step')
    tmp_folder.mkdir(exist_ok=True, parents=True)

    # Download the assembly graph
    logger.info('Reading assembly graph')
    graph_path = tmp_folder / 'assembly.gfa'
    bucket_resource.download_file(job.data_path('assembly/assembly.gfa'), graph_path.as_posix())
    assembly_graph = Gfa.from_file(graph_path.as_posix())

    # Build out sequences
    logger.info(f'Generating plasmid sequences from {str(assembly_graph)}')
    seq_info = plasmids_from_gfa(assembly_graph)
    if len(seq_info) == 0:
        raise ValueError(f'No plasmid sequences were generated from the GFA for {job.data_path()}')
    logger.info(f'{len(seq_info):d} valid paths were found. {str([i[1] for i in seq_info])}')

    # Download the template plasmid
    logger.info('Determining sequence directionality')
    template_path = tmp_folder / f'{job.template_name}.gb'
    bucket_resource.download_file(job.data_path(f'template/{job.template_name}.gb'), template_path.as_posix())
    template_seq: pyd.Dseqrecord = pyd.read(template_path.as_posix())

    # Determine if the sequences need to be reverse complemented
    temp_fwd_fragments: list[pyd.Dseq] = [template_seq.seq[i + 1:i + 10] for i in range(0, 51, 10)]
    temp_rev_fragments: list[pyd.Dseq] = [s.reverse_complement() for s in temp_fwd_fragments]

    # Create Sequence Records
    logger.info('Creating the assembly items')
    assembly_folder = tmp_folder / 'assemblies'
    assembly_folder.mkdir(exist_ok=True)

    assy_list = []
    all_records = []
    i = 0

    # Clear previous assemblies for this job
    prev_assemblies = session.query(PlasmidSeqAssembly).where(PlasmidSeqAssembly.run_id == job.id).all()
    pa: PlasmidSeqAssembly
    for pa in prev_assemblies:
        for ob in pa.features + pa.polymorphisms:
            session.delete(ob)
        session.delete(pa)

    for seq, path, prevalence in seq_info:
        i += 1
        cur_assy_obj = PlasmidSeqAssembly()
        cur_assy_obj.sample = job
        cur_assy_obj.assembly_name = f"{job.data_id}.{i:d}"
        cur_assy_obj.contig_path = path
        cur_assy_obj.length = len(seq)
        cur_assy_obj.min_prevalence = prevalence

        session.add(cur_assy_obj)
        session.commit()
        assy_list.append(cur_assy_obj)
        logger.info(f"Added {cur_assy_obj.assembly_name} to database, id={cur_assy_obj.id:d}")

        fwd_match_counts = sum(t_frag in seq for t_frag in temp_fwd_fragments)
        rev_match_counts = sum(t_frag in seq for t_frag in temp_rev_fragments)

        if fwd_match_counts > rev_match_counts:
            needs_rc = False
        elif rev_match_counts > fwd_match_counts:
            needs_rc = True
        elif fwd_match_counts == rev_match_counts:
            needs_rc = False  # May need to implement more logic here
        else:
            needs_rc = False  # Catch all for linting

        record_seq = seq.reverse_complement() if needs_rc else seq
        cur_record = pyd.Dseqrecord(record_seq, name=cur_assy_obj.assembly_name, circular=True,
                                    description=f'De novo assembly of {job.data_id}, path {i:d}')
        cur_record.write((assembly_folder / f"{cur_assy_obj.assembly_name}.gb").as_posix())
        all_records.append(cur_record)

    job = set_last_step(session, job, 'Split')
    job.assembly_count = len(assy_list)
    job.template_length = len(template_seq)
    session.add(job)
    session.commit()
    for pa in assy_list:
        session.refresh(pa)

    # Upload Assembly Files
    logger.info(f"Uploading {len(seq_info):d} plasmid maps to S3")
    SeqIO.write(all_records, assembly_folder / 'plasmids.fa', 'fasta')
    add_files(job.data_path('assembly_plasmids'), *assembly_folder.iterdir())

    return_obj = AssemblyList(assemblies=assy_list)
    logger.info(return_obj.json())
    return return_obj


@app.post('/transfer_annotations', response_model=PlasmidSeqRun)
async def transfer_annotations(assembly_obj: PlasmidSeqAssembly, bucket: S3Bucket, session: DSession):
    assembly_obj = session.get(PlasmidSeqAssembly, assembly_obj.id)
    job = assembly_obj.sample

    temp_dir = Path('/tmp') / job.data_path('analysis_step') / assembly_obj.assembly_name
    out_dir = temp_dir / 'outputs'
    out_dir.mkdir(exist_ok=True, parents=True)

    # Download template & assembly genbanks
    template_file = temp_dir / f"{job.template_name}.gb"
    bucket.download_file(job.data_path('template', job.template_name + '.gb'), template_file.as_posix())

    assembly_file = temp_dir / f"{assembly_obj.assembly_name}.gb"
    bucket.download_file(job.data_path('assembly_plasmids', assembly_file.name), assembly_file.as_posix())

    assembly_record = assembly_analysis_pipeline(template_file, assembly_file, assembly_obj, session)
    assembly_file = out_dir / assembly_obj.assembly_gb
    SeqIO.write([assembly_record], str(assembly_file), 'genbank')

    add_files(job.data_path(S3Folders.RESULT_SEQS), assembly_file)

    job = set_last_step(session, job, 'Analysis')
    return job


@app.post('/combine', response_model=PlasmidSeqRunList)
async def combine_results(jobs: PlasmidSeqRunList, bucket: S3Bucket, session: DSession):
    job_runs = jobs.runs_of_type(PlasmidSeqRun)

    if len(job_runs) == 0:
        raise ValueError("There are no runs to combine data for.")

    zipfiles: defaultdict[str, list[Path]] = defaultdict(list)

    # Download Genbank Files

    experiment_id = job_runs[0].experiment_id
    result_dir = Path('/tmp') / experiment_id / 'results'
    temp_dir = result_dir / 'full'
    temp_dir.mkdir(parents=True, exist_ok=True)

    assembly_data_dicts = []
    feature_data_dicts = []
    poly_data_dicts = []
    error_dicts = []
    all_experiment_runs: List[PlasmidSeqRun] = (session.query(PlasmidSeqRun)
                                                .where(PlasmidSeqRun.experiment_id == experiment_id).all())
    completed_jobs: List[PlasmidSeqRun] = []
    for cur_job in all_experiment_runs:
        # cur_job = session.get(PlasmidSeqRun, cur_job.id)
        if cur_job.last_step == 'Error':
            error_dicts.append(cur_job.dict(exclude=dict(id=True, last_step=True, basespace_href=True,
                                                         assembly_count=True, experiment_id=True,
                                                         template_length=True, assembly_coverage_mean=True,
                                                         assembly_coverage_std=True)))
            continue

        cur_job_dict = dict(template=cur_job.template_name, sample=cur_job.data_id,
                            assembly_count=cur_job.assembly_count, expected_bp=cur_job.template_length,
                            coverage_mean=cur_job.assembly_coverage_mean, coverage_std=cur_job.assembly_coverage_std)
        logger.info(cur_job_dict)
        has_valid_assemblies = False
        for cur_assembly in cur_job.assemblies:
            logger.info(f'Gathering information on {cur_assembly.assembly_name}')
            cur_gb_file = temp_dir / cur_assembly.assembly_gb
            gb_file_path = cur_job.data_path(S3Folders.RESULT_SEQS, cur_assembly.assembly_gb)
            try:
                bucket.download_file(gb_file_path, cur_gb_file.as_posix())
            except ClientError:
                if cur_assembly.length < cur_job.template_length/2:
                    error_dicts.append(dict(data_id=cur_job.data_id, template_name=cur_job.template_name,
                                            error_type='ValueError', error_path='/transfer_annotations',
                                            error_message=f"{cur_assembly.assembly_name} is too short to properly align"
                                                          f" / bp:{cur_assembly.length:d}, "
                                                          f"expected:{cur_job.template_length}"))
                    continue
                logger.warning(f'Could not find {gb_file_path} on S3. Retrying.')
                try:
                    bucket.download_file(gb_file_path, cur_gb_file.as_posix())
                except ClientError:
                    error_dicts.append(dict(data_id=cur_job.data_id, template_name=cur_job.template_name,
                                            error_type='ClientError', error_path='/combine',
                                            error_message=f"The genbank file for {cur_assembly.assembly_name} was not "
                                                          f"found on S3"))
                    continue
            zipfiles['genbanks'].append(cur_gb_file)

            cur_assembly_dict = cur_assembly.dict(exclude=dict(id=True, run_id=True, contig_path=True))
            assembly_name_dict = {'assembly': cur_assembly_dict.pop('assembly_name')}
            cur_assembly_dict['polymorphism_count'] = 0

            if len(cur_assembly.features) == 0:
                error_dicts.append(dict(data_id=cur_job.data_id, template_name=cur_job.template_name,
                                        error_type='Unknown', error_path='/transfer_annotations',
                                        error_message=f"No features were logged for {cur_assembly.assembly_name}. This"
                                                      f"usually indicates an error during annotation."))
                continue

            cur_feature_mods = []
            for cur_feature in cur_assembly.features:
                cur_feature_dict = cur_feature.dict(exclude=dict(id=True, assembly_id=True))
                cur_feature_dict['modified'] = (cur_feature.wt_feature_name != cur_feature.assembly_feature_name
                                                or cur_feature.deleted)
                if cur_feature_dict['modified']:
                    cur_feature_mods.append(cur_feature.wt_feature_name)
                feature_data_dicts.append({**cur_job_dict, **assembly_name_dict, **cur_feature_dict})
            cur_assembly_dict['modified_features'] = ', '.join(cur_feature_mods)

            for cur_poly in cur_assembly.polymorphisms:
                cur_poly_dict = cur_poly.dict(exclude=dict(id=True, assembly_id=True))
                cur_poly_dict['features'] = ', '.join(f.wt_feature_name for f in cur_poly.features)
                poly_data_dicts.append({**cur_job_dict, **assembly_name_dict, **cur_poly_dict})
                cur_assembly_dict['polymorphism_count'] += 1

            assembly_data_dicts.append({**cur_job_dict, **assembly_name_dict, **cur_assembly_dict})
            has_valid_assemblies = True

        if has_valid_assemblies:
            completed_jobs.append(cur_job)
        else:
            err_job = session.get(PlasmidSeqRun, cur_job.id)
            err_job.last_step = 'Error'
            err_job.error_type = 'Assembly'
            err_job.error_message = 'The job contained no valid assemblies'
            err_job.error_path = '/combine'
            session.add(err_job)

    # Generate Summary dataframes
    feature_data_df = pd.DataFrame(feature_data_dicts)
    feature_summary = feature_data_df.groupby('assembly').sum()[['deleted', 'modified']]
    poly_data_df = pd.DataFrame(poly_data_dicts)
    poly_summary = poly_data_df.groupby(['assembly', 'poly_type'])['sample'].count().reset_index()
    poly_pivot = poly_summary.pivot(columns='poly_type', index='assembly', values='sample').fillna(0)

    assembly_data_df = pd.DataFrame(assembly_data_dicts).sort_values(['template', 'assembly']).set_index('assembly')
    assembly_data_df['bp_diff'] = abs(assembly_data_df['expected_bp'] - assembly_data_df['length'])
    assembly_data_df.loc[feature_summary.index, 'del_feature_count'] = feature_summary['deleted']
    assembly_data_df.loc[feature_summary.index, 'mod_feature_count'] = feature_summary['modified']
    assembly_data_df.loc[poly_pivot.index, 'indel_count'] = poly_pivot.get('Insertion', 0) + poly_pivot.get('Deletion', 0)

    assembly_data_df = assembly_data_df.reset_index(drop=False)
    assembly_data_df = assembly_data_df.reindex(['template', 'sample', 'assembly_count', 'assembly',
                                                 'expected_bp', 'length', 'bp_diff', 'min_prevalence',
                                                 'coverage_mean', 'coverage_std', 'polymorphism_count',
                                                 'indel_count', 'del_feature_count', 'mod_feature_count',
                                                 'modified_features'], axis=1)  # final col order

    # write the data
    results_xlsx = result_dir / f'{experiment_id}_full_results.xlsx'
    logger.info(f'Writing combined data to {results_xlsx.as_posix()}')
    with pd.ExcelWriter(results_xlsx) as xlsx:
        assembly_data_df.to_excel(xlsx, 'Summary', index=False)
        feature_data_df.drop_duplicates().to_excel(xlsx, 'Features', index=False)
        poly_data_df.drop_duplicates().to_excel(xlsx, 'Polymorphisms', index=False)
        pd.DataFrame(error_dicts).to_excel(xlsx, 'Errors', index=False)

    zip_path = result_dir / f'{experiment_id}_NGS_analysis.zip'
    with ZipFile(zip_path, 'w') as zf:
        zf.write(results_xlsx, results_xlsx.name)
        for z_folder, z_file_list in zipfiles.items():
            for z_file in z_file_list:
                zf.write(z_file, f'{z_folder}/{z_file.name}')

    add_files(f'{experiment_id}/results', results_xlsx, zip_path)

    for j in completed_jobs:
        j = session.get(PlasmidSeqRun, j.id)
        j.last_step = 'Complete'
        session.add(j)

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


@app.get('/zipfile/{expt_id}/templates')
async def download_templates(expt_id: str, bucket: S3Bucket, session: DSession):
    tmp_path = Path('/tmp') / expt_id / 'templates' / 'genbanks'
    tmp_path.mkdir(parents=True, exist_ok=True)
    zip_path = tmp_path.parent / f'{expt_id}_template_gbs.zip'

    all_experiment_runs: List[PlasmidSeqRun] = (session.query(PlasmidSeqRun)
                                                .where(PlasmidSeqRun.experiment_id == expt_id).all())

    for cur_run in all_experiment_runs:
        template_path = cur_run.data_path(S3Folders.TEMPLATE, cur_run.template_gb)
        print(f'{template_path=}')
        template_gb_path = tmp_path / cur_run.template_gb
        bucket.download_file(template_path, template_gb_path)

        # Fix the record
        await clean_genbank_file(template_gb_path)

    with ZipFile(zip_path, 'w') as zf:
        for cur_file in tmp_path.iterdir():
            print(f'{cur_file.name=}')
            zf.write(cur_file, cur_file.name)

    print(f'{zip_path=}')
    return FileResponse(zip_path, filename=zip_path.name)


if __name__ == "__main__":
    uvicorn.run("app:handler")
else:
    handler = Mangum(app)
    handler = logger.inject_lambda_context(handler, clear_state=True)
