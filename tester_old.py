import json
from pathlib import Path
from subprocess import run

import requests
from sqlmodel import Session, SQLModel

from AppContainer.app.db_model import PlasmidSeqRun, make_engine, ExperimentStartModel, PlasmidSeqRunList
from AppContainer.app.helpers import upload_files, add_files

test_json = {
    "body": "{}",
    "resource": "/",
    "path": "/up",
    "httpMethod": "POST",
    "requestContext": {}
}
filter_json = {
    "body": "{}",
    "resource": "/",
    "path": "/trim/OW4QRniRLSGten9_JZFkv",
    "httpMethod": "GET",
    "requestContext": {}
}
assemble_json = {
    "body": "{}",
    "resource": "/",
    "path": "/assemble/Qmivz9EPOHku8c99k509K",
    "httpMethod": "GET",
    "requestContext": {}
}
gff_json = {
    "body": "{}",
    "resource": "/",
    "path": "/to/gff",
    "httpMethod": "POST",
    "requestContext": {}
}
align_json = {
    "body": "{}",
    "resource": "/",
    "path": "/transfer_annotations/qT22qXTOSZ1Q0uMr7cAWi",
    "httpMethod": "POST",
    "requestContext": {}
}


def full_pipeline(lambda_endpoint: str, fastq_dir: Path, template_path: Path,
                  fastq_folder: str = None, filter_folder: str = None, assembly_folder: str = None):
    fastq_folder = fastq_folder or upload_files(*fastq_dir.iterdir())

    # Run Filter
    if not filter_folder:
        resp1 = requests.get(lambda_endpoint, json=dict(
            body='{}',
            resource='/',
            path=f'/trim/{fastq_folder}',
            httpMethod='GET',
            requestContext={}
        ))
        print(resp1.content)
        filter_folder = json.loads(resp1.json()['body'])['filtered_folder']

    # Assemble
    if not assembly_folder:
        resp2 = requests.get(lambda_endpoint, json=dict(
            body='{}',
            resource='/',
            path=f'/assemble/{filter_folder}',
            httpMethod='GET',
            requestContext={}
        ))
        print(resp2.content)
        assembly_folder = json.loads(resp2.json()['body'])['assembly_folder']

    # Align
    resp3 = requests.get(lambda_endpoint, json=dict(
        body=json.dumps(dict(value=template_path.read_text())),
        resource='/',
        path=f'/transfer_annotations/{assembly_folder}',
        httpMethod='POST',
        requestContext={}
    ))
    print(resp3.content)


def init():
    fastq_path = Path(r'C:\BaseSpace\Expt1295-393806414\FASTQ_Generation_2023-07-22_08_39_25Z-683268586'
                      r'\GBFP-1295-0047_L001-ds.13fed8fe790f47b0bb67165fe0a4f612')

    gb_path = Path(r"C:\Users\RobertWarden-Rothman\GRO Biosciences\Projects - Foundry\Collaborations"
                   r"\1295 - CX - intein golden gate assembly and transformations\base plasmids\pGRO-K1196.gb")

    template_name = fastq_path.name.split('_')[0]
    cur_pseq = PlasmidSeqRun(data_id=template_name, experiment_id='20230929-testing', template_name=template_name)

    add_files(f'{cur_pseq.experiment_id}/{cur_pseq.data_id}/fastq', *fastq_path.iterdir())
    add_files(f'{cur_pseq.experiment_id}/{cur_pseq.data_id}/template', gb_path)

    engine = make_engine()
    SQLModel.metadata.create_all(engine)
    with Session(engine) as session:
        session.add(cur_pseq)
        session.commit()


if __name__ == '__main__':
    init()
    base_ip = '192.168.135.66'
    full_ip = f"http://{base_ip}:9000"
    lambda_uri = f"{full_ip}/2015-03-31/functions/function/invocations"
    print(lambda_uri)

    resp1 = requests.get(lambda_uri)


    print(r1_data)
    # full_pipeline(lambda_uri, fastq_path, gb_path)
