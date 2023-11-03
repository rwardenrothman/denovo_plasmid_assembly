from typing import Type, TypeVar

import requests
from sqlmodel import Session, SQLModel

from AppContainer.app.db_model import ExperimentStartModel, PlasmidSeqRunList, PlasmidSeqRun, AssemblyList

SM = TypeVar('SM', bound=SQLModel)


def parse_resp_body(response_class: Type[SM], response: requests.Response) -> SM:
    return response_class.parse_raw(response.json()['body'])


if __name__ == '__main__':
    base_ip = '192.168.135.66'
    full_ip = f"http://{base_ip}:9000"
    lambda_uri = f"{full_ip}/2015-03-31/functions/function/invocations"
    print(lambda_uri)

    start_model = ExperimentStartModel(name='1311_b')

    resp1 = requests.post(lambda_uri, json=start_model.lambda_json('start_experiment'))

    r1_data = parse_resp_body(PlasmidSeqRunList, resp1)
    for r in r1_data.runs:
        if r.data_id == 'Uox17_F9_MW':
            test_job = r
            break

    # for cur_step in ['setup', 'trim', 'assemble', 'transfer_annotations']:
    for cur_step in ['setup', 'trim', 'assemble', 'make_plasmids']:
        print(f"\nBeginning Step {cur_step}")
        resp = requests.post(lambda_uri, json=test_job.lambda_json(cur_step))
        print(resp.json())
        for c_type in [PlasmidSeqRun, AssemblyList, PlasmidSeqRunList]:
            try:
                test_job = parse_resp_body(c_type, resp)
                break
            except AttributeError:
                continue


    # full_pipeline(lambda_uri, fastq_path, gb_path)
