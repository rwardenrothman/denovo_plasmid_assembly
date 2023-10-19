import json
import os
from typing import Optional, List, Tuple

import boto3
from sqlmodel import SQLModel, Field, create_engine, Relationship
from nanoid import generate as get_nanoid


class LambdaMixin(SQLModel):

    def lambda_json(self, function):
        return dict(
            body=self.json(),
            resource=f'/{function}',
            path=f'/{function}',
            httpMethod='POST',
            requestContext={}
        )


class PlasmidSeqRun(LambdaMixin, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    data_id: str
    experiment_id: str
    last_step: Optional[str]
    template_name: Optional[str]
    basespace_href: Optional[str]
    error_type: Optional[str]
    error_message: Optional[str]
    error_path: Optional[str]

    def data_path(self, folder: str = None) -> str:
        path = f"{self.experiment_id}/{self.data_id}"
        if folder:
            path += '/' + folder
        return path


class PlasmidSeqRunList(LambdaMixin):
    runs: List[PlasmidSeqRun]

class ExperimentStartModel(LambdaMixin):
    name: str


def make_engine():
    print('Entered make engine')
    engine = create_engine(os.environ['DB_URI'])
    print('Engine Created')
    return engine