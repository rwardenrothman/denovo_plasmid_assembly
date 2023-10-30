import json
import os
from typing import Optional, List, Tuple, Union

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
    assembly_count: int = Field(default=0)
    error_type: Optional[str]
    error_message: Optional[str]
    error_path: Optional[str]

    assemblies: List["PlasmidSeqAssembly"] = Relationship(back_populates='sample')

    def data_path(self, folder: str = None) -> str:
        path = f"{self.experiment_id}/{self.data_id}"
        if folder:
            path += '/' + folder
        return path


class PlasmidSeqAssembly(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    run_id: int = Field(foreign_key='plasmidseqrun.id')
    assembly_name: str
    contig_path: str
    length: int

    sample: PlasmidSeqRun = Relationship(back_populates='assemblies')
    features: List["AssemblyFeature"] = Relationship(back_populates='assembly')
    polymorphisms: List["AssemblyPolymorphism"] = Relationship(back_populates='assembly')


class FeaturePolymorphism(SQLModel, table=True):
    feature_id: Optional[int] = Field(default=None, foreign_key='assemblyfeature.id', primary_key=True)
    polymorphism_id: Optional[int] = Field(default=None, foreign_key='assemblypolymorphism.id', primary_key=True)


class AssemblyFeature(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    assembly_id: int = Field(foreign_key='plasmidseqassembly.id')
    wt_feature_name: str
    assembly_feature_name: str
    deleted: bool = Field(default=False)
    frameshift_residue: Optional[int]

    assembly: PlasmidSeqAssembly = Relationship(back_populates='features')
    polymorphisms: List["AssemblyPolymorphism"] = Relationship(back_populates='features',
                                                               link_model=FeaturePolymorphism)


class AssemblyPolymorphism(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    assembly_id: int = Field(foreign_key='plasmidseqassembly.id')
    wt_nt_start: int
    wt_nt_end: int
    assembly_nt_start: int
    assembly_nt_end: int
    cds_effect: Optional[str]

    assembly: PlasmidSeqAssembly = Relationship(back_populates='features')
    features: List["AssemblyFeature"] = Relationship(back_populates='polymorphisms', link_model=FeaturePolymorphism)


class MapStatusInfo(SQLModel):
    Status: Optional[str]
    Error: Optional[str]
    Cause: Optional[str]


class PlasmidSeqRunList(LambdaMixin):
    runs: List[Union[PlasmidSeqRun, MapStatusInfo]]


class ExperimentStartModel(LambdaMixin):
    name: str


def make_engine():
    print('Entered make engine')
    engine = create_engine(os.environ['DB_URI'])
    print('Engine Created')
    return engine


if __name__ == '__main__':
    SQLModel.metadata.create_all(make_engine())
