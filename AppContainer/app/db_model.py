import json
import os
from typing import Optional, List, Tuple, Union, TypeVar, Annotated

import boto3
from sqlmodel import SQLModel, Field, create_engine, Relationship
from pydantic import validator
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
    template_length: Optional[int]
    basespace_href: Optional[str]
    assembly_count: int = Field(default=0)
    assembly_coverage_mean: Optional[float]
    assembly_coverage_std: Optional[float]
    error_type: Optional[str]
    error_message: Optional[str]
    error_path: Optional[str]

    assemblies: List["PlasmidSeqAssembly"] = Relationship(back_populates='sample')

    def data_path(self, *folders: str) -> str:
        path = f"{self.experiment_id.replace(' ', '_')}/{self.data_id}"
        for folder in folders:
            path += '/' + folder
        return path

    @property
    def template_gb(self) -> str:
        return f"{self.template_name}.gb"

    @validator('template_name', 'data_id', 'experiment_id')
    @classmethod
    def strip_template_name(cls, v: str) -> str:
        return v.strip()

class PlasmidSeqAssembly(LambdaMixin, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    run_id: int = Field(foreign_key='plasmidseqrun.id')
    assembly_name: str
    contig_path: str
    length: int
    min_prevalence: float

    sample: PlasmidSeqRun = Relationship(back_populates='assemblies')
    features: List["AssemblyFeature"] = Relationship(back_populates='assembly')
    polymorphisms: List["AssemblyPolymorphism"] = Relationship(back_populates='assembly')

    @property
    def assembly_gb(self) -> str:
        return f"{self.assembly_name}.gb"


class AssemblyList(LambdaMixin):
    assemblies: List[PlasmidSeqAssembly]


class FeaturePolymorphism(SQLModel, table=True):
    feature_id: Optional[int] = Field(default=None, foreign_key='assemblyfeature.id', primary_key=True)
    polymorphism_id: Optional[int] = Field(default=None, foreign_key='assemblypolymorphism.id', primary_key=True)


class AssemblyFeature(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    assembly_id: int = Field(foreign_key='plasmidseqassembly.id')
    feature_type: str
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
    poly_type: str
    wt_nt_start: Optional[int]
    wt_nt_end: Optional[int]
    assembly_nt_start: Optional[int]
    assembly_nt_end: Optional[int]
    cds_effect: Optional[str]

    assembly: PlasmidSeqAssembly = Relationship(back_populates='polymorphisms')
    features: List["AssemblyFeature"] = Relationship(back_populates='polymorphisms', link_model=FeaturePolymorphism)


class MapStatusInfo(SQLModel):
    Error: Optional[str]
    Cause: Optional[str]
    RunID: Optional[int]


RunListTypes = Union[PlasmidSeqRun, MapStatusInfo, "PlasmidSeqRunList"]
RLT = TypeVar("RLT", bound=RunListTypes)


class PlasmidSeqRunList(LambdaMixin):
    runs: List[RunListTypes]

    def runs_of_type(self, run_type: type[RLT]) -> list[RLT]:
        out_list = []
        for r in self.runs:
            if isinstance(r, run_type):
                out_list.append(r)
            elif isinstance(r, PlasmidSeqRunList):
                out_list.extend(r.runs_of_type(run_type))
        return out_list


class ExperimentStartModel(LambdaMixin):
    name: str


def make_engine():
    print('Entered make engine')
    engine = create_engine(os.environ['DB_URI'])
    print('Engine Created')
    return engine


if __name__ == '__main__':
    SQLModel.metadata.create_all(make_engine())
