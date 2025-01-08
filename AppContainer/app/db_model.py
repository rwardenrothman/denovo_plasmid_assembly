import json
import os
from typing import Optional, List, Tuple, Union, TypeVar, Annotated

import boto3
from sqlmodel import SQLModel, Field, create_engine, Relationship
from pydantic import validator
from nanoid import generate as get_nanoid


class LambdaMixin(SQLModel):
    """
    Provides functionality for generating JSON payloads in a Lambda-compatible format.

    This mixin class is designed to complement SQLModel by adding a method to create
    JSON payloads in the structure expected by AWS Lambda. It constructs and returns
    a dictionary in the format of a mock HTTP POST request, enabling usage in
    Lambda-like environments or testing scenarios.
    """

    def lambda_json(self, function):
        """
        Generates a dictionary mimicking an AWS Lambda Proxy Integration event with a JSON body.

        This method creates a simulated Lambda event structure intended for testing
        or local development. The returned dictionary includes fields such as the body,
        resource path, HTTP method, and a minimal request context. The JSON body is
        produced by calling the instance's `json` method.

        Args:
            function (str): Name of the Lambda handler or function being invoked.

        Returns:
            dict: A dictionary representing an AWS Lambda Proxy Integration event with a
            JSON body.
        """
        return dict(
            body=self.json(),
            resource=f'/{function}',
            path=f'/{function}',
            httpMethod='POST',
            requestContext={}
        )


class PlasmidSeqRun(LambdaMixin, table=True):
    """
    Represents a plasmid sequencing run, storing relevant metadata and enabling
    data operations.

    This class is designed to store information about a plasmid sequencing run,
    including details about the experiment, associated assemblies, and potential
    errors encountered during processing. It provides methods for formatting paths
    and template genomes while enforcing validation on specific fields.

    Attributes:
        id: Optional unique identifier for the sequence run.
        data_id: Identifier for the associated data.
        experiment_id: Identifier for the associated experiment.
        last_step: Optional last recorded step of the sequencing process.
        template_name: Optional name of the sequencing template.
        template_length: Optional length of the sequencing template.
        basespace_href: Optional hyperlink to the Basespace resource.
        assembly_count: Count of assemblies associated with the sequencing run.
        assembly_coverage_mean: Optional mean coverage for associated assemblies.
        assembly_coverage_std: Optional standard deviation of coverage for assemblies.
        error_type: Optional type of error encountered in the sequencing run.
        error_message: Optional error message for the encountered issue.
        error_path: Optional path to the resource that caused the error.
        assemblies: List of related plasmid sequence assemblies.
    """
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
        """
        Generates a formatted file path based on experiment ID and data ID.

        This method constructs a specific file path by combining the experiment ID
        and data ID as base components, and optionally appending additional folder
        names provided as arguments. Spaces in the experiment ID are replaced with
        underscores to ensure valid path formatting.

        Parameters:
            folders: str
                Optional folder names to be appended to the base path. Multiple
                folder names can be provided as separate arguments.

        Returns:
            str: A formatted file path string composed of the base experiment and
            data ID, optionally followed by the appended folder names.
        """
        path = f"{self.experiment_id.replace(' ', '_')}/{self.data_id}"
        for folder in folders:
            path += '/' + folder
        return path

    @property
    def template_gb(self) -> str:
        """
            Generates a formatted string representing the template name with a '.gb' file
            extension.

            This property provides a convenient way to access the full template name
            with the '.gb' suffix appended.

            Return:
                str: The formatted template name with a '.gb' file extension.
        """
        return f"{self.template_name}.gb"

    @validator('template_name', 'data_id', 'experiment_id')
    @classmethod
    def strip_template_name(cls, v: str) -> str:
        """
        Validates and processes the given string attributes for the specified fields.

        This method ensures that the provided string values for the given fields
        are appropriately stripped of leading and trailing whitespace. This can
        be used as a validation method in data models or for preprocessing
        input data fields.

        Parameters:
            v: The string value to be processed.

        Returns:
            The processed string value with leading and trailing whitespace removed.
        """
        return v.strip()

class PlasmidSeqAssembly(LambdaMixin, table=True):
    """
    Represents a plasmid sequence assembly.

    This class models a plasmid sequence assembly including its associated attributes,
    relationships, and a method to get its GenBank formatted file name. It is a table-enabled
    class with relationships to other relevant entities, making it suitable for use in a database
    context.

    Attributes:
        id (Optional[int]): Unique identifier for the plasmid sequence assembly. Defaults to None.
        run_id (int): Foreign key that links to a specific plasmid sequence run.
        assembly_name (str): Name of the assembly.
        contig_path (str): Path to the contig file for the assembly.
        length (int): Length of the assembly in base pairs.
        min_prevalence (float): Minimum prevalence threshold for identifying features.

    Relationships:
        sample (PlasmidSeqRun): Associated plasmid sequencing run connected through a back_populates
            relationship with the 'assemblies' attribute.
        features (List[AssemblyFeature]): List of assembly features identified in the plasmid
            sequence, linked via a back_populates relationship with the 'assembly' attribute.
        polymorphisms (List[AssemblyPolymorphism]): List of assembly polymorphisms linked to the
            plasmid sequence, connected via a back_populates relationship with the 'assembly'
            attribute.

    Properties:
        assembly_gb (str): Returns the GenBank formatted file name for the assembly. The name is
            derived by appending a ".gb" extension to the assembly_name attribute.
    """
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
        """The assembly GenBank formatted file name."""
        return f"{self.assembly_name}.gb"


class AssemblyList(LambdaMixin):
    """
    Represents a collection of plasmid sequence assemblies.

    This class is designed to encapsulate and manage a list of PlasmidSeqAssembly
    objects. It may provide utilities for operating on these assemblies, enabling
    structured handling and manipulation of the data related to plasmid sequencing
    processes.

    Attributes:
        assemblies (List[PlasmidSeqAssembly]): List of plasmid sequence assemblies
        managed within this collection.

    :meta private:
    """
    assemblies: List[PlasmidSeqAssembly]


class FeaturePolymorphism(SQLModel, table=True):
    """
    Represents a many-to-many relationship between features and polymorphisms
    within a database for modeling genetic associations.

    This class is designed as a table representation to establish a relational link
    between `Feature` and `Polymorphism` entities while also supporting SQLModel
    functionalities. It includes primary key constraints to ensure the uniqueness
    of the connection.

    Attributes:
    feature_id: Optional int
        The identifier linking to a feature in the `assemblyfeature` table.
    polymorphism_id: Optional int
        The identifier linking to a polymorphism in the `assemblypolymorphism` table.
    """
    feature_id: Optional[int] = Field(default=None, foreign_key='assemblyfeature.id', primary_key=True)
    polymorphism_id: Optional[int] = Field(default=None, foreign_key='assemblypolymorphism.id', primary_key=True)


class AssemblyFeature(SQLModel, table=True):
    """
    Represents a feature within a plasmid sequence assembly.

    This class is used to define the details of a specific feature associated
    with a plasmid sequence assembly. It includes information about the
    feature type, its reference name, its name in the assembly, and any
    potential changes made to it such as deletions or frameshift residues.
    Relationships to the parent assembly and associated polymorphisms are
    also included for data integrity and query purposes.

    Attributes:
        id (Optional[int]): The unique identifier for the AssemblyFeature.
            This is the primary key in the database table.
        assembly_id (int): The foreign key referencing the 'plasmidseqassembly.id'.
        feature_type (str): The type of the feature (e.g., gene, CDS, etc.).
        wt_feature_name (str): The name of the feature in the wild-type sequence.
        assembly_feature_name (str): The name of the feature as identified in
            the assembly.
        deleted (bool): Indicates whether the feature was deleted. Defaults to False.
        frameshift_residue (Optional[int]): The residue position associated
            with a frameshift mutation, if applicable.
        assembly (PlasmidSeqAssembly): The parent plasmid sequence assembly object
            associated with the feature.
        polymorphisms (List[AssemblyPolymorphism]): A list of polymorphisms linked
            to the feature through a many-to-many relationship using the
            FeaturePolymorphism link model.
    """
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
    """
    A class representing assembly polymorphisms and their relationships.

    This class defines the properties and relationships of polymorphisms in an assembly. It is
    a SQLModel table that includes fields for referencing assemblies, polymorphism types,
    coordinate ranges, and potential effects on coding sequences (CDS). The relationships to
    assemblies and associated features are also defined, enabling efficient querying and
    management of polymorphism data.

    Attributes:
        id (Optional[int]): Unique identifier for the polymorphism entry; primary key.
        assembly_id (int): Foreign key referencing the associated assembly.
        poly_type (str): Type of polymorphism (e.g., "SNP", "insertion", "deletion").
        wt_nt_start (Optional[int]): Start coordinate of the wild-type nucleotide sequence.
        wt_nt_end (Optional[int]): End coordinate of the wild-type nucleotide sequence.
        assembly_nt_start (Optional[int]): Start coordinate of the polymorphism in the assembly.
        assembly_nt_end (Optional[int]): End coordinate of the polymorphism in the assembly.
        cds_effect (Optional[str]): Description of the polymorphism's effect on coding sequences.

    Relationships:
        assembly (PlasmidSeqAssembly): Relationship linking to the associated
            plasmid sequence assembly.
        features (List[AssemblyFeature]): Relationship linking to associated assembly features,
            managed with the FeaturePolymorphism link model.
    """
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
    """Represents mapping status information in a structured format.

    This class defines a model for mapping status information used within applications
    that leverage structured data for processing and integration. It inherits from the
    SQLModel, making it suitable for use with SQL-based databases. It contains attributes
    describing the status of a specific process or operation, alongside its associated
    error and metadata information.

    Attributes:
        Error: Optional description of the error as a string.
        Cause: Optional detailed cause of the error as a string.
        RunID: Optional identifier for the associated run as an integer.
    """
    Error: Optional[str]
    Cause: Optional[str]
    RunID: Optional[int]


RunListTypes = Union[PlasmidSeqRun, MapStatusInfo, "PlasmidSeqRunList"]
"""Type alias for run list types."""
RLT = TypeVar("RLT", bound=RunListTypes)
"""Type variable for run list types."""


class PlasmidSeqRunList(LambdaMixin):
    """
    Represents a list of plasmid sequencing runs.

    The PlasmidSeqRunList class is used to organize and manipulate a collection of
    sequencing runs, which could be of various types. It provides functionality
    to filter and retrieve runs based on their type, traversing nested instances
    of PlasmidSeqRunList if necessary. This class extends the functionality of the
    LambdaMixin class for additional utility.

    Attributes:
        runs (List[RunListTypes]): A list containing sequencing runs, which may
            include nested PlasmidSeqRunList instances or other run types.
    """
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
    """
    Represents the data model for starting an experiment.

    This class is designed to encapsulate the essential data
    and structure required to initiate an experiment. It serves
    as a container for relevant attributes and integrates with
    LambdaMixin for additional functionality.

    Attributes:
        name (str): The name of the experiment being initiated.
    """
    name: str


def make_engine():
    print('Entered make engine')
    engine = create_engine(os.environ['DB_URI'])
    print('Engine Created')
    return engine


if __name__ == '__main__':
    SQLModel.metadata.create_all(make_engine())
