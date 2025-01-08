import re
from pathlib import Path
from typing import Optional, List, Tuple

import boto3
from Bio.Data import CodonTable
from gfapy.line.edge import Link
from nanoid import generate

from gfapy import Gfa
from networkx import Graph, MultiGraph
import networkx.algorithms as nxa
import pydna.all as pyd
from pydna.dseq import Dseq


def upload_files(*file_paths: Path) -> Optional[str]:
    """
    Uploads given files to a specific AWS S3 bucket.

    This function generates a unique nanoid as folder name, uploads both files
    into that folder on an AWS S3 bucket and returns the unique nanoid.

    Args:
        file_paths: The path of the second file.

    Returns:
        str: The unique nanoid used as folder name.
            Returns None if an error occurred during file upload.
    """
    # Create a session using your AWS credentials
    s3 = boto3.client('s3', region_name='us-east-1')

    # Generate a unique nanoid
    folder_name = generate()

    # Create file name with the newly created folder
    file_names_by_in_path = {str(p): f"{folder_name}/{p.name}" for p in file_paths}

    # Upload files to the S3 bucket
    try:
        for cur_path, cur_name in file_names_by_in_path.items():
            s3.upload_file(cur_path, 'foundry-plasmid-seq', cur_name)
    except Exception as e:
        print(e)
        return None

    # Return the unique nanoid
    return folder_name


def add_files(folder: str, *file_paths: Path, recursive=True) -> bool:
    """
    :param folder: A string indicating the name of the folder in the S3 bucket where the files will be uploaded.
    :param file_paths: Any number of `Path` objects representing the paths of the files to be uploaded.

    :return: A boolean value indicating whether the files were successfully uploaded to the S3 bucket.

    """
    # Create a session using your AWS credentials
    s3 = boto3.client('s3', region_name='us-east-1')

    # Create file name with the newly created folder
    file_names_by_in_path = {str(p): f"{folder}/{p.name}" for p in file_paths}

    # Upload files to the S3 bucket
    for cur_path, cur_name in file_names_by_in_path.items():
        try:
            s3.upload_file(cur_path, 'foundry-plasmid-seq', cur_name)
        except IsADirectoryError:
            if recursive:
                child_folder = Path(cur_path)
                add_files(cur_name, *child_folder.iterdir(), recursive=True)
    return True


def filter_invalid_paths(path: List[str]) -> bool:
    """
        Filters a list of string paths to determine if all required paths are present.

        Given a list of paths, this function ensures that for each node within the list,
        all its required child paths (both "L" and "R" suffixed versions) are also present
        in the list. The function works by identifying all potential parent nodes, deriving
        their left and right child paths, and confirming that these paths exist in the
        provided list of paths.

        Parameters
        ----------
        path : List[str]
            List of string paths where each string represents a path node.

        Returns
        -------
        bool
            A boolean value indicating whether all required paths are present.
    """
    nodes = set(n[:-1] for n in path)
    needed_nodes = {f"{n}L" for n in nodes}
    needed_nodes |= {f"{n}R" for n in nodes}
    return len(needed_nodes - set(path)) == 0


def plasmids_from_gfa(gfa_data: Gfa) -> list[tuple[Dseq, str, float]]:
    """
    :param gfa_data: Gfa object that represents the GFA data containing segment names, edges, and sequences.
    :return: List of pydna.Dseq objects representing the plasmids derived from the GFA data.

    This method takes in a GFA data object and constructs plasmids based on the segment names, edges, and sequences
    present in the GFA data. It creates a networkx graph and adds nodes and edges based on the segment names and
    edges in the GFA data. It also creates a dictionary to store the sequences of the segments.

    Next, it iterates over the edges in the GFA data and adds them to the graph. Then, it filters out invalid paths
    and finds all simple cycles in the graph.

    Finally, it constructs the full sequences of the plasmids by concatenating the sequences of the segments in each
    path. It returns a list of pydna.Dseq objects that represent the full sequences of the plasmids derived from the
    GFA data.
    """
    # First, check if this is even necessary
    if len(gfa_data.segment_names) == 1:
        segment_name = gfa_data.segment_names[0]
        segment = gfa_data.segment(segment_name)
        return [(pyd.Dseq(segment.sequence), f'{segment_name}+', 1.0)]


    graph = MultiGraph()
    # graph.add_nodes_from(g.segment_names)
    sequences = {}
    prevalence = {}
    for n in gfa_data.segment_names:
        graph.add_edge(f"{n}L", f"{n}R")
        cur_seq = pyd.Dseq(gfa_data.segment(n).sequence, linear=True)
        sequences[f"{n}L"] = cur_seq
        sequences[f"{n}R"] = cur_seq.reverse_complement()
        prevalence[f"{n}L"] = prevalence[f"{n}R"] = gfa_data.segment(n).get('dp')

    cur_edge: Link
    for cur_edge in gfa_data.edges:
        print(cur_edge.from_end, cur_edge.to_end)
        graph.add_edge(str(cur_edge.from_end), str(cur_edge.to_end))

    all_paths = list(filter(filter_invalid_paths, nxa.simple_cycles(graph)))

    full_sequences: List[pyd.Dseq] = []
    out_path_strings: List[str] = []
    min_prevalence: List[float] = []
    for path in all_paths:
        cur_seq = sequences[path[0]]
        for frag_start in path[2::2]:
            cur_seq = cur_seq + sequences[frag_start]
        full_sequences.append(cur_seq)
        out_path_strings.append(', '.join(path[0::2]).replace('L', '+').replace('R', '-'))
        min_prevalence.append(min(prevalence[f] for f in path))

    # Normalize min_prevalence
    if min_prevalence:
        max_prevalance = sum(min_prevalence)
        min_prevalence = [v/max_prevalance for v in min_prevalence]

    return list(zip(full_sequences, out_path_strings, min_prevalence))


def create_and_add_new_codon_table(table_name, base_table, modifications, start_codons=None, stop_codons=None):
    """
    Modify and create a new codon translation table based on an existing one.

    This function creates a new codon translation table by modifying an existing
    one. It applies specific modifications to codon-to-amino-acid translations and
    allows customization of start and stop codons. The new table is then added to
    the collection of available codon translation tables for unambiguous DNA.

    Args:
        table_name (str): The name of the new codon table.
        base_table (str): The name of the existing base codon table to modify.
        modifications (dict): A dictionary where keys are codons (str) and values
            are single-letter amino acid codes (str) representing modifications to
            the codon-to-amino-acid mapping.
        start_codons (list[str], optional): A list of codons to use as start codons
            for the new table. Defaults to the start codons of the base table.
        stop_codons (list[str], optional): A list of codons to use as stop codons
            for the new table. Defaults to the stop codons of the base table.
    """
    # Create a new table based on the given table
    new_table = dict(CodonTable.unambiguous_dna_by_name[base_table].forward_table)

    # Apply the modifications to the new table
    for codon, aa in modifications.items():
        new_table[codon] = aa

    # Use the base table's start/stop codons if new ones aren't specified
    if start_codons is None:
        start_codons = CodonTable.unambiguous_dna_by_name[base_table].start_codons

    if stop_codons is None:
        stop_codons = CodonTable.unambiguous_dna_by_name[base_table].stop_codons

    # Create a new CodonTable with the modified translations
    new_codon_table = CodonTable.CodonTable(forward_table=new_table,
                                            start_codons=start_codons,
                                            stop_codons=stop_codons)
    new_codon_table.names = [table_name]

    # Add the new table to the available unambiguous DNA tables
    CodonTable.unambiguous_dna_by_name[table_name] = new_codon_table


async def clean_genbank_file(genbank_path: Path, feature_blacklist: List[str] = None):
    feature_blacklist = feature_blacklist or ['Concatenated_se', 'Editing_History', 'extracted_regio', 'primer_bind',
                                              'restriction_sit', 'source', 'primer_bind_rev']
    genbank_fix_re1 = re.compile(r'^(\s{5}\w{15})\w(\d|complement|join)', re.MULTILINE)
    genbank_fix_re2 = re.compile(r'(\n\s{21}/[\w\s]+)\n\s{21}([\w\s]+)=')

    # Fix the record
    genbank_path.write_text(
        genbank_fix_re1.sub(r'\1 \2', genbank_path.read_text())
    )
    genbank_path.write_text(
        genbank_fix_re2.sub(r'\1_\2=', genbank_path.read_text())
    )

    record: pyd.Dseqrecord = pyd.read(str(genbank_path))
    record.features = [f for f in record.features if f.type not in feature_blacklist]
    record.write(str(genbank_path))



if __name__ == '__main__':
    plasmids_from_gfa(Gfa.from_file(r'C:\Users\RobertWarden-Rothman\PycharmProjects\denovo_plasmid_assembly\tests\GBFP-1384-0254\007_final_clean.gfa'))
