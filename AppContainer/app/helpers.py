from itertools import product
from pathlib import Path
from typing import Optional, List

import boto3
from gfapy.line.edge import Link
from nanoid import generate

from gfapy import Gfa
from networkx import MultiDiGraph
import networkx.algorithms as nxa
import pydna.all as pyd


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


def add_files(folder: str, *file_paths: Path) -> bool:
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
    try:
        for cur_path, cur_name in file_names_by_in_path.items():
            s3.upload_file(cur_path, 'foundry-plasmid-seq', cur_name)
    except Exception as e:
        print(e)
        return False

    return True


def plasmids_from_gfa(gfa_data: Gfa) -> List[pyd.Dseq]:
    graph = MultiDiGraph()
    # graph.add_nodes_from(g.segment_names)
    sequences = {}
    for n in gfa_data.segment_names:
        graph.add_edge(f"{n}L", f"{n}R")
        cur_seq = pyd.Dseq(gfa_data.segment(n).sequence, linear=True)
        sequences[f"{n}L"] = cur_seq
        sequences[f"{n}R"] = cur_seq.reverse_complement()

    cur_edge: Link
    for cur_edge in gfa_data.edges:
        print(cur_edge.from_end, cur_edge.to_end)
        graph.add_edge(str(cur_edge.from_end), str(cur_edge.to_end))

    roots = [n for n, d in graph.in_degree if d == 0]
    leaves = [n for n, d in graph.out_degree if d == 0]

    all_paths = []
    for cur_root, cur_leaf in product(roots, leaves):
        for cur_path in list(nxa.all_simple_paths(graph, cur_root, cur_leaf)):
            if str(cur_path[0])[:-1] != str(cur_path[1])[:-1]:
                cur_path = [f"{str(cur_path[0])[:-1]}R"] + cur_path
            if str(cur_path[-1])[:-1] != str(cur_path[-2])[:-1]:
                cur_path = cur_path + [f"{str(cur_path[-1])[:-1]}L"]
            print(cur_path)
            all_paths.append(cur_path)

    full_sequences = []
    for path in all_paths:
        # print(path)
        print(path[::2])

        cur_seq = sequences[path[0]]
        for frag_start in path[2::2]:
            cur_seq = cur_seq + sequences[frag_start]
        full_sequences.append(cur_seq)

    return full_sequences
