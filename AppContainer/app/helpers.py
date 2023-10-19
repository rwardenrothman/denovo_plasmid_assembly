from pathlib import Path
from typing import Optional

import boto3
from nanoid import generate


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
