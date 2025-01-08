import re
from pathlib import Path
from itertools import repeat, count, starmap, cycle
from collections import deque
from typing import Iterable, TypeVar, Optional, Tuple, Dict, Union

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.SeqFeature import FeatureLocation
from pandas import DataFrame
from pydna.dseqrecord import Dseqrecord
from pydna.seqfeature import SeqFeature
from sqlmodel import Session
from Bio.Align import PairwiseAligner, PairwiseAlignment, PairwiseAlignments
import pydna.all as pyd


if 'AppContainer' in __file__:
    from AppContainer.app.helpers import create_and_add_new_codon_table
    from AppContainer.app.db_model import PlasmidSeqAssembly, AssemblyFeature, AssemblyPolymorphism
else:
    from helpers import create_and_add_new_codon_table
    from db_model import PlasmidSeqAssembly, AssemblyFeature, AssemblyPolymorphism

create_and_add_new_codon_table('RC63', 'Standard', {'TAG': 'U'}, stop_codons=['TGA', 'TAA'])


T = TypeVar('T')


def threepeat(target: Iterable[T], max_iters: int) -> Iterable[T]:
    """
    Provides a generator function that replicates each item in the input iterable three times and stops after a maximum number of items
    have been yielded, based on the specified limit. The function supports input validation for the 'max' parameter.

    Parameters:
        target: Iterable[T]
            An iterable containing the items to be repeated.
        max_iters: int
            The maximum number of items to yield. Must be a positive integer.

    Returns:
        Iterable[T]:
            A generator that yields items from the input iterable, each repeated three times, and stops
            after yielding 'max' items.

    Raises:
        None
    """
    if not isinstance(max_iters, int) or max_iters < 1:
        return []
    cur_count = 0
    try:
        for i in target:
            for j in repeat(i, 3):
                yield j
                cur_count += 1
                if cur_count == max_iters:
                    raise StopIteration()
    except StopIteration:
        pass


def triplets(target: Iterable[T]) -> Iterable[T]:
    """
    Generates iterable triplets from a given input iterable.

    This function takes an input iterable and processes it to produce triplets,
    with each triplet containing up to three elements from the original iterable.
    The generator continues until all elements from the input iterable have been
    processed.

    Parameters:
    target: Iterable
        The input iterable to be split into triplets.

    Yields:
    Iterable
        A triplet (up to three elements) chunk of the input iterable, returned
        as another iterable.

    """
    while target:
        triplet, target = target[:3], target[3:]
        yield triplet


def find_frameshift(seq: str) -> Optional[int]:
    """
    Find the position of a frameshift in a DNA sequence.

    This function identifies the first codon containing a gap symbol ('-'),
    indicating a frameshift in the DNA sequence. It returns the 1-based
    index of the frameshift location or None if no frameshift is found.

    Parameters:
    seq: str
        The input DNA sequence to analyze. It is expected to consist of
        nucleotide characters ('A', 'T', 'C', 'G') and may contain '-'
        for gaps.

    Returns:
    Optional[int]
        The 1-based index of the codon containing a frameshift. Returns
        None if no frameshift is identified.
    """
    for i, cur_codon in enumerate(triplets(seq)):
        if '-' in cur_codon:
            return i + 1
    return None


def get_best_alignment(alignments: PairwiseAlignments) -> PairwiseAlignment:
    """
        Selects the best alignment from a collection of alignments.

        This function takes a collection of pairwise alignments and determines the
        best alignment based on sorting. It handles memory or overflow errors by
        falling back to returning the next alignment in the iterator when sorting
        fails.

        Args:
            alignments: A PairwiseAlignments object, an iterable of alignments.

        Returns:
            A PairwiseAlignment object representing the best alignment determined
            either by sorting or by fallback to the first available alignment in the
            iterable.

        Raises:
            MemoryError: If sorting fails to process due to memory constraints.
            OverflowError: If sorting encounters an arithmetic overflow.
    """
    try:
        return sorted(alignments)[0]
    except (MemoryError, OverflowError):
        return next(alignments)


def align_sequences(template_path: Path, assembly_path: Path,
                    codon_table='RC63') -> tuple[DataFrame, Dseqrecord, dict[int, SeqFeature]]:
    """
    Aligns a template DNA sequence to an assembled DNA sequence while considering circular templates,
    detecting polymorphisms, and associating sequence features for downstream analysis.

    The function performs the following key steps:
    - Reads in the template and assembly DNA sequences.
    - Aligns the sequences using pairwise alignment, accounting for sequence circularity.
    - Determines and applies sequence shifts to optimize alignment.
    - Maps base pair position correspondence between the template and assembly sequences.
    - Identifies and annotates polymorphisms such as SNPs, insertions, and deletions.
    - Matches features from the template sequence to aligned regions, including CDS-specific details.
    - Resolves nucleotide-to-amino acid mappings for CDS features with consideration of codon positions.

    Parameters:
        template_path (Path): The file path to the template DNA sequence in GenBank format.
        assembly_path (Path): The file path to the assembled DNA sequence in GenBank format.
        codon_table (str): The name of the codon table used for translating sequences during analysis.
                           Defaults to 'RC63'.

    Returns:
        tuple[DataFrame, Dseqrecord, dict[int, SeqFeature]]: A tuple consisting of:
          - A DataFrame mapping the correspondence of base pair positions, nucleotides, and polymorphisms across
            the template and assembly sequences, with feature-related annotations.
          - The assembly sequence record after alignment and adjustment for proper translation.
          - A dictionary mapping feature IDs in the template sequence to the extracted SeqFeature objects,
            excluding those with missing or invalid labels.

    Raises:
        AssemblyTooShortError: Raised when the assembly sequence is significantly shorter than half the length
                               of the template sequence, preventing meaningful alignment.
    """
    codon_table_obj = CodonTable.unambiguous_dna_by_name[codon_table]

    # Read in files
    template_record: pyd.Dseqrecord = pyd.read(str(template_path))
    assembly_record: pyd.Dseqrecord = pyd.read(str(assembly_path))
    if len(assembly_record) < len(template_record)/2:
        raise AssemblyTooShortError(f"{assembly_record.name} is too short to properly align")

    # Align circular templates
    aligner = PairwiseAligner()
    aligner.open_gap_score = -.05 * len(template_record)
    aligner.query_end_gap_score = 0.
    aligner.target_end_gap_score = 0.
    alignments = aligner.align(str(template_record.seq), str(assembly_record.seq) * 2)
    alignment = get_best_alignment(alignments)

    # Determine sequence shift and re-align
    shift = dict(alignment.indices.T)[0]
    assembly_record = assembly_record.shifted(shift)

    final_aligner = PairwiseAligner()
    final_aligner.open_gap_score = -5
    final_alignments = final_aligner.align(str(template_record.seq), str(assembly_record.seq))
    best_alignment = get_best_alignment(final_alignments)

    bp_map = pd.DataFrame(best_alignment.indices.T, columns=['template_pos', 'assembly_pos'])

    # Add nucleotides
    # bp_map['template_nt'] = bp_map['assembly_nt'] = None
    bp_map = bp_map.merge(pd.DataFrame(enumerate(template_record.seq), columns=['template_pos', 'template_nt']),
                          'left', 'template_pos')
    bp_map = bp_map.merge(pd.DataFrame(enumerate(assembly_record.seq), columns=['assembly_pos', 'assembly_nt']),
                          'left', 'assembly_pos')

    # Find polymorphisms
    bp_map.loc[bp_map.template_nt != bp_map.assembly_nt, 'poly_type'] = 'SNP'
    bp_map.loc[bp_map.template_nt.isna(), 'poly_type'] = 'Insertion'
    bp_map.loc[bp_map.assembly_nt.isna(), 'poly_type'] = 'Deletion'
    bp_map.replace('nan', np.nan, inplace=True)

    # Give the polymorphisms an id number
    polys_only = bp_map[bp_map['poly_type'].notna()].copy()
    diffs0 = polys_only.index.to_series().diff()
    diffs1 = abs(polys_only['template_pos'].diff().replace(0, 1))
    diffs2 = abs(polys_only['assembly_pos'].diff().replace(0, 1))
    diffs = (diffs0 != 1) | (diffs1 != 1) | (diffs2 != 1)
    poly_ids = (diffs).cumsum()
    bp_map.loc[polys_only.index, 'poly_id'] = poly_ids
    bp_map['template_aa'] = None
    bp_map['assembly_aa'] = None


    # Add features to the table
    features_by_id: dict[int, SeqFeature] = dict(enumerate(template_record.sorted_features()))
    unknown_feature_ids = []

    for feature_id, cur_feature in features_by_id.items():
        cur_feature_name = cur_feature.qualifiers.get('label', [None])[0]
        if not cur_feature_name:
            unknown_feature_ids.append(feature_id)
        feature_indices = []
        for cur_loc in cur_feature.location.parts:
            start_index = bp_map[bp_map['template_pos'] == cur_loc.start].index.min()
            end_index = bp_map[bp_map['template_pos'] == cur_loc.end].index.max()
            if pd.isna(end_index):
                end_index = bp_map.index.max()+1
            feature_indices.extend(range(start_index, end_index))

        template_indices: list[int] = list(bp_map.loc[feature_indices]['template_nt'].dropna().index.values)
        assembly_indices: list[int] = list(bp_map.loc[feature_indices]['assembly_nt'].dropna().index.values)

        if cur_feature.location.strand == -1:
            feature_indices.reverse()
            template_indices.reverse()
            assembly_indices.reverse()

        bp_map.loc[feature_indices, 'feature_id'] = feature_id
        bp_map.loc[feature_indices, 'feature_name'] = cur_feature_name
        bp_map.loc[feature_indices, 'feature_type'] = cur_feature.type
        bp_map.loc[feature_indices, 'feature_strand'] = cur_feature.strand

        template_len = len(template_indices)
        assembly_len = len(assembly_indices)
        if cur_feature.type == 'CDS':
            bp_map.loc[template_indices, 'template_res_id'] = list(threepeat(count(1), template_len))

            wt_seq = pyd.Dseq(''.join(bp_map.loc[template_indices, 'template_nt']))
            if cur_feature.strand == -1:
                wt_seq = wt_seq.complement()
            wt_aa = wt_seq.translate(codon_table_obj) + ('X' * template_len)
            bp_map.loc[template_indices, 'template_aa'] = list(threepeat(wt_aa, template_len))

            t_codon_pos = [1, 2, 3] * template_len
            bp_map.loc[template_indices, 'template_codon_pos'] = t_codon_pos[:template_len]

            bp_map.loc[assembly_indices, 'assembly_res_id'] = list(threepeat(count(1), assembly_len))

            asy_seq = pyd.Dseq(''.join(bp_map.loc[assembly_indices, 'assembly_nt'].fillna('-')))
            if cur_feature.strand == -1:
                asy_seq = asy_seq.complement()
            asy_aa = asy_seq.translate(codon_table_obj) + ('X' * assembly_len)
            bp_map.loc[assembly_indices, 'assembly_aa'] = list(threepeat(asy_aa, assembly_len))

            a_codon_pos = [1, 2, 3] * assembly_len
            bp_map.loc[assembly_indices, 'assembly_codon_pos'] = a_codon_pos[:assembly_len]

            bp_map.loc[feature_indices, 'reading_frame'] = (bp_map.loc[feature_indices, 'assembly_codon_pos'] -
                                                            bp_map.loc[feature_indices, 'template_codon_pos']) % 3

            if cur_feature.strand == -1:
                bp_map.loc[feature_indices, 'reading_frame'] = bp_map.loc[feature_indices, 'reading_frame'].bfill()
                bp_map.loc[feature_indices, 'template_res_id'] = bp_map.loc[feature_indices, 'template_res_id'].bfill()
                bp_map.loc[feature_indices, 'assembly_res_id'] = bp_map.loc[feature_indices, 'assembly_res_id'].bfill()
            else:
                bp_map.loc[feature_indices, 'reading_frame'] = bp_map.loc[feature_indices, 'reading_frame'].ffill()
                bp_map.loc[feature_indices, 'template_res_id'] = bp_map.loc[feature_indices, 'template_res_id'].ffill()
                bp_map.loc[feature_indices, 'assembly_res_id'] = bp_map.loc[feature_indices, 'assembly_res_id'].ffill()

        else:
            bp_map.loc[template_indices, 'template_res_id'] = list(range(1, template_len + 1))
            bp_map.loc[assembly_indices, 'assembly_res_id'] = list(range(1, assembly_len + 1))
            bp_map.loc[template_indices, 'template_aa'] = np.nan
            bp_map.loc[assembly_indices, 'assembly_aa'] = np.nan
            bp_map.loc[template_indices, 'template_codon_pos'] = np.nan
            bp_map.loc[assembly_indices, 'assembly_codon_pos'] = np.nan
            bp_map.loc[feature_indices, 'reading_frame'] = np.nan

    for fid in unknown_feature_ids:
        del features_by_id[fid]

    return bp_map, assembly_record, features_by_id


def seq(x):
    """
        Concatenates the non-NA/null values of the given input into a single string.

        This function takes an input, drops all NA/null elements, and joins the
        remaining elements into a single string.

        Parameters
        ----------
        x: pandas.Series
            Input series containing elements to be processed.

        Returns
        -------
        str
            A string resulting from the concatenation of non-NA/null elements
            of the input series.
    """
    return ''.join(x.dropna())


RF_DELIMITER = '!rf'


def find_effect(row_index, row_data):
    """
    Finds and constructs a textual representation of the effect of a given row's genetic
    mutation or polymorphism. The function analyzes mutation types such as SNP (Single
    Nucleotide Polymorphism), Deletion, and Insertion to derive a concise description of
    the mutation effect based on the provided row data. In case of invalid data, proper
    default behaviors are applied, and None is returned if the effect cannot be determined.

    Args:
        row_index (int): The index of the row being processed.
        row_data (pd.Series): A series containing mutation data, which includes fields
            such as poly_type, reading_frame, template_res_id, assembly_nt, template_nt,
            template_aa, assembly_aa, and feature_strand.

    Returns:
        Optional[str]: A string describing the mutation effect if it can be computed
        from the supplied data. Returns None if the effect cannot be determined.
    """
    effect_text = None
    poly_type = row_data[('poly_type', 'first')]
    reading_frame = row_data[('reading_frame', 'first')]
    reading_frame = None if pd.isna(reading_frame) else int(reading_frame)

    try:
        res_min, res_max = int(row_data[('template_res_id', 'min')]), int(row_data[('template_res_id', 'max')])
        t_nt, a_nt = row_data[('template_nt', 'seq')], row_data[('assembly_nt', 'seq')]
        t_aa, a_aa = row_data[('template_aa', 'seq')][::3], row_data[('assembly_aa', 'seq')][::3],
        if row_data[('feature_strand', 'min')] == -1:
            t_aa, a_aa = t_aa[::-1], a_aa[::-1]

        if poly_type == 'SNP':
            if t_aa == a_aa:
                effect_text = 'None'
            else:
                effect_text = f'{t_aa}{res_min:d}{a_aa}'

        elif poly_type == 'Deletion':
            effect_text = f"del{res_min:d}"
            if res_max != res_min:
                effect_text += f".{res_max:d}"
            if reading_frame is not None:
                effect_text += f"{RF_DELIMITER}{reading_frame:d}"

        elif poly_type == 'Insertion':
            effect_text = f"ins{res_min:d}"
            effect_text += a_aa or a_nt
            if reading_frame is not None:
                effect_text += f"{RF_DELIMITER}{reading_frame:d}"

    except ValueError as e:
        pass
    return effect_text


def get_polymorphism_features(bp_map: pd.DataFrame) -> pd.DataFrame:
    """
    Extract polymorphism features from a genome base-pair mapping.

    This function processes a DataFrame containing genome base-pair mappings,
    identifying positions where polymorphisms occur, calculates relevant
    differences for those positions, and organizes them into a summarized
    DataFrame for further analysis. It fills missing positions forward based
    on the `template_pos` and `assembly_pos` columns and extracts details
    about the polymorphism effects using the provided utility function.

    Parameters
    ----------
    bp_map : pd.DataFrame
        Input DataFrame that contains genome base-pair mappings. It is
        expected to have required columns, such as 'template_pos',
        'assembly_pos', 'poly_id', 'poly_type', and other necessary fields
        for identifying and describing polymorphisms.

    Returns
    -------
    pd.DataFrame
        A DataFrame summarizing polymorphism features. The output includes
        details about polymorphism type, genomic locations, sequence changes,
        and potential functional effects. Grouped and aggregated polymorphic
        information is provided in a multi-index format.

    Raises
    ------
    ValueError
        Raised if required columns are missing from the input DataFrame.
    """
    bp_map = bp_map.copy()
    bp_map['template_pos'] = bp_map['template_pos'].replace(-1, None).ffill()
    bp_map['assembly_pos'] = bp_map['assembly_pos'].replace(-1, None).ffill()

    diffs = bp_map[bp_map['poly_id'].notnull()].copy()
    polys = diffs.groupby('poly_id').agg({'poly_type': ['first'],
                                          'feature_strand': ['min'],
                                          'template_pos': ['min', 'max'],
                                          'assembly_pos': ['min', 'max'],
                                          'template_nt': [seq],
                                          'assembly_nt': [seq],
                                          'template_aa': [seq],
                                          'assembly_aa': [seq],
                                          'template_res_id': ['min', 'max'],
                                          'reading_frame': ['first']
                                          })
    polys[('Effect', 'text')] = list(starmap(find_effect, polys.iterrows()))
    return polys


def list_unique(x: pd.Series) -> str:
    """
    Returns a string containing unique non-null values in a pandas Series, joined by commas.

    This function processes a given pandas Series, removes null values, extracts
    the unique values, converts them to strings, and concatenates them into a
    single string separated by commas.

    Arguments:
        x (pd.Series): The pandas Series from which unique non-null values
        are extracted.

    Returns:
        str: A comma-separated string of unique non-null values from the
        input Series.
    """
    return ', '.join(x.dropna().unique().astype(str))


def get_features_with_polymorphisms(bp_map: pd.DataFrame) -> pd.DataFrame:
    """
        Extracts features containing sequence polymorphisms from a given base-pair mapping dataframe.

        This function processes a given DataFrame to identify features which have polymorphisms in
        their amino acid sequences, along with additional associated data such as the minimum and
        maximum positions of these features in the template and assembly sequences. It groups and
        aggregates the results for each feature and returns a DataFrame containing the processed
        information.

        Parameters:
            bp_map (pd.DataFrame): DataFrame containing mappings of base-pair positions, features,
            and other related information.

        Returns:
            pd.DataFrame: A DataFrame where each row corresponds to a feature with polymorphisms
            along with its aggregated attributes.
    """
    bp_map = bp_map.copy()
    bp_map['aa_change'] = bp_map.fillna('nan').template_aa != bp_map.fillna('nan').assembly_aa
    bp_map['template_pos_min'] = bp_map['template_pos_max'] = bp_map['template_pos'].replace(-1, None).ffill()
    bp_map['assembly_pos_min'] = bp_map['assembly_pos_max'] = bp_map['assembly_pos'].replace(-1, None).ffill()
    feature_polys = (bp_map.groupby(['feature_id', 'feature_name', 'feature_type', 'feature_strand'])
                     .agg({'poly_id': list_unique, 'poly_type': list_unique, 'aa_change': 'max',
                           'assembly_pos_min': 'min', 'assembly_pos_max': 'max'})
                     .reset_index().set_index('feature_id'))
    return feature_polys


def feature_from_row(row_index: int, row_data) -> SeqFeature:
    """
    Constructs and returns a SeqFeature object representing a polymorphism based on
    the given row data and its index.

    This function processes the row data to extract the type of polymorphism, its start and
    end positions in the assembly, and additional metadata. The function adjusts the positions
    for deletions as per the standard, initializes a SeqFeature with the resulting feature
    location and other metadata, and then populates any sequence-specific qualifiers found
    in the input data.

    Arguments:
        row_index: int
            Index of the row used to uniquely identify the feature.
        row_data: dict[tuple[str, str], Union[int, str]]
            Data for a single row, represented as a dictionary where keys are
            tuples specifying a category and subcategory, and values specify
            the data corresponding to that category.

    Returns:
        SeqFeature
            A constructed SeqFeature object based on the processed row data.

    Raises:
        KeyError
            If mandatory keys are not found in the row data.
    """
    poly_type = row_data[('poly_type', 'first')]
    start_pos, end_pos = row_data[('assembly_pos', 'min')], row_data[('assembly_pos', 'max')] + 1
    if poly_type == 'Deletion':
        start_pos += 1
    feature_loc = FeatureLocation(int(start_pos), int(end_pos))
    row_data: dict[tuple[str, str], Union[int, str]] = dict(row_data)
    feature = SeqFeature(feature_loc, 'Polymorphism', strand=0, id=f'Poly_{int(row_index):d}',
                         qualifiers={'Polymorphism Type': [poly_type], 'Effect': [row_data[('Effect', 'text')]]})
    for k, v in row_data.items():
        if k[1] == 'seq' and v:
            feature.qualifiers[k[0]] = [v]
    return feature


MUT_RE = re.compile(r'([A-Z*]+)(\d+)([A-Z*]+)')
INS_RE = re.compile(r'ins(\d+)([A-Z*]+)!rf(\d)')
DEL_RE = re.compile(r'del(\d+)\.?(\d*)!rf(\d)')
FS_RE = re.compile(r'fs(\d+)[+-](\d)')


def organize_mutations(feature_name: str) -> str:
    """
    Organizes a mutation string into a standardized format.

    This function processes a string containing feature name and associated
    mutations. It identifies, parses, and categorizes mutations into various
    types such as substitution, insertion, deletion, or nonsense mutations.
    Additionally, it handles reading frame adjustments in insertion and deletion
    mutations if applicable. The output is a reformatted string with mutations
    organized in a consistent structure.

    Parameters:
        feature_name: str
            A string representing the base feature name followed by a series
            of mutations separated by underscores.

    Returns:
        str
            A reformatted string where mutations are organized and structured
            in a specific format.
    """
    name_base, *mutations = feature_name.split('_')

    muts_by_residue: dict[int, tuple[str, str]] = {}
    for cur_mut in mutations:
        if m := MUT_RE.search(cur_mut):
            wt, res, mut = m.groups()
            ires = int(res)
            if ires in muts_by_residue:
                inverse_mut = f"{mut}{res}{wt}"
                if muts_by_residue[ires][1] == inverse_mut:
                    del muts_by_residue[ires]
                    continue
                elif m2 := MUT_RE.search(muts_by_residue[ires][1]):
                    wt2, re2, mut2 = m2.groups()
                    if wt == mut2:
                        wt = wt2
                    elif mut == wt2:
                        mut = mut2

            if '*' in mut:
                # remove any excess wt/mut residues
                stop_loc = mut.find('*')
                wt = wt[:stop_loc+1]
                mut = mut[:stop_loc+1]

                muts_by_residue[ires] = ('NONSENSE', f"{wt}{res}{mut}")
            else:
                muts_by_residue[ires] = ('SENSE', f"{wt}{res}{mut}")
        elif m := INS_RE.search(cur_mut):
            res, mut, frame = m.groups()
            muts_by_residue[int(res)] = ('INS', cur_mut)
        elif m := DEL_RE.search(cur_mut):
            res, stop_res, frame = m.groups()
            muts_by_residue[int(res)] = ('DEL', cur_mut)
        else:
            name_base += '_' + cur_mut

    out_parts = [name_base]
    reading_frame = 0
    for res_id in sorted(muts_by_residue.keys()):
        m_type, m_str = muts_by_residue[res_id]
        if m_type in ['INS', 'DEL'] and RF_DELIMITER in m_str:
            m_str, rf_str = m_str.split(RF_DELIMITER, 1)
            rf_int = int(rf_str)
            if rf_int != reading_frame:
                shift = (rf_int - reading_frame) % 3
                if reading_frame == 0 or rf_int == 0:
                    out_parts.append(f'fs{res_id:d}+{shift:d}')
                reading_frame = rf_int
            elif reading_frame == 0:
                out_parts.append(m_str)

        elif reading_frame == 0:
            out_parts.append(m_str)

            # if m_str[-1] == '*':  # truncation, no need to continue
            #     break

    return '_'.join(out_parts)


def assembly_analysis_pipeline(template_path: Path, assembly_path: Path, assembly_obj: PlasmidSeqAssembly,
                               session: Session, codon_table='RC63') -> pyd.Dseqrecord:
    """
    This function performs an end-to-end analysis pipeline for plasmid sequence assembly. It aligns
    the provided template and assembly sequences, identifies polymorphisms, maps features, and
    updates the corresponding database models with assembly-related features and polymorphisms.

    Args:
        template_path (Path): The file path to the template DNA sequence.
        assembly_path (Path): The file path to the assembly DNA sequence.
        assembly_obj (PlasmidSeqAssembly): An instance of the PlasmidSeqAssembly object representing
            the assembly process.
        session (Session): A database session used for persisting assembly-related data.
        codon_table (str): The codon translation table to be used for sequence analysis. Defaults
            to 'RC63'.

    Returns:
        pyd.Dseqrecord: A Biopython SeqRecord object with annotated features resulting from the
            assembly and polymorphism analysis.
    """
    bp_map, assembly_record, features_by_id = align_sequences(template_path, assembly_path, codon_table)

    # Add polymorphisms to the map
    poly_features = get_polymorphism_features(bp_map)
    polymorphism_features = [feature_from_row(i, v) for i, v in poly_features.iterrows() if
                             pd.notna(v[('assembly_pos', 'min')]) and pd.notna(v[('assembly_pos', 'max')])]
    assembly_record.features.extend(polymorphism_features)

    # Create polymorhphism models
    polymorphism_col_map = {
        ('template_pos', 'min'): 'wt_nt_start',
        ('template_pos', 'max'): 'wt_nt_end',
        ('assembly_pos', 'min'): 'assembly_nt_start',
        ('assembly_pos', 'max'): 'assembly_nt_end',
        ('Effect', 'text'): 'cds_effect',
        ('poly_type', 'first'): 'poly_type'
    }

    apms_by_poly_id: dict[float, AssemblyPolymorphism] = {}
    for poly_id, poly_data in poly_features.iterrows():
        apm_dict = {'assembly': assembly_obj}
        for k, v in zip(poly_data.index, poly_data):
            if k in polymorphism_col_map:
                apm_dict[polymorphism_col_map[k]] = v
        apm = AssemblyPolymorphism(**apm_dict)
        session.add(apm)
        apms_by_poly_id[poly_id] = apm

    # Add Features to the map and database
    feature_data = get_features_with_polymorphisms(bp_map)
    for cur_f_id, cur_wt_feature in features_by_id.items():
        try:
            wt_name = feature_data.loc[cur_f_id, 'feature_name']
        except KeyError:
            error_string = (f"Feature ID {cur_f_id:d} not found in feature data for {assembly_obj.assembly_name}.\n\n"
                            f"Feature:\n{str()}\n\nfeature_data:\n{feature_data.to_string()}")
            continue

        feature_type = feature_data.loc[cur_f_id, 'feature_type']
        f_start = feature_data.loc[cur_f_id, 'assembly_pos_min']
        f_end = feature_data.loc[cur_f_id, 'assembly_pos_max']
        af = AssemblyFeature(assembly=assembly_obj, wt_feature_name=wt_name, assembly_feature_name=wt_name,
                             deleted=f_start==f_end, feature_type=feature_type)
        if f_end > f_start:  # Add the feature to the map
            assy_loc = FeatureLocation(start=int(f_start), end=int(f_end)+1,
                                       strand=int(feature_data.loc[cur_f_id, 'feature_strand']))
            assy_feature = SeqFeature(assy_loc, feature_type)
            assy_feature.qualifiers.update(cur_wt_feature.qualifiers)
            af_poly_ids = [float(i) for i in feature_data.loc[cur_f_id, 'poly_id'].split(', ') if i]
            for pid in af_poly_ids:
                apm_for_current_poly = apms_by_poly_id[pid]
                af.polymorphisms.append(apm_for_current_poly)
                if apm_for_current_poly.cds_effect and apm_for_current_poly.cds_effect != 'None':
                    af.assembly_feature_name += '_' + apm_for_current_poly.cds_effect

            af.assembly_feature_name = organize_mutations(af.assembly_feature_name)
            assy_feature.qualifiers['label'] = [af.assembly_feature_name]
            assembly_record.features.append(assy_feature)

            if m := FS_RE.search(af.assembly_feature_name):
                res_id, shift = m.groups()
                af.frameshift_residue = int(res_id)
        session.add(af)

    return assembly_record


if __name__ == '__main__':
    print(organize_mutations('uox13_CYBJN_F6_Y200*_Q217*'))
    # names = [
    #     'GFP_M1C_D17G_C1M',
    #     'GFP_M1C_D17G_C1F',
    # ]
    #
    # for n in names:
    #     print(organize_mutations(n))
    #
    # t_file = Path(r"C:\Users\RobertWarden-Rothman\AppData\Roaming\JetBrains\PyCharm2023.2\docker\p_seq\tmp\1384_repick_1\GBFP-1384-0176\analysis_step\GBFP-1384-0176.1\GBFP-1384-0176.gb")
    # a_file = Path(r"C:\Users\RobertWarden-Rothman\AppData\Roaming\JetBrains\PyCharm2023.2\docker\p_seq\tmp\1384_repick_1\GBFP-1384-0176\analysis_step\GBFP-1384-0176.1\GBFP-1384-0176.1.gb")
    # o_file = Path(r"C:\Users\RobertWarden-Rothman\AppData\Roaming\JetBrains\PyCharm2023.2\docker\p_seq\tmp\1384_repick_1\GBFP-1384-0176\analysis_step\GBFP-1384-0176.1\GBFP-1384-0176.1_a.gb")
    #
    # from AppContainer.app.app import engine
    # from AppContainer.app.db_model import PlasmidSeqRun
    #
    # with Session(engine) as cur_session:
    #     assy_obj = cur_session.get(PlasmidSeqAssembly, 667)
    #     a_record = assembly_analysis_pipeline(t_file, a_file, assy_obj, cur_session)
    #     a_record.write(str(o_file))


class AssemblyTooShortError(ValueError):
    """
    Exception raised when an assembly sequence is too short.

    This exception is a specialized type of ValueError. It is used when the
    provided assembly sequence does not meet a predefined minimum length
    requirement.
    """
    pass
