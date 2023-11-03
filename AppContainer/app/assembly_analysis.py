from pathlib import Path
from itertools import repeat, count
from collections import deque
from typing import Iterable, TypeVar, Optional

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from pydna.seqfeature import SeqFeature
from sqlmodel import Session
from Bio.Align import PairwiseAligner
import pydna.all as pyd

if 'AppContainer' in __file__:
    from AppContainer.app.helpers import create_and_add_new_codon_table
else:
    from helpers import create_and_add_new_codon_table

create_and_add_new_codon_table('RC63', 'Standard', {'TAG': 'U'}, stop_codons=['TGA', 'TAA'])


T = TypeVar('T')


def threepeat(target: Iterable[T], max: int = 0) -> Iterable[T]:
    cur_count = 0
    try:
        for i in target:
            for j in repeat(i, 3):
                yield j
                cur_count += 1
                if cur_count == max:
                    raise StopIteration()
    except StopIteration:
        pass


def triplets(target: Iterable[T]) -> Iterable[T]:
    while target:
        triplet, target = target[:3], target[3:]
        yield triplet


def find_frameshift(seq: str) -> Optional[int]:
    for i, cur_codon in enumerate(triplets(seq)):
        if '-' in cur_codon:
            return i + 1
    return None


def align_sequences(template_path: Path, assembly_path: Path, codon_table='RC63'):
    codon_table_obj = CodonTable.unambiguous_dna_by_name[codon_table]

    # Read in files
    template_record: pyd.Dseqrecord = pyd.read(str(template_path))
    assembly_record: pyd.Dseqrecord = pyd.read(str(assembly_path))

    # Align circular templates
    aligner = PairwiseAligner()
    aligner.internal_open_gap_score = -10
    alignments = aligner.align(str(template_record.seq), str(assembly_record.seq) * 2)
    alignment = sorted(alignments)[0]

    # Determine sequence shift and re-align
    shift = dict(alignment.indices.T)[0]
    assembly_record = assembly_record.shifted(shift)

    final_alignments = aligner.align(template_record.seq, assembly_record.seq)
    best_alignment = sorted(final_alignments)[0]

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
    diffs = polys_only['assembly_pos'].diff()
    poly_ids = (diffs != 1).cumsum()
    bp_map.loc[polys_only.index, 'poly_id'] = poly_ids

    # Add features to the table
    features_by_id: dict[int, SeqFeature] = dict(enumerate(template_record.sorted_features()))

    for feature_id, cur_feature in features_by_id.items():
        feature_indices = []
        for cur_loc in cur_feature.location.parts:
            start_index = bp_map[bp_map['template_pos'] == cur_loc.start].index.min()
            end_index = bp_map[bp_map['template_pos'] == cur_loc.end].index.max()
            feature_indices.extend(range(start_index, end_index))

        template_indices: list[int] = list(bp_map.loc[feature_indices]['template_nt'].dropna().index.values)
        assembly_indices: list[int] = list(bp_map.loc[feature_indices]['assembly_nt'].dropna().index.values)

        if cur_feature.location.strand == -1:
            feature_indices.reverse()
            template_indices.reverse()
            assembly_indices.reverse()

        bp_map.loc[feature_indices, 'feature_id'] = feature_id
        bp_map.loc[feature_indices, 'feature_name'] = cur_feature.qualifiers['label'][0]
        bp_map.loc[feature_indices, 'feature_type'] = cur_feature.type

        template_len = len(template_indices)
        assembly_len = len(assembly_indices)

        if cur_feature.type == 'CDS':
            bp_map.loc[template_indices, 'template_res_id'] = list(threepeat(count(1), template_len))

            wt_seq = pyd.Dseq(''.join(bp_map.loc[template_indices, 'template_nt']))
            if cur_feature.strand == -1:
                wt_seq = wt_seq.complement()
            wt_aa = wt_seq.translate(codon_table_obj) + ('X' * template_len)
            bp_map.loc[template_indices, 'template_aa'] = list(threepeat(wt_aa, template_len))

            bp_map.loc[assembly_indices, 'assembly_res_id'] = list(threepeat(count(1), assembly_len))

            asy_seq = pyd.Dseq(''.join(bp_map.loc[assembly_indices, 'assembly_nt'].fillna('-')))
            if cur_feature.strand == -1:
                asy_seq = asy_seq.complement()
            asy_aa = asy_seq.translate(codon_table_obj) + ('X' * assembly_len)
            bp_map.loc[assembly_indices, 'assembly_aa'] = list(threepeat(asy_aa, assembly_len))

        else:
            bp_map.loc[template_indices, 'template_res_id'] = list(range(1, template_len + 1))
            bp_map.loc[assembly_indices, 'assembly_res_id'] = list(range(1, assembly_len + 1))

    print(alignment)


if __name__ == '__main__':
    t_file = Path(r"C:\Users\RobertWarden-Rothman\PycharmProjects\denovo_plasmid_assembly\tests\Uox17_F9_MW.gb")
    a_file = Path(r"C:\Users\RobertWarden-Rothman\PycharmProjects\denovo_plasmid_assembly\tests\Uox17_F9_MW.1.gb")

    from AppContainer.app.app import get_session

    align_sequences(t_file, a_file)
