from pathlib import Path
from itertools import repeat, count, starmap
from collections import deque
from typing import Iterable, TypeVar, Optional, Tuple, Dict

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.SeqFeature import FeatureLocation
from pandas import DataFrame
from pydna.dseqrecord import Dseqrecord
from pydna.seqfeature import SeqFeature
from sqlmodel import Session
from Bio.Align import PairwiseAligner
import pydna.all as pyd


if 'AppContainer' in __file__:
    from AppContainer.app.helpers import create_and_add_new_codon_table
    from AppContainer.app.db_model import PlasmidSeqAssembly, AssemblyFeature, AssemblyPolymorphism
else:
    from helpers import create_and_add_new_codon_table
    from db_model import PlasmidSeqAssembly, AssemblyFeature, AssemblyPolymorphism

create_and_add_new_codon_table('RC63', 'Standard', {'TAG': 'U'}, stop_codons=['TGA', 'TAA'])


T = TypeVar('T')


def threepeat(target: Iterable[T], max: int) -> Iterable[T]:
    if not isinstance(max, int) or max < 1:
        return []
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


def align_sequences(template_path: Path, assembly_path: Path,
                    codon_table='RC63') -> tuple[DataFrame, Dseqrecord, dict[int, SeqFeature]]:
    codon_table_obj = CodonTable.unambiguous_dna_by_name[codon_table]

    # Read in files
    template_record: pyd.Dseqrecord = pyd.read(str(template_path))
    assembly_record: pyd.Dseqrecord = pyd.read(str(assembly_path))
    if len(assembly_record) < len(template_record)/2:
        raise AssemblyTooShortError(f"{assembly_record.name} is too short to properly align")

    # Align circular templates
    aligner = PairwiseAligner()
    aligner.open_gap_score = -.1 * len(template_record)
    aligner.query_end_gap_score = 0.
    aligner.target_end_gap_score = 0.
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

            bp_map.loc[assembly_indices, 'assembly_res_id'] = list(threepeat(count(1), assembly_len))

            asy_seq = pyd.Dseq(''.join(bp_map.loc[assembly_indices, 'assembly_nt'].fillna('-')))
            if cur_feature.strand == -1:
                asy_seq = asy_seq.complement()
            asy_aa = asy_seq.translate(codon_table_obj) + ('X' * assembly_len)
            bp_map.loc[assembly_indices, 'assembly_aa'] = list(threepeat(asy_aa, assembly_len))

        else:
            bp_map.loc[template_indices, 'template_res_id'] = list(range(1, template_len + 1))
            bp_map.loc[assembly_indices, 'assembly_res_id'] = list(range(1, assembly_len + 1))
            bp_map.loc[template_indices, 'template_aa'] = np.nan
            bp_map.loc[template_indices, 'assembly_aa'] = np.nan

    return bp_map, assembly_record, features_by_id


def seq(x):
    return ''.join(x.dropna())


def find_effect(row_index, row_data):
    effect_text = None
    poly_type = row_data[('poly_type', 'first')]

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

        elif poly_type == 'Insertion':
            effect_text = f"del{res_min:d}"
            effect_text += a_aa or a_nt
    except ValueError as e:
        pass
    return effect_text


def get_polymorphism_features(bp_map: pd.DataFrame) -> pd.DataFrame:
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
                                          'template_res_id': ['min', 'max']
                                          })
    polys[('Effect', 'text')] = list(starmap(find_effect, polys.iterrows()))
    return polys


def list_unique(x: pd.Series) -> str:
    return ', '.join(x.dropna().unique().astype(str))


def get_features_with_polymorphisms(bp_map: pd.DataFrame) -> pd.DataFrame:
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
    poly_type = row_data[('poly_type', 'first')]
    start_pos, end_pos = row_data[('assembly_pos', 'min')], row_data[('assembly_pos', 'max')] + 1
    if poly_type == 'Deletion':
        start_pos += 1
    feature_loc = FeatureLocation(int(start_pos), int(end_pos))
    row_data: dict[tuple[str, str], int | str] = dict(row_data)
    feature = SeqFeature(feature_loc, 'Polymorphism', strand=0, id=f'Poly_{int(row_index):d}',
                         qualifiers={'Polymorphism Type': [poly_type], 'Effect': [row_data[('Effect', 'text')]]})
    for k, v in row_data.items():
        if k[1] == 'seq' and v:
            feature.qualifiers[k[0]] = [v]
    return feature


def assembly_analysis_pipeline(template_path: Path, assembly_path: Path, assembly_obj: PlasmidSeqAssembly,
                               session: Session, codon_table='RC63') -> pyd.Dseqrecord:
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
        wt_name = feature_data.loc[cur_f_id, 'feature_name']
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

            assy_feature.qualifiers['label'] = [af.assembly_feature_name]
            assembly_record.features.append(assy_feature)
        session.add(af)

    return assembly_record


if __name__ == '__main__':
    t_file = Path(r"C:\Users\RobertWarden-Rothman\AppData\Roaming\JetBrains\PyCharm2023.2\docker\p_seq\tmp\1384_repicks\GBFP-1384-0183\analysis_step\GBFP-1384-0183.1\GBFT-1384-0022.gb")
    a_file = Path(r"C:\Users\RobertWarden-Rothman\AppData\Roaming\JetBrains\PyCharm2023.2\docker\p_seq\tmp\1384_repicks\GBFP-1384-0183\analysis_step\GBFP-1384-0183.1\GBFP-1384-0183.1.gb")
    o_file = Path(r"C:\Users\RobertWarden-Rothman\AppData\Roaming\JetBrains\PyCharm2023.2\docker\p_seq\tmp\1384_repicks\GBFP-1384-0183\analysis_step\GBFP-1384-0183.1\GBFP-1384-0183.1_a.gb")

    from AppContainer.app.app import engine
    from AppContainer.app.db_model import PlasmidSeqRun

    with Session(engine) as cur_session:
        assy_obj = cur_session.get(PlasmidSeqAssembly, 667)
        a_record = assembly_analysis_pipeline(t_file, a_file, assy_obj, cur_session)
        a_record.write(str(o_file))


class AssemblyTooShortError(ValueError):
    pass
