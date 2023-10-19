import re
from collections import defaultdict
from copy import copy
from typing import Dict, List, Tuple, Union, Optional
from Bio import SeqIO, pairwise2
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dataclasses import dataclass, asdict
import pandas as pd


@dataclass(frozen=True)
class Mutation:
    position_start: int
    position_end: int
    mutation_type: str
    sequence_change: str


@dataclass(frozen=True)
class Annotation:
    position_start: int
    position_end: int
    annotation_type: str
    annotation_name: Optional[str] = None


@dataclass(frozen=True)
class AminoAcidEffect:
    annotation: Annotation
    first_aa: int
    original_aa: str
    mutated_aa: str
    mutation_details: List[Mutation]


@dataclass(frozen=True)
class FrameshiftMutation:
    mutation: Mutation
    annotation: Annotation


@dataclass(frozen=True)
class AADifference:
    position: int
    template_aa: str
    assembly_aa: str

    @classmethod
    def from_string(cls, difference_string: str) -> 'AADifference':
        """
        Parse a string in format "{template_aa}{position:d}{assembly_aa}" and return a Difference instance.
        """
        match = (re.match(r'(\D*)(\d+)(\D*)', difference_string) or
                 re.match(r'(.?)(\d+)(/+[12])', difference_string))
        if match:
            template_aa, position, assembly_aa = match.groups()
            return cls(int(position), template_aa, assembly_aa)
        else:
            raise ValueError(f"Could not parse '{difference_string}' as a difference.")

    @staticmethod
    def make_genotype(base_genotype: str, mutations: List['AADifference']):
        for m in mutations:
            base_genotype += m
            if '/' in m.assembly_aa:
                break
        return base_genotype

    def __add__(self, other) -> Union[str, 'AADifference']:
        if other == 0:
            return self  # make it compatible with sum()

        if isinstance(other, AADifference):
            return '_'.join(map(str, sorted([self, other])))

        if isinstance(other, str):
            gene_parts = []
            mut_parts = [self]
            for cur_part in other.split('_'):
                try:
                    mut_parts.append(AADifference.from_string(cur_part))
                except ValueError:
                    gene_parts.append(cur_part)
            return '_'.join(map(str, gene_parts + sorted(mut_parts)))

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        return self.__add__(other)

    def __lt__(self, other: 'AADifference') -> bool:
        """Less than function to make class sortable, sorting by position."""
        if other is None:
            return False
        return self.position < other.position

    def __str__(self):
        if '/' in self.assembly_aa:
            return f"{self.position:d}/{self.assembly_aa}"
        else:
            return f"{self.template_aa}{self.position:d}{self.assembly_aa}"

    def __repr__(self):
        return f"(AADifference: {str(self)})"


def create_and_add_new_codon_table(table_name, base_table, modifications, start_codons=None, stop_codons=None):
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


create_and_add_new_codon_table('RC63', 'Standard', {'TAG': 'U'}, stop_codons=['TGA', 'TAA'])


def map_positions(alignments: List[Tuple[str]]) -> Dict[int, int]:
    """
    Map positions between the template and the assembly sequences.

    :param alignments: List of resulting alignments from the Biopython pairwise alignment function.
    :return: A dictionary representing the mapping from positions in the template sequence to the assembly sequence.
    """
    position_map = dict()
    for alignment in alignments:
        t_seq, a_seq, _, _, _ = alignment
        t_actual_pos = 0
        a_actual_pos = 0
        for i in range(len(t_seq)):
            if t_seq[i] != '-':
                if a_seq[i] != '-':
                    position_map[t_actual_pos] = a_actual_pos
                t_actual_pos += 1
            if a_seq[i] != '-':
                a_actual_pos += 1
        return position_map


def propagate_annotations(template_seq_record: SeqRecord, assembly_seq_record: SeqRecord, position_map: Dict[int, int]) -> SeqRecord:
    """
    Propagates annotations from the template sequence record to the assembly sequence record.

    :param template_seq_record: SeqRecord object of the template sequence.
    :param assembly_seq_record: SeqRecord object of the assembly sequence.
    :param position_map: Mapping of positions from the template sequence to the assembly sequence.
    :return: Assembly SeqRecord object with updated features.
    """
    for feature in template_seq_record.features:
        # Check if the feature is within the limits of the assembly sequence
        if feature.location.start in position_map and feature.location.end in position_map:
            # Map the start and end positions of the feature to the assembly sequence
            start = position_map[feature.location.start]
            end = position_map[feature.location.end]
            new_feature = copy(feature)
            new_feature.location = FeatureLocation(start, end, new_feature.location.strand)
            assembly_seq_record.features.append(new_feature)
    return assembly_seq_record


def find_indels(map_positions: Dict[int, int]) -> List[Mutation]:
    """
    Find the positions of insertions and deletions.

    :param map_positions: Mapping of positions from the template sequence to the assembly sequence.
    :return: A list of Mutation object representing the insertions and deletions.
    """
    indels = []

    # Get sorted positions
    positions = sorted(map_positions.items())

    for i in range(1, len(positions)):
        t_gap = positions[i][0] - positions[i - 1][0]
        a_gap = positions[i][1] - positions[i - 1][1]

        # Check if there are gaps in the template sequence
        if t_gap > 1:
            # This is a deletion (a_gap should be 1 if it's aligned)
            indels.append(Mutation(positions[i - 1][0] + 1, positions[i][0] - 1, "Deletion", ""))

        # Check if there are gaps in the assembly sequence
        if a_gap > 1:
            # This is an insertion (t_gap should be 1 if it's aligned)
            indels.append(Mutation(positions[i - 1][1] + 1, positions[i][1] - 1, "Insertion", ""))

    return indels


def find_mutated_regions(position_map: Dict[int, int], template_seq_record: SeqRecord, assembly_seq_record: SeqRecord) -> List[Mutation]:
    """
    Function to identify mutations in the assembly sequence compared to the template sequence

    Args:
    position_map (Dict[int, int]): Mapping of positions from the template sequence to the assembly sequence
    template_seq_record (SeqRecord): SeqRecord of the template sequence
    assembly_seq_record (SeqRecord): SeqRecord of the Assembly sequence

    Returns:
    List: List of Mutation objects representing the mutated regions in the assembly sequence
    """

    mutated_regions = []  # List to hold the mutated regions
    region_start = None  # Store the start position of a mutated region
    original_sequence = ''  # Keep track of the original sequence
    mutated_sequence = ''  # Keep track of the mutated sequence

    # Iterate over the sorted position mapping
    for t_position, a_position in sorted(position_map.items()):

        # Get the nucleotide bases at the current positions in the template and assembly sequences
        t_base = template_seq_record.seq[t_position]
        a_base = assembly_seq_record.seq[a_position]

        # Check if the nucleotide bases are different
        if t_base != a_base:

            # Start a new mutated region if we haven't started one
            if region_start is None:
                region_start = t_position

            # Update the original and mutated sequence strings
            original_sequence += t_base
            mutated_sequence += a_base

        else:

            # If we are in a mutated region and encounter a matching base, we end the mutated region
            if region_start is not None:
                mutated_regions.append(Mutation(region_start, t_position - 1, 'Mutation', original_sequence+' -> '+mutated_sequence))

                # Reset the mutated region start position and the sequences
                region_start = None
                original_sequence = ''
                mutated_sequence = ''

    # If we have started a mutated region and haven't ended it, we end it
    if region_start is not None:
        mutated_regions.append(Mutation(region_start, t_position, 'Mutation', original_sequence+' -> '+mutated_sequence))

    return mutated_regions


def translate_and_compare(template_seq_record: SeqRecord,
                          assembly_seq_record: SeqRecord,
                          indels: List[Mutation],
                          codon_table: str = "Standard") -> Dict[str, List[AADifference]]:
    """
    Compare the translation of each CDS feature in the template sequence record and the assembly sequence record.
    Return a dictionary with differences between the two sequences.

    :param template_seq_record: The template sequence record containing CDS features.
    :param assembly_seq_record: The assembly sequence record containing CDS features.
    :param indels: A list of mutations.
    :param codon_table: The codon table to use for translation. Default is "Standard".
    :return: A dictionary with differences between the template and assembly sequences.

    """
    codon_table_obj = CodonTable.unambiguous_dna_by_name[codon_table]
    differences = defaultdict(list)

    for tmp_feature in filter(lambda f: f.type == 'CDS', template_seq_record.features):
        tmp_qualifier = tmp_feature.qualifiers.get('label', [None])[0]
        for asy_feature in filter(lambda f: f.type == 'CDS' and f.qualifiers.get('label', [None])[0] == tmp_qualifier,
                                  assembly_seq_record.features):
            if tmp_feature.location.strand not in [1, 0]:
                tmp_feature_seq = tmp_feature.location.extract(template_seq_record.seq).reverse_complement()
            else:
                tmp_feature_seq = tmp_feature.location.extract(template_seq_record.seq)

            if asy_feature.location.strand not in [1, 0]:
                ass_feature_seq = asy_feature.location.extract(assembly_seq_record.seq).reverse_complement()
            else:
                ass_feature_seq = asy_feature.location.extract(assembly_seq_record.seq)

            tmp_feature_seq_aa = tmp_feature_seq.translate(table=codon_table_obj)
            ass_feature_seq_aa = ass_feature_seq.translate(table=codon_table_obj)

            # Check if this CDS contains an indel resulting in a frameshift mutation
            frameshift = next((indel for indel in indels
                               if indel.position_start >= tmp_feature.location.start
                               and indel.position_end <= tmp_feature.location.end
                               and (fs_shift := (indel.position_end - indel.position_start) % 3) != 0), None)

            if frameshift is not None:
                if tmp_feature.location.strand not in [0, 1]:
                    start_pos = (tmp_feature.location.end - frameshift.position_end) // 3
                else:
                    start_pos = (frameshift.position_start - tmp_feature.location.start) // 3
                differences[tmp_qualifier].append(AADifference(start_pos, '/', f'+{fs_shift:d}'))
                continue

            if tmp_feature_seq_aa != ass_feature_seq_aa:
                alignment = pairwise2.align.globalxs(tmp_feature_seq_aa, ass_feature_seq_aa, -2, -.1)[0]

                for i, (p, a) in enumerate(zip(alignment.seqA, alignment.seqB)):
                    if p != '-' and a != '-' and p != a:
                        differences[tmp_qualifier].append(AADifference(i + 1, p, a))

    return differences


def align_sequences(template_file: str, assembly_file: str, codon_table='RC63') -> Tuple[SeqRecord, Dict[int, int],
    pd.DataFrame, pd.DataFrame]:
    """
    :param template_file: Path to the template (genbank) file.
    :param assembly_file: Path to the assembly (fasta) file.
    :param codon_table: Codon table to use for translation. Defaults to 'RC63'.
    :return: Tuple containing the aligned assembly sequence, position map, mutation dataframe, and CDS dataframe.

    This method aligns the assembly sequence to the template sequence, propagates annotations from the template to the
    assembly, identifies indels and mutations between the two sequences, adds indels and mutations as features to the
    assembly sequence, calculates CD effects, updates CDS features in the assembly sequence with genotype information,
    and returns the aligned assembly sequence, position map, mutation dataframe, and CDS dataframe.

    Example usage:
        template_file = 'template.gb'
        assembly_file = 'assembly.fasta'
        codon_table = 'RC63'
        aligned_sequence, position_map, mutation_df, cds_df = align_sequences(template_file, assembly_file, codon_table)
    """
    # Reading the template (genbank) and assembly (fasta) files
    template_seq_record = SeqIO.read(template_file, "genbank")
    assembly_seq_record = next(SeqIO.parse(assembly_file, "fasta"))
    assembly_seq_record_rc = assembly_seq_record.reverse_complement()

    # Extracting sequences and converting to string
    template_seq = str(template_seq_record.seq)
    assembly_seq = str(assembly_seq_record.seq) * 2
    assembly_seq_rc = str(assembly_seq_record_rc.seq) * 2

    # Performing the alignment (Using global alignment in this case)
    alignments = pairwise2.align.globalxs(template_seq, assembly_seq, -2, -.1)
    alignments_rc = pairwise2.align.globalxs(template_seq, assembly_seq_rc, -2, -.1)

    # Check for a reverse-complement assembly sequence
    if alignments_rc[0][2] > alignments[0][2]:
        alignments = alignments_rc
        assembly_seq_record = assembly_seq_record_rc
        assembly_seq = assembly_seq_rc

    # Rotate the assembly to match the template
    shift = len(alignments[0][0]) - len(alignments[0][0].lstrip('-'))
    if shift > 0:
        assembly_seq_record = assembly_seq_record[shift:] + assembly_seq_record[:shift]
        assembly_seq = str(assembly_seq_record.seq)

        # Realign
        alignments = pairwise2.align.globalxs(template_seq, assembly_seq, -2, -.1)

    # Map positions between the template and assembly sequences
    position_map = map_positions(alignments)

    # Propagate annotations from the template sequence to the assembly sequence
    assembly_seq_record = propagate_annotations(template_seq_record, assembly_seq_record, position_map)
    assembly_seq_record.annotations.update(dict(topology='circular', molecule_type='DNA', data_file_division='SYN'))
    assembly_seq_record.name = f"pseq_{template_seq_record.name}"[:16].replace(' ', '_')

    # find indels & add to map
    indels = find_indels(position_map)
    mutations = find_mutated_regions(position_map, template_seq_record, assembly_seq_record)
    differences = indels + mutations

    for c_pm in differences:
        if c_pm.mutation_type == 'Deletion':
            start_position = position_map[c_pm.position_start - 1] + 1
            end_position = position_map[c_pm.position_end + 1]
        else:
            start_position = position_map[c_pm.position_start]
            end_position = position_map[c_pm.position_end] + 1
        new_feature = SeqFeature(FeatureLocation(start_position, end_position), 'Polymorphism',
                                 qualifiers={'Polymorphism Type': [c_pm.mutation_type],
                                             'Change': [c_pm.sequence_change]})
        assembly_seq_record.features.append(new_feature)

    # alter genotypes
    cds_effects = translate_and_compare(template_seq_record, assembly_seq_record, indels, 'RC63')
    cur_cds: SeqFeature
    for cur_cds in (f for f in assembly_seq_record.features if f.type == 'CDS'):
        f_label = cur_cds.qualifiers.get('label', [None])[0]
        if not f_label or f_label not in cds_effects:
            continue

        cur_cds.qualifiers['template_label'] = [f_label]
        cur_cds.qualifiers['additional_mutations'] = [sum(cds_effects[f_label])]
        cur_cds.qualifiers['label'] = [AADifference.make_genotype(f_label, cds_effects[f_label])]

    # make polymorphism dataframe
    mut_df = pd.DataFrame(map(asdict, differences))
    cds_df_list = []
    for cds_name, diffs in cds_effects.items():
        cur_df = pd.DataFrame(map(asdict, diffs))
        cur_df.insert(0, 'CDS', cds_name)
        cds_df_list.append(cur_df)
    cds_df = pd.concat(cds_df_list) if cds_df_list else pd.DataFrame()

    return assembly_seq_record, position_map, mut_df, cds_df


if __name__ == '__main__':

    align_sequences(
        r"C:\Users\RobertWarden-Rothman\AppData\Roaming\JetBrains\PyCharm2023.2\docker\p_seq\tmp\gff_conversion\pGRO-K1196.gbk",
        r"C:\Users\RobertWarden-Rothman\AppData\Roaming\JetBrains\PyCharm2023.2\docker\p_seq\tmp\Qmivz9EPOHku8c99k509K\assembly2.fasta"
    )
