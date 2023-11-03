#! /usr/bin/env python3
#
# denovo_plasmid_assembly.py
# de novo assemblies for circular plasmids from paired-end Illumina sequences.
#
# Copyright 2019 Mark F. Rogers (rogers@genofab.com)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from optparse import OptionParser
import gzip
import os
import shutil
import subprocess
import sys
import time

# Quality filtering defaults
DEFAULT_WINDOW = 9

# DNA complements plus nucleotide ambiguity codes (IUPAC):
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y': 'R', 'R': 'Y', 'W': 'W', 'S': 'S',
              'K': 'M', 'M': 'K', 'D': 'H', 'V': 'B', 'H': 'D', 'B': 'V', 'X': 'X', 'N': 'N',
              'a': 'T', 'c': 'G', 'g': 'C', 't': 'A', 'y': 'R', 'r': 'Y', 'w': 'W', 's': 'S',
              'k': 'M', 'm': 'K', 'd': 'H', 'v': 'B', 'h': 'D', 'b': 'V', 'x': 'X', 'n': 'N'}

FASTX_INDICATORS = '>@'
FASTA_INDICATOR = '>'
FASTQ_EXTS = ['.fastq', '.fq']
GZIP_EXTS = ['.gz', '.gzip']

# Format-specific constants
FASTA_HEADER_CHAR = '>'

# External program values:
TRIMMOMATIC_PATH='./Trimmomatic-0.39/trimmomatic-0.39.jar'
TRIM_FORMAT = "java -jar %s PE -phred33"
TRIMMOMATIC_DEFAULT_QUALITY = 35
TRIMMOMATIC_MIN_LENGTH = 50

# Example: fastp -i $FWD -o fwd_filtered.fastq -I $REV -O rev_filtered.fastq -5 -3 -W 9 -M 38 -l $MIN_LEN -e 38 -w 8
FASTP_FWD = 'fwd_filtered.fastq'
FASTP_REV = 'rev_filtered.fastq'
FASTP_DEFAULT_QUALITY = 38
FASTP_MIN_LENGTH = 150

UNICYCLER_LOG = 'run_unicycler.log'
UNICYCLER_ASSEMBLY = 'assembly.fasta'
FINAL_ASSEMBLY = 'final_assembly.fasta'
UNICYCLER_OPTIONS = '--keep 3 --no_rotate'

# CDF estimated from bootstrapped combinations of 'primary' and 'contaminant' sequences with contaminant
# proportions ranging from 5% to 50% of total reads.  For assemblies with a single sequence, reads were compared
# to the assembly to find exact matches and mismatches.  All data sets yield some proportion of mismatches, but
# contaminated data sets yield higher proportions than clean data.  These (value,percentile) pairs reflect the
# estimated distribution for contaminated data sets.
CONTAMINATED_CDF = [(11.15,1), (11.43,2), (11.6,3), (11.8,4), (12.0,5), (12.2,6), (12.3,7), (12.4,8), (12.5,9), (12.6,10),
                    (12.7,11), (12.8,13), (12.9,14), (13.0,15), (13.1,16), (13.12,17), (13.2,18), (13.3,19), (13.4,20), (13.47,21),
                    (13.5,22), (13.6,23), (13.64,24), (13.7,25), (13.8,27), (13.9,28), (13.98,29), (14.1,30), (14.10,31), (14.2,32),
                    (14.3,34), (14.4,35), (14.5,37), (14.6,38), (14.7,39), (14.73,40), (14.8,41), (14.9,42), (15.0,44), (15.1,45),
                    (15.2,46), (15.3,47), (15.4,48), (15.5,49), (15.6,50), (15.7,51), (15.70,52), (15.8,53), (15.9,54), (16.0,55),
                    (16.09,56), (16.1,57), (16.2,58), (16.3,59), (16.4,61), (16.5,62), (16.6,63), (16.7,64), (16.8,65), (16.9,66),
                    (17.0,67), (17.1,68), (17.2,69), (17.21,70), (17.3,71), (17.4,72), (17.5,73), (17.56,74), (17.63,75), (17.7,76),
                    (17.8,77), (17.9,78), (18.0,79), (18.14,80), (18.3,81), (18.4,82), (18.55,83), (18.7,84), (18.8,85), (18.9,86),
                    (19.0,87), (19.3,88), (19.49,89), (19.80,90), (20.07,91), (20.34,92), (20.7,93), (20.90,94), (21.39,95), (21.91,96),
                    (22.4,97), (23.4,98), (24.68,99), (29.78,100)]


def assure_directory(d, verbose=False):
    """Guarantee that the given directory exists, or raise an exception if there is a problem with it."""
    if not os.path.isdir(d):
        if verbose:
            sys.stderr.write('Creating directory %s\n' % d)
        os.makedirs(d)
    return d


def compute_mismatch_percentage(reference_file, queries, max_readlen, verbose=False):
    """Finds reads in the query files that fail to match the reference sequence exactly."""
    # We use an expanded reference sequence to account for wrap-around
    # (circular assemblies), so this method only needs to detect reads not found
    # anywhere in the reference sequence -- simpler and faster than Smith-Waterman.
    expanded_assembly = load_expanded_assembly(reference_file, max_readlen)
    reverse_assembly = reverse_complement(expanded_assembly)
    mismatches = 0
    for q in queries:
        if q in expanded_assembly or q in reverse_assembly:
            continue
        mismatches += 1

    return mismatches * 100.0 / len(queries)


def denovo_assembly(fwd_fastq, rev_fastq, output_dir, verbose=False):
    """Runs Unicycler to produce a de novo assembly from the given paired-end FASTQ files
    and places the results in the output directory."""
    cmd = 'unicycler -1 %s -2 %s -o %s %s' % (fwd_fastq, rev_fastq, output_dir, UNICYCLER_OPTIONS)
    run_command(cmd, verbose=verbose)
    return os.path.join(output_dir, UNICYCLER_ASSEMBLY)


def estimate_likelihood(misalign_pct, cdf):
    """Returns the estimated likelihood (area under cdf) for the given misalignment percentage."""
    for (pct, percentile) in cdf:
        if misalign_pct < pct:
            return percentile
    return 100


def ezopen(fileName) :
    """Open files even if they're gzipped."""
    if not (os.path.exists(fileName) and os.path.isfile(fileName)):
        raise ValueError('file does not exist at %s' % fileName)

    handle = gzip.open(fileName, mode='rt')
    try :
        line = handle.readline()
        handle.close()
        return gzip.open(fileName, mode='rt')
    except :
        return open(fileName)


def fastp(fwd_fastq, rev_fastq, min_length, window_size, min_quality, dest_dir, verbose=False):
    """Runs fastp to filter FASTQ reads using a sliding window and a minimum final length."""

    fwd_base = find_file_prefix(fwd_fastq)
    rev_base = find_file_prefix(rev_fastq)

    fwd_filtered = os.path.join(dest_dir, FASTP_FWD)
    rev_filtered = os.path.join(dest_dir, FASTP_REV)

    # Example: fastp -i forward.fastq -o fwd_filtered.fastq -I reverse.fastq -O rev_filtered.fastq -5 -3 -W 9 -M 38 -l 150 -e 38 -w 8
    cmd = "fastp -i %s -o %s -I %s -O %s -5 -3 -W %d -M %d -l %d -e %d" % \
            (fwd_fastq, fwd_filtered, rev_fastq, rev_filtered, window_size, min_quality, min_length, min_quality)
    run_command(cmd, verbose=verbose)

    return fwd_filtered, rev_filtered


def find_file_prefix(fastq_file):
    """Finds the file prefix: for a file named foo.fastq, foo.fastq.gz, foo.fq or foo.fq.gz this returns just 'foo'.
    This is inelegant but should work for most FASTQ/gzipped FASTQ files."""
    result = os.path.basename(fastq_file)
    for gz_ext in GZIP_EXTS:
        if result.endswith(gz_ext):
            result = result[:-len(gz_ext)]

    for fq_ext in FASTQ_EXTS:
        if result.endswith(fq_ext):
            result = result[:-len(fq_ext)]
    return result


def load_expanded_assembly(assembly_file, max_readlen):
    """Loads a reference assembly and returns an expanded version based on the extension length."""
    extension_length = max_readlen - 1
    result = None
    for n, s, ignore in read_fastx(assembly_file):
        if result is None:
            if len(s) < extension_length:
                raise ValueError('Assembly in %s is shorter than queries!' % assembly_file)
            result = s + s[:extension_length]
        else:
            raise ValueError('Found multiple sequences in assembly file')
    return result


def load_queries(file_list):
    """Load all the query sequences from the given list.  Returns the sequences along with their maximum size."""
    max_length = 0
    queries = []
    for query_file in file_list:
        try:
            for n, s, q in read_fastx(query_file):
                max_length = max(max_length, len(s))
                queries.append(s)
        except IndexError:
            raise ValueError('Unable to load sequences from %s\n' % query_file)
    return queries, max_length


def post_process_assembly(input_fasta, output_fasta, prefix):
    """Augments headers in a FASTA file using the given prefix
    and returns the number of sequences found in the file."""
    records = [s.strip() for s in open(input_fasta)]
    fragment_count = 0
    total_bases = 0
    ostr = open(output_fasta, 'w')
    for s in records:
        if s.startswith(FASTA_HEADER_CHAR):
            ostr.write('%s%s-%s\n' % (FASTA_HEADER_CHAR, prefix, s[1:]))
            fragment_count += 1
        else:
            ostr.write('%s\n' % s)
            total_bases += len(s)
    ostr.close()
    return fragment_count, total_bases


def read_count_warnings(queries, opts):
    """Convenience method for reporting possibly too few or too many reads to run an assembly."""
    half_size = len(queries)/2
    if half_size < 1000:
        print('\n  ** Filtering with minimum quality %d and window size %d yields %d filtered read pairs.' % (opts.min_quality, opts.filter_window, half_size))
        print('  ** Consider decreasing min-quality (--min-quality) or increasing the filter window (-w)\n')
    elif half_size > 10000:
        print('\n  ** Filtering with minimum quality %d and window size %d yields %d filtered read pairs.' % (opts.min_quality, opts.filter_window, half_size))
        print('  ** Consider increasing min-quality (--min-quality) or decreasing the filter window (-w)\n')


def read_fastx(fastx_path):
    """ read a FASTA/FASTQ sequence file """
    # Note: converts sequences to upper-case to avoid 'a' vs. 'A' issues
    def read_fasta(f):
        """ read a fasta file """
        seq_id = ''
        sequence = ''
        for l in f:
            if l.startswith('>'):
                if sequence:
                    yield seq_id, sequence, ''
                seq_id = l.strip()[1:].split()[0]
                sequence = ''
            else:
                sequence += l.strip()

        yield seq_id, sequence.upper(), None

    def read_fastq(f):
        """ read a fastq file """
        # Each record consists of four lines:
        #    sequence id / sequence / sequence id (duplicate) / quality encoding
        for l in f:
            seq_id = l.strip()[1:].split()[0]
            sequence = f.readline()
            s3 = f.readline()  # ignore
            phred_score = f.readline()
            yield seq_id.strip(), sequence.strip().upper(), phred_score.strip()

    # test if fasta or fastq
    fastx_stream = ezopen(fastx_path)
    first_line = fastx_stream.readline()
    assert(type(first_line) == str)
    fastx_stream.close()

    first_char = first_line[0]
    if first_char not in FASTX_INDICATORS:
        raise ValueError('%s is not a FASTA/FASTQ file' % fastx_path)

    fastx_stream = ezopen(fastx_path)

    # iterate over sequences in the file
    if first_char == FASTA_INDICATOR:
        for sId, sSeq, ignore in read_fasta(fastx_stream):
            yield sId, sSeq, ignore
    else:
        for sId, sSeq, sQual in read_fastq(fastx_stream):
            yield sId, sSeq, sQual


def reverse_complement(s):
    """Returns the reverse-complement of a DNA sequence."""
    return ''.join([COMPLEMENT[x] for x in s[::-1]])


def run_command(s, err_stream=subprocess.DEVNULL, out_stream=subprocess.DEVNULL, cwd=None, verbose=False):
    """Simple wrapper method for running shell commands."""
    if verbose:
        sys.stderr.write(time_string('%s\n' % s))

    status = subprocess.call(s, shell=True, stderr=err_stream, stdout=out_stream, cwd=cwd, timeout=14*60)

    if status != 0:
        sys.stderr.write("Failed to execute the following command:\n%s\n" % s)

    return status


def time_string(s, format_string='%X', lf=False):
    """Returns the input string with a user-readable timestamp prefix."""
    timestamp = time.strftime(format_string, time.localtime())
    result = '%s %s' % (timestamp, s)
    if lf:
        result += '\n'
    return result


def trimmomatic(fwd_fastq, rev_fastq, min_length, window_size, min_quality, dest_dir, jar_path, verbose=False):
    """Runs Trimmomatic to filter FASTQ reads using a sliding window and a minimum final length."""

    fwd_base = find_file_prefix(fwd_fastq)
    rev_base = find_file_prefix(rev_fastq)

    fwd_paired = os.path.join(dest_dir, '%s_paired.fastq' % fwd_base)
    fwd_unpaired = os.path.join(dest_dir, '%s_unpaired.fastq' % fwd_base)
    rev_paired = os.path.join(dest_dir, '%s_paired.fastq' % rev_base)
    rev_unpaired = os.path.join(dest_dir, '%s_unpaired.fastq' % rev_base)

    options_string = 'SLIDINGWINDOW:%d:%.3f MINLEN:%d' % (window_size, min_quality, min_length)

    trim_prefix = TRIM_FORMAT % jar_path

    cmd = '%s %s %s %s %s %s %s %s' % (trim_prefix, fwd_fastq, rev_fastq, fwd_paired, fwd_unpaired, rev_paired, rev_unpaired, options_string)
    run_command(cmd, verbose=verbose)

    return fwd_paired, rev_paired


def validate_file(path):
    """Standard method for validating file paths."""
    if not path:
        raise Exception("'%s' is not a valid file path; exiting." % path)

    if not os.path.exists(path):
        raise Exception("File '%s' not found; exiting." % path)


USAGE = """%prog fwd-fastq rev-fastq assembly-name [options]

Produces a de novo assembly for circular plasmids from paired-end Illumina sequences.

Users must provide 'fwd-fastq' and 'rev-fastq' FASTQ files that contain the forward
and reverse reads, respectively.

Users must also provide a name for the assembly that will be used for the sequence
header(s) in the output assembly file.

Prerequisites:

    Trimmomatic (http://www.usadellab.org/cms/index.php?page=trimmomatic, also requires Java)
    - or -
    fastp (https://github.com/OpenGene/fastp)

    AND

    Unicycler (https://github.com/rrwick/Unicycler) """


# Establish command-line options:
if __name__ == '__main__':
    parser = OptionParser(usage=USAGE)
    parser.add_option('-d', dest='output_dir', default='denovo_assembly', help='Output directory [default: %default]')
    parser.add_option('-m', dest='min_length', default=None, help='Minimum filtered read length [default: Trimmomatic=50, fastp=150]', type='int')
    parser.add_option('-w', dest='filter_window', default=DEFAULT_WINDOW, help='Filtering sliding window size [default: %default]', type='int')
    parser.add_option('-v', dest='verbose', default=False, help='Verbose mode [default: %default]', action='store_true')
    parser.add_option('--min-quality', dest='min_quality', default=None, help='Minimum filtering quality [default: Trimmomatic=35, fastp=38]', type='float')
    parser.add_option('--jar', dest='jar_path', default=TRIMMOMATIC_PATH, help='(Trimmomatic) Path to JAR file [default: %default]')
    parser.add_option('--fastp', dest='fastp', default=False, help='Use fastp trimming [default: %default]', action='store_true')
    opts, args = parser.parse_args(sys.argv[1:])


    if len(args) != 3 :
        parser.print_help()
        sys.exit(1)

    unicycler_path = shutil.which('unicycler')
    if not unicycler_path:
        parser.print_help()
        sys.stderr.write('\n** Cannot find unicycler in your path!  Is it installed?\n\n')
        sys.exit(1)

    for f in args[:2]:
        validate_file(f)

    [fwd_fastq, rev_fastq, assembly_name] = args

    assure_directory(opts.output_dir)

    start_time = time.time()

    if opts.fastp:
        fastp_path = shutil.which('fastp')
        if not unicycler_path:
            parser.print_help()
            sys.stderr.write('\n** Cannot find fastp in your path!  Is it installed?\n\n')
            sys.exit(1)

        print(time_string('Using fastp for trimming: no contamination analysis will be performed.'))
        if opts.min_quality is None:
            opts.min_quality = FASTP_DEFAULT_QUALITY

        if opts.min_length is None:
            opts.min_length = FASTP_MIN_LENGTH

        fwd_filtered, rev_filtered = fastp(fwd_fastq, rev_fastq, opts.min_length,
                                              opts.filter_window,
                                              opts.min_quality,
                                              opts.output_dir,
                                              verbose=opts.verbose)
    else:
        java_path = shutil.which('java')
        if not unicycler_path:
            parser.print_help()
            sys.stderr.write('\n** Cannot find java in your path!  Is it installed?\n\n')
            sys.exit(1)

        if not os.path.exists(opts.jar_path):
            raise ValueError('Trimmomatic JAR file %s not found' % opts.jar_path)

        if opts.min_quality is None:
            opts.min_quality = TRIMMOMATIC_DEFAULT_QUALITY

        if opts.min_length is None:
            opts.min_length = TRIMMOMATIC_MIN_LENGTH

        fwd_filtered, rev_filtered = trimmomatic(fwd_fastq, rev_fastq, opts.min_length,
                                              opts.filter_window,
                                              opts.min_quality,
                                              opts.output_dir,
                                              opts.jar_path,
                                              verbose=opts.verbose)

    # Finally, load the filtered queries into an array
    queries, max_readlen = load_queries([fwd_filtered, rev_filtered])

    filter_endtime = time.time()
    filter_time = filter_endtime - start_time

    try:
        assembly_fasta = denovo_assembly(fwd_filtered, rev_filtered, opts.output_dir, verbose=opts.verbose)
        assembly_endtime = time.time()
        assembly_time = assembly_endtime - filter_endtime
    except OSError as ose:
        read_count_warnings(queries, opts)
        raise ose

    validate_file(assembly_fasta)
    final_fasta = os.path.join(opts.output_dir, FINAL_ASSEMBLY)
    fragments, assembly_len = post_process_assembly(assembly_fasta, final_fasta, assembly_name)

    if fragments == 1:
        print(time_string('Assembly contains %d contig of length %d.' % (fragments, assembly_len)))
    else:
        print(time_string('** Assembly contains %d fragments totalling %d bases, indicating possible contamination.' % (fragments, assembly_len)))
        total_time = time.time() - start_time

        read_count_warnings(queries, opts)

        print(time_string('Total time %.3f seconds; filtering %.3f seconds; assembly %.3f seconds.' % (total_time, filter_time, assembly_time)))
        sys.exit(0)

    if opts.fastp:
        # Note: mismatch analysis was not performed for fastp, so contamination likelihood is unavailable
        total_time = time.time() - start_time
        print(time_string('Total time %.3f seconds; filtering %.3f seconds; assembly %.3f seconds.' % (total_time, filter_time, assembly_time)))
        print(time_string('The final assembly has been written to %s' % final_fasta))
        sys.exit(0)

    # Run mismatch analysis using filtered subset:
    percent = compute_mismatch_percentage(assembly_fasta, queries, max_readlen, verbose=opts.verbose)
    contamination_likelihood = estimate_likelihood(percent, CONTAMINATED_CDF)

    end_time = time.time()
    mismatch_time = end_time - assembly_endtime
    total_time = end_time - start_time

    print(time_string('%.1f%% realignment mismatches (%d%% estimated likelihood of contamination)' % (percent, contamination_likelihood)))
    print(time_string('Total time %.3f seconds; filtering %.3f seconds; assembly %.3f seconds; realignment %.3g seconds.' % (total_time, filter_time, assembly_time, mismatch_time)))
    print(time_string('The final assembly has been written to %s' % final_fasta))