import os
import sys
import argparse
import logging
import struct
from .core.barcode import Barcode, BarcodeTag, EMPTY_BARCODE, BarcodeStat, PairedBarcodeStat
from .core.fastq import FastqReader, FastqRecord, FastqFileStat
from .core import util


class Context:
    BARCODE_SEPARATOR = '='
    BARCODE_PAIR_STAT_FILE_NAME = 'bpseq.tsv'
    BARCODE_UP_STAT_FILE_NAME = 'up.bstat.tsv'
    BARCODE_DN_STAT_FILE_NAME = 'dn.bstat.tsv'
    BARCODES_FILE_SUFFIX = '.barcodes'
    LOG_FILE_NAME = 'bpseq.log'

    output_dir = None
    up_tag = None
    dn_tag = None
    primer_position_shifts = None
    min_barcode_quality = None
    file_stat = None
    sim_ratio_threshold = None
    chim_ratio_threshold = None

    @staticmethod
    def build_context(args):
        Context.up_tag = BarcodeTag(args.up_primer1_sequence, args.up_primer1_pos,
                                    args.up_primer2_sequence, args.up_primer2_pos)

        Context.dn_tag = BarcodeTag(args.dn_primer1_sequence, args.dn_primer1_pos,
                                    args.dn_primer2_sequence, args.dn_primer2_pos)
        Context.fastq_source = args.input
        Context.output_dir = args.output
        Context.primer_position_shifts = args.shift
        Context.min_barcode_quality = args.min_barcode_quality
        Context.sim_ratio_threshold = args.sim_ratio_threshold
        Context.chim_ratio_threshold = args.chim_ratio_threshold

    @staticmethod
    def log_fname():
        return os.path.join(Context.output_dir, Context.LOG_FILE_NAME)

    @staticmethod
    def barcode_pair_stat_fname():
        return os.path.join(Context.output_dir, Context.BARCODE_PAIR_STAT_FILE_NAME)

    @staticmethod
    def barcode_up_stat_fname():
        return os.path.join(Context.output_dir, Context.BARCODE_UP_STAT_FILE_NAME)

    @staticmethod
    def barcode_dn_stat_fname():
        return os.path.join(Context.output_dir, Context.BARCODE_DN_STAT_FILE_NAME)

    @staticmethod
    def barcodes_fname(fastq_fname):
        base_file_name = os.path.basename(fastq_fname)
        return os.path.join(Context.output_dir, base_file_name + Context.BARCODES_FILE_SUFFIX)


def parse_args():

    parser = argparse.ArgumentParser(
        description='''
        The bpseq program processes the results of BPSeq PCR assay to extract the pairs of 
        aassociated barcodes. Each read is expected to have two associated 20-mer barocode 
        tags: up and down tag.

        The expected structure of a sequence read:

        ---[dn_primer1]<dn_barcode>[dn_primer2]----[up_primer1]<up_barcode>[dn_primer2]---
        ...[........... down tag  ............]----[........... up tag  ..............]...
        
        Bpseq program first extracts up and down barcodes given the expected coordinates and 
        sequences of primers. The --primer-position-shifts paramter allows flexibility in the 
        coordiates of primers.
        
        All extracted barcodes are filtered by the quality of their sequences controlled
        by --min-barcode-quality parameter. 
        
        Bpseq program implements two additional types of filters to minimize the number of 
        erroneous pairs of barcodes: i) sequence errors introduced by PCR, ii) chimeric 
        constructs created by PCR

        To minimize the amount of erroneous barcodes caused by sequence errors, similar 
        barcodes were identified for each extracted barcode (up and down barcodes were processed 
        separately). Two barcodes are considered to be similar if their sequnece is different 
        by exactly one nucleotide. A given barcode is considered to be true (recommended) 
        barcode if all idnetified similar barcodes are less frequent, and the ratio of frequencies 
        of a given barcode and the most abundant similar barcode is greater than a threshold
        defined by --sim-ratio-threshold parameter.
        
        To minimize the amount of erroneous barcode pairs caused by chimeras, barcode pairs 
        associated with each up (down) barcode were collected. The presence of the same  
        barcode in multiple barcode pairs is a sign of chimeric sequence. To distinguish the 
        true barcode pair from the chimeric one, the frequency of a barcode pair was compared  
        with frequncy of all associated barcode pairs. A given barcode pair is considered to  
        be non-chimeric (recommended) if all associated barcode pairs are less frequent, and 
        the ratio of frequencies of a given barcode pair and the most abundant associated
        barcode pair is greater than a theshold defined by --chim-ratio-threshold parameter.
        
        Overall, a given barcode pair is considered to be true barcode pair if up and down 
        barcodes have high sequence quality, both up and down barcodes are recommended as 
        non-erroneous, and barcode pair is recommended as non-chimeric.

        Bpseq program produces four types of files:
        1. *.barcodes - file that has all extracted barcode pairs (for each source fastq file)
        2. up.bstat.tsv - statistics for the up barcodes
        3. dn.bstat.tsv - statistics for the down barcodes
        4. bpseq.tsv - the main output file with barcode paris, statistics  and recomendation 
        flags

        Examples to run the bpseq program:

        python -m dubseq.bpseq -i /path/to/fastq/files -o /output/dir

        ''',
        formatter_class=util.RawDescriptionArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input',
                        dest='input',
                        help='''fastq file (can be gzipped) or directory with fastq files 
                        (can be gzipped).''',
                        type=str,
                        required=True
                        )

    parser.add_argument('-o', '--output',
                        dest='output',
                        help='output directory',
                        type=str,
                        required=True
                        )

    parser.add_argument('-q', '--min-barcode-quality',
                        dest='min_barcode_quality',
                        help='The minimal quality of the barcode sequence',
                        default=20,
                        type=int
                        )

    parser.add_argument('--up-primer1-pos',
                        dest='up_primer1_pos',
                        help='Position of primer1 of up tag',
                        default=90,
                        type=int
                        )

    parser.add_argument('--up-primer1-sequence',
                        dest='up_primer1_sequence',
                        help='Sequence of primer1 of up tag',
                        default='CAGCGTACG',
                        type=str
                        )

    parser.add_argument('--up-primer2-pos',
                        dest='up_primer2_pos',
                        help='Position of primer2 of up tag',
                        default=119,
                        type=int
                        )

    parser.add_argument('--up-primer2-sequence',
                        dest='up_primer2_sequence',
                        help='Sequence of primer2 of up tag',
                        default='AGAGACCTC',
                        type=str
                        )

    parser.add_argument('--dn-primer1-pos',
                        dest='dn_primer1_pos',
                        help='Position of primer1 of down tag',
                        default=13,
                        type=int
                        )

    parser.add_argument('--dn-primer1-sequence',
                        dest='dn_primer1_sequence',
                        help='Sequence of primer1 of down tag',
                        default='GTCTCGTAG',
                        type=str
                        )

    parser.add_argument('--dn-primer2-pos',
                        dest='dn_primer2_pos',
                        help='Position of primer2 of down tag',
                        default=42,
                        type=int
                        )

    parser.add_argument('--dn-primer2-sequence',
                        dest='dn_primer2_sequence',
                        help='Sequence of primer2 of down tag',
                        default='CGATGAAT',
                        type=str
                        )

    parser.add_argument('-s', '--primer-position-shifts',
                        dest='shift',
                        help='Alowed shifts of the primer position',
                        default=[0, -1, 1, -2, 2],
                        type=int,
                        nargs='+'
                        )

    parser.add_argument('-m', '--sim-ratio-threshold',
                        dest='sim_ratio_threshold',
                        help='''The minimum ratio of frequencies of a given and the most abundant similar 
                        barcodes to consider the given barcode to be true (not caused by sequence errors 
                        introduced by PCR)
                        ''',
                        default=2,
                        type=int
                        )

    parser.add_argument('-c', '--chim-ratio-threshold',
                        dest='chim_ratio_threshold',
                        help='''The minimum ratio of frequencies of a given barcode pair and the most abundant 
                        associated barcode pair to consider the given barcode pair to be non-chimeric
                        ''',
                        default=2,
                        type=int
                        )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


def check_args(args):
    pass


def main():
    BarcodeStat.SIM_RATIO_THRESHOLD = Context.sim_ratio_threshold
    PairedBarcodeStat.CHIM_RATIO_THRESHOLD = Context.chim_ratio_threshold

    up_barcodes = {}
    dn_barcodes = {}
    barcodes12 = {}

    logging.info("Processing fastq files started")
    util.process_fastq_files(Context.fastq_source,
                             process_fastq_file, barcodes12, up_barcodes, dn_barcodes)

    print('Analyze similar barcodes: up tags')
    BarcodeStat.find_similar_barcodes(up_barcodes)
    update_pair_counts(up_barcodes, barcodes12, 0)
    PairedBarcodeStat.save_barcode_stats(
        Context.barcode_up_stat_fname(), up_barcodes)
    # process_barcode_stat(up_barcodes)

    print('Analyze similar barcodes: down tags')
    BarcodeStat.find_similar_barcodes(dn_barcodes)
    update_pair_counts(dn_barcodes, barcodes12, 1)
    PairedBarcodeStat.save_barcode_stats(
        Context.barcode_dn_stat_fname(), dn_barcodes)
    # process_barcode_stat(dn_barcodes)

    print('Export results')
    save_barcode_pair_stat(barcodes12, up_barcodes, dn_barcodes)

    logging.info("Done!")


def init_logger():
    with open(Context.log_fname(), 'w') as f:
        f.write("Parameters:\n")
        for arg, value in vars(args).items():
            f.write("\t%s=%s\n" % (arg, value))
        f.write("Report columns:\n")
        f.write("\t%s\n" % FastqFileStat.header(sep='\n\t'))
        f.write("\n\n")

    logging.basicConfig(
        filename=Context.log_fname(),
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p")


def process_fastq_file(fastq_fname, barcodes12, up_barcodes, dn_barcodes):
    ''' Process fastq file.
        1. Extracts barcode pairs
        2. Collect stats on barcode frequency and similar barcodes (that differs by 1 nucleotide)
    '''
    print("Doing file:%s" % fastq_fname)

    fastq_file_stat = FastqFileStat()
    print("\tExtracting barcodes...")
    extract_barcodes(fastq_fname, barcodes12, up_barcodes,
                     dn_barcodes, fastq_file_stat)

    # Log the file stat
    logging.info("%s\t%s" % (fastq_fname, str(fastq_file_stat)))


def extract_barcodes(fastq_fname, barcodes12, up_barcodes, dn_barcodes, fastq_file_stat):
    ''' Extracts barcodes from fastq file based on the barcode_tag and primer_position_shifts.
        The list of barcodes is stored in the output_dir.
        The barcodes with good qaulity are collected in the barcode2stat for further processing.
    '''

    # open a file to store the extracted barcodes
    with open(Context.barcodes_fname(fastq_fname), 'w') as f:
        # write a header
        f.write("seq_id\t%s\t%s\n" % (
            Barcode.header(prefix='up_'),
            Barcode.header(prefix='dn_'),
        ))

        try:
            # Process each record from the fastq file
            reader = FastqReader(fastq_fname)
            record = FastqRecord()
            i = 0
            while reader.next_record(record):
                i += 1
                if i % 100000 == 0:
                    print('.', end='', flush=True)
                if i % 1000000 == 0:
                    print(i)

                fastq_file_stat.total_reads_inc()

                # try to extract a barcode from the sequence
                up_barcode = Context.up_tag.extract_barcode(
                    record, Context.primer_position_shifts)

                dn_barcode = Context.dn_tag.extract_barcode(
                    record, Context.primer_position_shifts)

                if up_barcode or dn_barcode:
                    # store the extracted barcode pair
                    record_id = record.id.split(' ')[0]
                    f.write("%s\t%s\t%s\n" %
                            (record_id,
                             up_barcode if up_barcode else EMPTY_BARCODE,
                             dn_barcode if dn_barcode else EMPTY_BARCODE))

                    # both barcodes should be present for the downstresam analysis
                    if not (up_barcode and dn_barcode):
                        continue

                    # store high quality barcode pairs in barcode2stat dictionary for the downstream analysis
                    if up_barcode.min_quality >= Context.min_barcode_quality \
                            and dn_barcode.min_quality >= Context.min_barcode_quality:

                        fastq_file_stat.barcode_extracted_reads_inc()

                        up_key = up_barcode.sequence
                        dn_key = dn_barcode.sequence
                        bpair_key = up_key + Context.BARCODE_SEPARATOR + dn_key
                        bpair_reads_count = barcodes12.get(bpair_key)

                        up_stat = None
                        dn_stat = None
                        if bpair_reads_count:
                            barcodes12[bpair_key] = bpair_reads_count + 1
                            up_stat = up_barcodes[up_key]
                            dn_stat = dn_barcodes[dn_key]
                        else:
                            barcodes12[bpair_key] = 1

                            up_stat = up_barcodes.get(up_key)
                            if not up_stat:
                                up_stat = PairedBarcodeStat()
                                up_barcodes[up_key] = up_stat

                            dn_stat = dn_barcodes.get(dn_key)
                            if not dn_stat:
                                dn_stat = PairedBarcodeStat()
                                dn_barcodes[dn_key] = dn_stat

                        up_stat.reads_count_inc()
                        dn_stat.reads_count_inc()

        finally:
            print('')
            reader.close()


def update_pair_counts(barcode_stat, barcodes12, barcode_index):
    for bpair, bpair_count in barcodes12.items():
        barcode = bpair.split(Context.BARCODE_SEPARATOR)[barcode_index]
        barcode_stat[barcode].add_pair_reads_count(bpair_count)


def save_barcode_pair_stat(barcodes12, up_barcodes, dn_barcodes):

    # store barcode stat
    with open(Context.barcode_pair_stat_fname(), 'w') as f:
        # write a header
        f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
            'up_barcode', 'dn_barcode', 'barcode_pair_recommended', 'reads_count',
            PairedBarcodeStat.header(prefix='up_'),
            PairedBarcodeStat.header(prefix='dn_')
        ))

        i = 0
        count = len(barcodes12)
        for bpair, bpair_reads_count in barcodes12.items():
            i += 1
            if i % 100000 == 0:
                print('.', end='', flush=True)
            if i % 1000000 == 0:
                print('%s of %s' % (i, count))

            up_barcode, dn_barcode = bpair.split(Context.BARCODE_SEPARATOR)
            up_stat = up_barcodes[up_barcode]
            dn_stat = dn_barcodes[dn_barcode]

            bpair_recommended = True
            if not up_stat.sim_recommended():
                bpair_recommended = False
            if not dn_stat.sim_recommended():
                bpair_recommended = False

            if not up_stat.chim_recommended():
                bpair_recommended = False
            if not dn_stat.chim_recommended():
                bpair_recommended = False

            if bpair_reads_count < up_stat.pair_reads_count_max:
                bpair_recommended = False

            if bpair_reads_count < dn_stat.pair_reads_count_max:
                bpair_recommended = False

            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                up_barcode, dn_barcode,
                '+' if bpair_recommended else '-',
                bpair_reads_count,
                up_stat, dn_stat
            ))


if __name__ == '__main__':
    args = parse_args()
    check_args(args)
    Context.build_context(args)
    init_logger()

    main()
