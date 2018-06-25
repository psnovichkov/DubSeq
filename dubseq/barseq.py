import os
import sys
import argparse
import logging
from .core.barcode import Barcode, BarcodeTag, BarcodeStat
from .core.fastq import FastqReader, FastqRecord, FastqFileStat
from .core import util


class Context:
    BARCODES_FNAME_SUFFIX = '.barcodes'
    BARCODE_STAT_FNAME_SUFFIX = '.bstat.tsv'
    LOG_FILE_NAME = 'barseq.log'

    output_dir = None

    barcode_tag = None
    fastq_source = None
    primer_position_shifts = None
    min_barcode_quality = None
    sim_ratio_threshold = None

    @staticmethod
    def build_context(args):

        pre_seaquence = 'CAGCGTACG'
        post_sequence = 'AGAGACC'
        pre_pos = 14
        position_shifts = [-2, -1, 0, 1, 2]

        if args.n25:
            # 11:14
            pre_pos = 11
            position_shifts = [0, 1, 2, 3]
        elif args.bs3:
            # 18:21
            pre_pos = 18
            position_shifts = [0, 1, 2, 3]
        else:
            pre_pos = args.pos1
            pre_seaquence = args.sequence1
            post_sequence = args.sequence2
            position_shifts = args.shift

        Context.barcode_tag = BarcodeTag(
            pre_seaquence, pre_pos,
            post_sequence, pre_pos + len(pre_seaquence) + 20)

        Context.primer_position_shifts = position_shifts

        Context.fastq_source = args.input
        Context.output_dir = args.output
        Context.min_barcode_quality = args.min_barcode_quality
        Context.sim_ratio_threshold = args.sim_ratio_threshold

    @staticmethod
    def barcodes_fname(fastq_file_name):
        base_file_name = os.path.basename(fastq_file_name)
        return os.path.join(Context.output_dir, base_file_name + Context.BARCODES_FNAME_SUFFIX)

    @staticmethod
    def bstat_fname(fastq_file_name):
        base_file_name = os.path.basename(fastq_file_name)
        return os.path.join(Context.output_dir, base_file_name + Context.BARCODE_STAT_FNAME_SUFFIX)

    @staticmethod
    def log_fname():
        return os.path.join(Context.output_dir, Context.LOG_FILE_NAME)


def parse_args():

    parser = argparse.ArgumentParser(
        description='''
        The narseq program processes the results of BarSeq assays to extract the 
        up barcodes and caluclate the number of reads supporting each barcode.

        The expected structure of a sequence read:

        ---[primer1]<barcode>[primer2]--
        
        Barseq program first extracts barcodes given the expected coordinates and sequences 
        of primers. The --primer-position-shifts paramter allows flexibility in the coordiates
        of primers.
        
        There are two predefined modes for newer multiplexing designs:
        1. n25 mode (--n25) means 11:14 nt before the pre-sequence, corresponding to a read with
	            2:5 Ns, GTCGACCTGCAGCGTACG, N20, AGAGACC
        2. bs3 mode (--bs3) eans 1:4 + 6 + 11 = 18:21 nt before the pre-sequence, corresponding to
	            1:4 Ns, index2, GTCGACCTGCAGCGTACG, N20, AGAGACC

        All extracted barcodes are filtered by the quality of their sequences controlled
        by --min-barcode-quality parameter. 
        
        To minimize the amount of erroneous barcodes cased by sequence errors introduced by PCR, 
        similar barcodes were identified for each extracted and mapped barcode. Two barcodes are
        considered to be similar if their sequnece is different by exactly one nucleotide. A given 
        barcode is considered to be true (recommended) barcode if all idnetified similar barcodes 
        are less frequent, and the ratio of frequencies of a given barcode and the most abundant
        similar barcode is greater than a threshold defined by --sim-ratio-threshold parameter.
        

        Bagseq program produces the following types of files:
        1. *.barcodes - file that has all extracted barcodes (for each source fastq file)
        2. *.bstat.tsv - file with statistics of extracted barcodes (for each source fastq file)


        Examples to run the barseq program:

        python -m dubseq.barseq -i /path/to/fastq/files -o /output/dir


        ''',
        formatter_class=util.RawDescriptionArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input',
                        dest='input',
                        help="""fastq file (can be gzipped) or directory with fastq files (can be gzipped).
                        If it is a directory, each file will be processed separately""",
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

    parser.add_argument('--n25',
                        dest='n25',
                        help=''' means 11:14 nt before the pre-sequence, corresponding to a read with
	                            2:5 Ns, GTCGACCTGCAGCGTACG, N20, AGAGACC ''',
                        action='store_true')

    parser.add_argument('--bs3',
                        dest='bs3',
                        help=''' means 1:4 + 6 + 11 = 18:21 nt before the pre-sequence, corresponding to
	                        1:4 Ns, index2, GTCGACCTGCAGCGTACG, N20, AGAGACC ''',
                        action='store_true')

    parser.add_argument('--primer1-pos',
                        dest='pos1',
                        help='Position of primer1',
                        default=14,
                        type=int
                        )

    parser.add_argument('--primer1-sequence',
                        dest='sequence1',
                        help='Sequence of primer1',
                        default='CAGCGTACG',
                        type=str
                        )

    parser.add_argument('--primer2-sequence',
                        dest='sequence2',
                        help='Sequence of primer2',
                        default='AGAGACC',
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
                        help='''The minimum ratio of frequencies of a given barcode and the most abundant similar 
                        barcodes to consider the given barcode to be true (not caused by sequence errors 
                        introduced by PCR)
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

    # process file(s)
    logging.info("Processing fastq files started")

    util.process_fastq_files(
        Context.fastq_source,
        process_fastq_file)

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


def process_fastq_file(fastq_fname):
    ''' Process fastq file.
        1. Extracts barcodes
        2. Collect stats on barcode frequency and similar barcodes (that differs by 1 nucleotide)
    '''
    print("Doing file:%s" % fastq_fname)

    # Hashtable to accumulate barcodes extracted from all fastq files
    # it is a hashtable: {barcode sequnece => BarcodeStat}
    barcode_stats = {}

    # To collect a general stats on the processed fastq file (will be logged)
    fastq_file_stat = FastqFileStat()

    print("\tExtracting barcodes...")
    extract_barcodes(fastq_fname, barcode_stats, fastq_file_stat)

    print("\tProcessing stat...")
    # Identify similar barcodes and collect number of reads supporintg similar barcodes
    BarcodeStat.find_similar_barcodes(barcode_stats)

    # Store barcoe stat
    BarcodeStat.save_barcode_stats(
        Context.bstat_fname(fastq_fname), barcode_stats)

    # Log the file stat
    logging.info("%s\t%s" % (fastq_fname, str(fastq_file_stat)))


def extract_barcodes(fastq_fname, barcode_stats, fastq_file_stat):
    '''
        It will:
        1. try to extract a barcode
        2. store the extracted barcode in barcodes file
        3. store the high quality barcode in barcode_stats for the downstream analysis
    '''

    # open a file to accumulate extracted barcodes
    barcodes_fp = open(Context.barcodes_fname(fastq_fname), 'w')
    barcodes_fp.write("seq_id\t%s\n" % Barcode.header())

    try:
        # open fastq file reader
        reader = FastqReader(fastq_fname)
        try:
            record = FastqRecord()
            while reader.next_record(record):

                # count the total amount of the processed reads
                fastq_file_stat.total_reads_inc()

                # try to extract a barcode from the sequence
                barcode = Context.barcode_tag.extract_barcode(
                    record, Context.primer_position_shifts)

                if barcode:
                    # count the reads with extracted barcodes
                    fastq_file_stat.barcode_extracted_reads_inc()

                    # store the extracted barcode
                    record_id = record.id.split(' ')[0]
                    barcodes_fp.write("%s\t%s\n" % (record_id, barcode))

                if barcode and barcode.min_quality >= Context.min_barcode_quality:

                    # register barcode in the barcode_stats
                    barcode_key = barcode.sequence
                    barcode_stat = barcode_stats.get(barcode_key)
                    if not barcode_stat:
                        barcode_stat = BarcodeStat()
                        barcode_stats[barcode_key] = barcode_stat

                    barcode_stat.reads_count_inc()
        finally:
            reader.close()
    finally:
        barcodes_fp.close()

    # Log the file stat
    logging.info("%s\t%s" % (fastq_fname, str(fastq_file_stat)))


if __name__ == '__main__':
    args = parse_args()
    check_args(args)
    Context.build_context(args)
    init_logger()

    main()
