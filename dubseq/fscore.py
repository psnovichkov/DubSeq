import os
import sys
import argparse
import logging
from .core.barcode import Barcode, BarcodeTag, BarcodeStat
from .core.fastq import FastqReader, FastqRecord, FastqFileStat
from .core import util
from .core.fitness import BarseqLayout, Fitness


class Context:
    BARCODE_STAT_FNAME_SUFFIX = '.bstat.tsv'
    FSCORE_BASE_FILE_NAME = 'fscore_base.tsv'
    BARSEQ_LAYOUT_OUT_FILE_NAME = 'barseq_layout.tsv'
    LOG_FILE_NAME = 'fscore.log'

    barseq_layout_fname = None
    barseq_bstat_dir = None
    bpag_fname = None
    output_dir = None
    min_time0_read_count = None

    @staticmethod
    def build_context(args):
        Context.barseq_layout_fname = args.barseq_layout_fname
        Context.barseq_bstat_dir = args.input
        Context.bpag_fname = args.bpag_fname
        Context.output_dir = args.output
        Context.min_time0_read_count = args.min_time0_read_count

    # @staticmethod
    # def barcodes_fname(fastq_file_name):
    #     base_file_name = os.path.basename(fastq_file_name)
    #     return os.path.join(Context.output_dir, base_file_name + Context.BARCODES_FNAME_SUFFIX)

    @staticmethod
    def fscore_base_fname():
        return os.path.join(Context.output_dir, Context.FSCORE_BASE_FILE_NAME)

    @staticmethod
    def barseq_layout_out_fname():
        return os.path.join(Context.output_dir, Context.BARSEQ_LAYOUT_OUT_FILE_NAME)

    @staticmethod
    def log_fname():
        return os.path.join(Context.output_dir, Context.LOG_FILE_NAME)


def parse_args():

    parser = argparse.ArgumentParser(
        description='''
        The fscore program calculates the fitness score of DubSeq fragments in a given sample based on the 
        difference of read counts supporting a given fragment in the stress and time-zero conditions. 
        
        The reference set of the DubSeq fragments, specified by --bpag-fname parameter, can be obtained by the bpag 
        program that combines the results of the BPSeq and BAGSeq assays. Each fragment in this set is associated with 
        a barcode which frequency under a given condition can be later estiamted in a BarSeq assay. The directory 
        with restuls of the BarSeq assays procssed by the barseq program can be indicated with --input parameter.        
        The layout of the BarSeq experiment set is necessary to identiy time-zero samples and should be speicified with
        --barseq-layout-fname parameter.   
        
        In order to calculate the fragment fitness scores (fscores), the program first identifies a subset of fragments that 
        are well represented in time-zero samples. A given fragment is considered to be valid if at least in one time-zero 
        sample it is supported by more reads than a given threshold controlled by the --min-time0-read-count parameter.
        
        
        The fitness score of a valid fragment is caluclated as a normalized log2 ratio of read counts supporting the associated 
        barcode in the treatment sample and in time-zero samples. The fitness scores in a given sample are normalized so that 
        the median is set to zero.

        Bpseq program produces four types of files:
        1. ITNNN.fscore.tsv - fitness scores of fragmetns for a given sample (for each sample defined by IT number)
        2. fscore_base.tsv - subset of fragments considered to be valid for a given experiment set
        3. barseq_layout.tsv - copy of the BarSeq layout used to define the time-zero samples
        ''',
        formatter_class=util.RawDescriptionArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input',
                        dest='input',
                        help='path to the directory with bstat files produced by the barseq program',
                        type=str,
                        required=True
                        )

    parser.add_argument('-l', '--barseq-layout-fname',
                        dest='barseq_layout_fname',
                        help='path to a file with layout of barseq experiments',
                        type=str,
                        required=True
                        )

    parser.add_argument('-p', '--bpag-fname',
                        dest='bpag_fname',
                        help='path to a file with barcode pairs mapped to a genome using bpag program',
                        type=str,
                        required=True
                        )

    parser.add_argument('-o', '--output',
                        dest='output',
                        help='output directory',
                        type=str,
                        required=True
                        )

    parser.add_argument('-t', '--min-time0-read-count',
                        dest='min_time0_read_count',
                        help='The minimal required number of reads supporting a barcode in time zero',
                        default=10,
                        type=int
                        )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


def check_args(args):
    pass


def main():
    Fitness.MIN_TIME0_READ_COUNT = Context.min_time0_read_count

    barseq_layout = BarseqLayout(Context.barseq_layout_fname)
    barseq_layout.save(Context.barseq_layout_out_fname())
    Fitness.init(barseq_layout, Context.barseq_bstat_dir, Context.bpag_fname)
    Fitness.save_fscore_base(Context.fscore_base_fname())

    for index, item in enumerate(barseq_layout.all_items):
        print('Doing %s' % item.itnum)
        ss = Fitness.get_sample(index)
        ts = Fitness.get_tzero_sample()
        fs = Fitness.build_fscores(ss, ts)

        fscore_fname = os.path.join(
            Context.output_dir, item.itnum + '.fscore.tsv')

        print('\tsaving scores to %s' % fscore_fname)
        Fitness.save_fscores(fscore_fname, fs, ss, ts)

    # # process file(s)
    # logging.info("Processing fastq files started")

    # util.process_fastq_files(
    #     Context.fastq_source,
    #     process_fastq_file)

    # logging.info("Done!")


def init_logger():
    with open(Context.log_fname(), 'w') as f:
        f.write("Parameters:\n")
        for arg, value in vars(args).items():
            f.write("\t%s=%s\n" % (arg, value))
        # f.write("Report columns:\n")
        # f.write("\t%s\n" % FastqFileStat.header(sep='\n\t'))
        # f.write("\n\n")

    logging.basicConfig(
        filename=Context.log_fname(),
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p")


if __name__ == '__main__':
    args = parse_args()
    check_args(args)
    Context.build_context(args)
    init_logger()

    main()
