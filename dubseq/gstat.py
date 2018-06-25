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
    GSCORE_BASE_FILE_NAME = 'gscore_base.tsv'

    LOG_FILE_NAME = 'bstat.log'

    barseq_layout_fname = None
    barseq_bstat_dir = None
    bpag_fname = None
    genes_gff_fname = None
    output_dir = None
    min_time0_read_count = None

    ridge_alpha = None
    lasso_alpha = None
    enet_alpha = None
    enet_l1_ratio = None
    gscore_varience_alpha = None

    @staticmethod
    def build_context(args):
        Context.barseq_layout_fname = args.barseq_layout_fname
        Context.barseq_bstat_dir = args.input
        Context.bpag_fname = args.bpag_fname
        Context.output_dir = args.output
        Context.min_time0_read_count = args.min_time0_read_count
        Context.genes_gff_fname = args.genes_gff_fname
        Context.ridge_alpha = args.ridge_alpha
        Context.lasso_alpha = args.lasso_alpha
        Context.enet_alpha = args.enet_alpha
        Context.enet_l1_ratio = args.enet_l1_ratio
        Context.gscore_varience_alpha = args.gscore_varience_alpha

    @staticmethod
    def gscore_base_fname():
        return os.path.join(Context.output_dir, Context.GSCORE_BASE_FILE_NAME)

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
            The gstat program ...


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

    parser.add_argument('-g', '--genes-gff-fname',
                        dest='genes_gff_fname',
                        help='path to a gff file for the genome used to build dubseq library',
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

    parser.add_argument('--ridge_alpha',
                        dest='ridge_alpha',
                        help='''Regularization parameter alpha for the Ridge regression defining
                        the amount of regularization in the Ridge objective function:
                        ||Ax-y||^2_2 + alpha * ||x||^2_2 ''',
                        default=1.0,
                        type=float
                        )

    parser.add_argument('--lasso_alpha',
                        dest='lasso_alpha',
                        help='''Regularization parameter alpha for the Lasso regression defining
                        the amount of regularization in the Lasso objective function:
                        ||Ax-y||^2_2 + alpha * ||x||_1 ''',
                        default=3.35,
                        type=float
                        )

    parser.add_argument('--enet_alpha',
                        dest='enet_alpha',
                        help='''Regularization parameter alpha for the Elastic Net regression defining
                        the amount of regularization in the Elastic Net objective function:
                        ||Ax-y||^2_2 + alpha * i1_ratio * ||x||_1 + 0.5 * alpha * (1-r1_ratio) * ||x||^2_2  ''',
                        default=3.62,
                        type=float
                        )

    parser.add_argument('--enet_i1_ratio',
                        dest='enet_l1_ratio',
                        help='''Regularization parameter l1_ratio for the Elastic Net regression defining
                        the amount of regularization in the Elastic Net objective function:
                        ||Ax-y||^2_2 + alpha * i1_ratio * ||x||_1 + 0.5 * alpha * (1-r1_ratio) * ||x||^2_2  ''',
                        default=0.7,
                        type=float
                        )

    parser.add_argument('--gscore_varience_alpha',
                        dest='gscore_varience_alpha',
                        help=''' Correction to calculate a moderated version of a gene score variance
                         Var_moderated = Var + alpha   
                        ''',
                        default=0.02,
                        type=float
                        )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


def check_args(args):
    pass


def main():
    Fitness.MIN_TIME0_READ_COUNT = Context.min_time0_read_count
    Fitness.RIDGE_PARAM_ALPHA = Context.ridge_alpha
    Fitness.LASSO_PARAM_ALPHA = Context.lasso_alpha
    Fitness.ELASTIC_NET_PARAM_ALPHA = Context.enet_alpha
    Fitness.ELASTIC_NET_PARAM_L1_RATIO = Context.enet_l1_ratio
    Fitness.GSCORE_VARIENCE_ALPHA = Context.gscore_varience_alpha

    barseq_layout = BarseqLayout(Context.barseq_layout_fname)
    barseq_layout.save(Context.barseq_layout_out_fname())

    Fitness.init(barseq_layout, Context.barseq_bstat_dir,
                 Context.bpag_fname, Context.genes_gff_fname)
    Fitness.save_fscore_base(Context.fscore_base_fname())
    Fitness.save_gscore_base(Context.gscore_base_fname())

    # First calculate all basic stat for all conditions (including qvalue for an individual condition)
    gstats_set = []
    for index, item in enumerate(barseq_layout.all_items):
        print('Doing %s' % item.itnum)

        (gstats, gscore_var_eff, pi0) = Fitness.build_gstat(
            index, Fitness.SCORE_TYPE_C_NNLS)

        print('gscore_var_eff=%s, pi0=%s' % (gscore_var_eff, pi0))
        gstats_set.append(gstats)

    # calculate qvalues based on the whole set
    (gstats_set, pi0) = Fitness.calculate_set_qvalues(gstats_set)
    print('set pi0=%s' % pi0)

    # Save gstats
    for index, item in enumerate(barseq_layout.all_items):
        print('Saving gstat for %s' % item.itnum)

        gstat_fname = os.path.join(
            Context.output_dir, item.itnum + '.gstat.tsv')
        Fitness.save_gstat(gstat_fname, gstats_set[index])


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
