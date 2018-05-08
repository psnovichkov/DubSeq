import os
import sys
import argparse
import logging
from .core import util


class Context:
    BPAG_FILE_NAME = 'bpag.tsv'
    LOG_FILE_NAME = 'bpag.log'

    output_dir = None
    bpseq_fname = None
    bagseq_up_fname = None
    bagseq_dn_fname = None

    min_region_len = None
    max_region_len = None

    @staticmethod
    def build_context(args):
        Context.output_dir = args.output

        Context.bpseq_fname = args.bpseq_fname
        Context.bagseq_up_fname = args.bagseq_up_fname
        Context.bagseq_dn_fname = args.bagseq_dn_fname
        Context.min_region_len = args.min_region_len
        Context.max_region_len = args.max_region_len

    @staticmethod
    def bpag_fname():
        return os.path.join(Context.output_dir, Context.BPAG_FILE_NAME)

    @staticmethod
    def log_fname():
        return os.path.join(Context.output_dir, Context.LOG_FILE_NAME)


def parse_args():

    parser = argparse.ArgumentParser(
        description='''
        The bpag program combines the barcode pairs produced by bpseq program and
        up and down barcodes mapped to genomic DNA by bageseq program to generate
        a master list of barcode pairs associated with genomic regions.

        Bpag program iterates over all barcode pairs produced by bpseq program and
        checks the presence of the up and down barcodes mapped to a genomic location.
        A barcode pair is considered to be valid if both up and down barcodes are 
        mapped to the same contig and strands of hits are opposite. 

        The valid barcode pair is consideed to be recommended if the lengh of genomic 
        region defined by mapped locations of up and down barcodes is more than min
        and less than max expected lenght defined by --min-region-len and --max-region-len
        parameters.

        Bpag produces bpag.tsv file with barcode pair, cooridnates of the correspodning 
        genomic region, statsitics and recommnedation flag.

        Examples to run the bpag program:

        python -m dubseq.bpag -p /path/to/bpseq.tsv 
               -u  /path/to/bagseq_up.tsv  -d /path/to/bagseq_dn.tsv
               -o /output/dir                                    
        ''',
        formatter_class=util.RawDescriptionArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--bpseq-fname',
                        dest='bpseq_fname',
                        help='file name with barcode pairs produced by bpseq program',
                        type=str,
                        required=True
                        )

    parser.add_argument('-u', '--bagseq-up-fname',
                        dest='bagseq_up_fname',
                        help='file name with up barcodes mapped to a genomic location by bagseq program',
                        type=str,
                        required=True
                        )

    parser.add_argument('-d', '--bagseq-dn-fname',
                        dest='bagseq_dn_fname',
                        help='file name with down barcodes mapped to a genomic location by bagseq program',
                        type=str,
                        required=True
                        )

    parser.add_argument('-o', '--output',
                        dest='output',
                        help='output directory',
                        type=str,
                        required=True
                        )

    parser.add_argument('-n', '--min-region-len',
                        dest='min_region_len',
                        help='The minmal expected lenght of a genomic fragment',
                        default=100,
                        type=int
                        )

    parser.add_argument('-x', '--max-region-len',
                        dest='max_region_len',
                        help='The maximal expected lenght of a genomic fragment',
                        default=6000,
                        type=int
                        )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


def check_args(args):
    pass


class _BLocation:
    __slots__ = [
        'read_count', 'contig_id', 'pos', 'strand'
    ]

    @staticmethod
    def header(sep='\t', prefix=''):
        return sep.join(prefix + x for x in [
            'read_count', 'contig_id', 'pos', 'strand'
        ])

    def __init__(self, read_count, contig_id, pos, strand):
        self.read_count = read_count
        self.contig_id = contig_id
        self.pos = pos
        self.strand = strand

    def __str__(self):
        return '\t'.join(str(x) for x in [
            self.read_count,
            self.contig_id,
            self.pos,
            self.strand
        ]
        )


def main():

    logging.info("Processing fastq files started")

    up_barcodes_locations = load_barcode_locations(
        Context.bagseq_up_fname, reverse_complement=True)
    dn_barcodes_locations = load_barcode_locations(
        Context.bagseq_dn_fname, reverse_complement=False)

    BPAIR_RECOMMENDED_INDEX = 2
    UP_SIM_RECOMMENDED_INDEX = 5
    UP_CHIM_RECOMMENDED_INDEX = 9
    DN_SIM_RECOMMENDED_INDEX = 14
    DN_CHIM_RECOMMENDED_INDEX = 18

    with open(Context.bpag_fname(), 'w') as fout:
        fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
            'barcode_up',
            'barcode_dn',
            'bpair_count',
            'pos_from',
            'pos_end',
            'region_len',
            'recommended',
            _BLocation.header(prefix='up_'),
            _BLocation.header(prefix='dn_')
        ))
        with open(Context.bpseq_fname, 'r') as fin:
            fin.readline()
            for line in fin:
                vals = line.split('\t')
                up_barcode = vals[0]
                dn_barcode = vals[1]
                bpair_count = int(vals[3])

                good_pair = True
                for index in [BPAIR_RECOMMENDED_INDEX,
                              UP_SIM_RECOMMENDED_INDEX, UP_CHIM_RECOMMENDED_INDEX,
                              DN_SIM_RECOMMENDED_INDEX, DN_CHIM_RECOMMENDED_INDEX]:
                    if vals[index] == '-':
                        good_pair = False

                if good_pair:
                    process_pair(fout, up_barcode, dn_barcode, bpair_count,
                                 up_barcodes_locations, dn_barcodes_locations)

    logging.info("Done!")


def load_barcode_locations(bagseq_fname, reverse_complement=False):
    barcodes_locations = {}
    with open(bagseq_fname, 'r') as f:
        f.readline()
        for line in f:
            vals = line.split('\t')
            sim_recommended = vals[2]
            loc_recommended = vals[6]

            if sim_recommended == '+' and loc_recommended == '+':
                barcode = vals[0]
                if reverse_complement:
                    barcode = util.reverse_complement(barcode)

                barcodes_locations[barcode] = _BLocation(
                    int(vals[8]),
                    vals[10],
                    int(vals[11]),
                    vals[12]
                )

    return barcodes_locations


def process_pair(fout, up_barcode, dn_barcode, bpair_count,
                 up_barcodes_locations, dn_barcodes_locations):

    up_location = up_barcodes_locations.get(up_barcode)
    dn_location = dn_barcodes_locations.get(dn_barcode)

    if up_location and dn_location:

        if up_location.contig_id == dn_location.contig_id \
                and up_location.strand != dn_location.strand:

            pos_from = 0
            pos_end = 0

            if up_location.strand == '-':
                pos_from = dn_location.pos
                pos_end = up_location.pos
            else:
                pos_from = up_location.pos
                pos_end = dn_location.pos

            region_len = pos_end - pos_from + 1
            recommended = (
                region_len >= Context.min_region_len and region_len <= Context.max_region_len)

            fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                up_barcode, dn_barcode, bpair_count,
                pos_from, pos_end, region_len,
                '+' if recommended else '-',
                str(up_location), str(dn_location)
            ))


def init_logger():
    with open(Context.log_fname(), 'w') as f:
        f.write("Parameters:\n")
        for arg, value in vars(args).items():
            f.write("\t%s=%s\n" % (arg, value))

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
