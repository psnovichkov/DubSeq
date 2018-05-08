import os
import sys
import argparse
import logging
import subprocess
from .core.barcode import Barcode, BarcodeTag, EMPTY_BARCODE, BarcodeStat, BarcodeLocation, BarcodeHits
from .core.fastq import FastqReader, FastqRecord, FastqFileStat
from .core.blat import BlatReader, BlatRecord
from .core import util


class Context:
    BLAT_CMD = '/usr2/people/pnovichkov/tools/blatSuite/blat'
    LOG_FILE_NAME = 'bagseq.log'

    FASTA_FNAME_SUFFIX = '.fna'
    BLAT_FNAME_SUFFIX = '.psl'
    BARCODES_FNAME_SUFFIX = '.barcodes'
    BHIT_FNAME_SUFFIX = '.bhit.tsv'

    BAGSEQ_OUTPUT_FNAME_PATTERN = 'bagseq_%s.tsv'
    UP_TAG_FNAME_PREFIX = 'up'
    DOWN_TAG_FNAME_PREFIX = 'dn'

    up_tag_fastq_source = None
    up_tag = None

    down_tag_fastq_source = None
    down_tag = None

    output_dir = None
    genome_fna_fname = None
    min_barcode_quality = None
    primer_position_shifts = None
    min_genomic_sequence_length = None
    max_blat_gap_bases = None
    min_blat_block_size = None
    sim_ratio_threshold = None
    loc_ratio_threshold = None

    @staticmethod
    def build_context(args):
        Context.output_dir = args.output

        Context.up_tag_fastq_source = args.input_up
        Context.up_tag = BarcodeTag(args.up_sequence1, args.up_pos1,
                                    args.up_sequence2, args.up_pos2)

        Context.down_tag_fastq_source = args.input_down
        Context.down_tag = BarcodeTag(args.dn_sequence1, args.dn_pos1,
                                      args.dn_sequence2, args.dn_pos2)

        Context.genome_fna_fname = args.input_genome
        Context.primer_position_shifts = args.shift
        Context.min_barcode_quality = args.min_barcode_quality
        Context.min_genomic_sequence_length = args.min_sequence_length
        Context.max_blat_gap_bases = args.max_blat_gap_bases
        Context.min_blat_block_size = args.min_blat_block_size
        Context.sim_ratio_threshold = args.sim_ratio_threshold
        Context.loc_ratio_threshold = args.loc_ratio_threshold

    @staticmethod
    def fasta_fname(fname_prefix):
        return os.path.join(Context.output_dir, fname_prefix + Context.FASTA_FNAME_SUFFIX)

    @staticmethod
    def log_fname():
        return os.path.join(Context.output_dir, Context.LOG_FILE_NAME)

    @staticmethod
    def bhit_fname(fname_prefix):
        return os.path.join(Context.output_dir, fname_prefix + Context.BHIT_FNAME_SUFFIX)

    @staticmethod
    def bagseq_output_fname(fname_prefix):
        return os.path.join(Context.output_dir, Context.BAGSEQ_OUTPUT_FNAME_PATTERN % fname_prefix)

    @staticmethod
    def blat_fname(fname_prefix):
        return os.path.join(Context.output_dir, fname_prefix + Context.BLAT_FNAME_SUFFIX)

    @staticmethod
    def barcodes_fname(fastq_file_name):
        base_file_name = os.path.basename(fastq_file_name)
        return os.path.join(Context.output_dir, base_file_name + Context.BARCODES_FNAME_SUFFIX)


def parse_args():

    parser = argparse.ArgumentParser(
        description='''
        The bagseq program processes the results of BAGSeq assay to extract the 
        up and down barcodes and their associated genomic sequences, and then map the sequences 
        to a genome. The ultimate goal of the program is to build a reliable association
        between a given up (down) barcode and the location of the corresponding genomic fragment 
        on the genome.

        Up and down barcodes are processed separately 

        The expected structure of a sequence read:

        ---[primer1]<barcode>[primer2]--[genomic_fragment]
        
        Bagseq program first extracts barcodes and candidate genomic fragment given the 
        expected coordinates and sequences of primers. The --primer-position-shifts paramter 
        allows flexibility in the coordiates of primers.
        
        All extracted barcodes are filtered by the quality of their sequences controlled
        by --min-barcode-quality parameter. 
        
        The short extracted genomic fragments are filtered out by the --min-seqeunce-length
        parameter.  

        All extracted genomic fragments with associated barcodes are exported to the fasta 
        file. Blat program with default parameters is used to map the extracted fragments 
        to the genomic DNA to identify the location of genomic framgents.

        To ensure a high quality mapping, several additional filters are applied. The minimal 
        number of gaps are required, which is controlled by --max-blat-gap-bases parameter.
        At least one block of the blat alignment should be long enough, which is controlled 
        by --min-blat-block-size. It is also required that the extracted genomic fragment is 
        mapped to one location on a genome. Thus, mappings to repeat regions are ignored.

        Bagseq program implements two additional types of filters to minimize the number of 
        erroneous barcodes and mappings to genomic DNA: i) sequence errors introduced by PCR, 
        ii) association of the same barcode to different genomic fragments

        To minimize the amount of erroneous barcodes cased by sequence errors, similar barcodes 
        were identified for each extracted and mapped barcode. Two barcodes are considered to 
        be similar if their sequnece is different by exactly one nucleotide. A given barcode is 
        considered to be true (recommended) barcode if all idnetified similar barcodes are less 
        frequent, and the ratio of frequencies of a given barcode and the most abundant similar 
        barcode is greater than a threshold defined by --sim-ratio-threshold parameter.
        
        Sometimes, several fastq sequence reads have the same barcode but different genomic
        fragments that can be caused by a number of experimental reasons (for example, randomly cloning 
        different fragments between barcode pairs). As a result, the same barcode
        will be mapped to different locations on a genome (with different frequecies).        
        To minimize the amount of erroneous barcode mappings, the number of
        reads supporting different locations for the same barcode were collected. 
        To distinguish the true location from the false one, the frequency of the most abundant
        location (the number of supported reads) was compared with frequencies of all other 
        locations for the same barcode. The most abundant location is considered to  
        be true (recommended) if the ratio of frequencies of this location
        and all other locations is greater than a treshold defined by --loc-ratio-threshold 
        parameter.        

        Bagseq program produces the following types of files:
        1. *.barcodes - file that has all extracted barcodes (for each source fastq file)
        2. up.fna, dn.fna  - fasta files with extracted genomic fragments associated barcodes
        3. up.psl, dn.psl  - output of blat program
        4. up.bhit.tsv, dn.bhit.tsv  - hits extracted from blat output (multiple hits 
            for the same genomic fragment are combined into one record)
        5. bagseq_up.tsv, bagseq_dn.tsv  - the main output of the bagsaeq program that has 
            barcodes, their most supported locations on a genomic DNA, statsitics  and recomendation 
            flags


        Examples to run the bagseq program:

        python -m dubseq.bagseq -u /path/to/up/fastq/files -d /path/to/dn/fastq/files 
               -g /path/to/genomic/dna/fasta/file  -o /output/dir                                    
        ''',
        formatter_class=util.RawDescriptionArgumentDefaultsHelpFormatter)

    parser.add_argument('-u', '--input-up',
                        dest='input_up',
                        help="""Upstream barcodes: fastq file (can be gzipped) or directory with fastq files (can be gzipped).
                        """,
                        type=str,
                        required=True
                        )

    parser.add_argument('-d', '--input-down',
                        dest='input_down',
                        help="""Downstream barcodes: fastq file (can be gzipped) or directory with fastq files (can be gzipped).
                        """,
                        type=str,
                        required=True
                        )

    parser.add_argument('-g', '--input-genome',
                        dest='input_genome',
                        help=""" Fasta file with nucleotide sequnce of a genome (can be more than one contig).
                        """,
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
                        dest='up_pos1',
                        help='Position of primer1 of the up tag',
                        default=15,
                        type=int
                        )

    parser.add_argument('--up-primer1-sequence',
                        dest='up_sequence1',
                        help='Sequence of primer1 of the up tag',
                        default='GAGGTCTCT',
                        type=str
                        )

    parser.add_argument('--up-primer2-pos',
                        dest='up_pos2',
                        help='Position of primer2 of the up tag',
                        default=44,
                        type=int
                        )

    parser.add_argument('--up-primer2-sequence',
                        dest='up_sequence2',
                        help='Sequence of primer2 of the up tag',
                        default='CGTACGCTG',
                        type=str
                        )

    parser.add_argument('--dn-primer1-pos',
                        dest='dn_pos1',
                        help='Position of primer1 of the down tag',
                        default=14,
                        type=int
                        )

    parser.add_argument('--dn-primer1-sequence',
                        dest='dn_sequence1',
                        help='Sequence of primer1 of the down tag',
                        default='GTCTCGTAG',
                        type=str
                        )

    parser.add_argument('--dn-primer2-pos',
                        dest='dn_pos2',
                        help='Position of primer2 of the down tag',
                        default=43,
                        type=int
                        )

    parser.add_argument('--dn-primer2-sequence',
                        dest='dn_sequence2',
                        help='Sequence of primer2 of the down tag',
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

    parser.add_argument('-l', '--min-seqeunce-length',
                        dest='min_sequence_length',
                        help=''' The minimum length of the extracted genomic sequence 
                        associated with a barcode.''',
                        type=int,
                        default=15
                        )

    parser.add_argument('-z', '--min-blat-block-size',
                        dest='min_blat_block_size',
                        help='''The size of at least one block in the blat alignemnt 
                        should be greater than this threshold. 
                        ''',
                        default=15,
                        type=int
                        )

    parser.add_argument('-a', '--max-blat-gap-bases',
                        dest='max_blat_gap_bases',
                        help=''' The total number of gaps in the blat alignment should not be 
                        greater than this threshold.
                        ''',
                        default=1,
                        type=int
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

    parser.add_argument('-r', '--loc-ratio-threshold',
                        dest='loc_ratio_threshold',
                        help=''' The minimum ratio of frequencies of the most abundant barcode location and 
                        all other barcode locations to consider the  most abundant barcode location to be true 
                        (not caused by chimeras created by PCR)
                        ''',
                        default=1,
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
    BarcodeLocation.CHIM_RATIO_THRESHOLD = Context.loc_ratio_threshold

    logging.info("Do up tags")
    process_barcodes(Context.up_tag_fastq_source, Context.up_tag,
                     Context.UP_TAG_FNAME_PREFIX)

    logging.info("Do dn tags")
    process_barcodes(Context.down_tag_fastq_source, Context.down_tag,
                     Context.DOWN_TAG_FNAME_PREFIX)

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


def process_barcodes(fastq_source, tag, fname_prefix):

    # process the source fastq files (extract barcodes, collect barcode stat, extract genomic sequences)
    process_fastq_files(fastq_source, tag, fname_prefix)

    # do blat on the extracted genomic sequences
    run_blat(fname_prefix)

    # process blat results and generates "bloc" file
    process_blat_results(fname_prefix)


def run_blat(fname_prefix):
    cmd = [
        Context.BLAT_CMD,
        Context.genome_fna_fname,
        Context.fasta_fname(fname_prefix),
        Context.blat_fname(fname_prefix)]
    subprocess.Popen(cmd).wait()


def process_blat_results(fname_prefix):
    '''
        Process the psl file with blat resuls of mapping barcode-associated
        genomic fragments to the genomic contigs
    '''
    # hashtable to accumulate genomic locations for each barcode
    # hastable { barcode sequence => BarcodeLocation}
    barcode_locations = {}

    bhit_fp = open(Context.bhit_fname(fname_prefix), 'w')
    try:
        bhit_fp.write("%s\n" % BarcodeHits.header())

        reader = BlatReader(Context.blat_fname(fname_prefix))
        try:
            prev_q_name = ''
            barcode_hits = BarcodeHits()
            blat_record = BlatRecord()
            while reader.next_record(blat_record):
                if blat_record.q_gapbases > Context.max_blat_gap_bases:
                    continue
                if blat_record.t_gapbases > Context.max_blat_gap_bases:
                    continue

                if blat_record.q_name != prev_q_name:
                    if barcode_hits.barcode:
                        process_barcode_hits(
                            bhit_fp, barcode_hits, barcode_locations)
                    barcode_hits(blat_record)
                    prev_q_name = blat_record.q_name
                else:
                    barcode_hits.add_hit(blat_record)

            if barcode_hits.barcode:
                process_barcode_hits(
                    bhit_fp, barcode_hits, barcode_locations)

        finally:
            reader.close()
    finally:
        bhit_fp.close()

    # Identify similar barcodes and collect number of reads supporintg similar barcodes
    BarcodeStat.find_similar_barcodes(barcode_locations)

    # Save barcode locations
    BarcodeLocation.save_barcodes_locations(
        Context.bagseq_output_fname(fname_prefix), barcode_locations)


def process_barcode_hits(bhit_fp, barcode_hits, barcode_locations):
    # store hits
    bhit_fp.write("%s\n" % str(barcode_hits))

    # collect location
    if barcode_hits.hits_count == 1 and max(barcode_hits.blat_max_block_sizes) >= Context.min_blat_block_size:
        barcode_location = barcode_locations.get(
            barcode_hits.barcode)
        if not barcode_location:
            barcode_location = BarcodeLocation()
            barcode_locations[barcode_hits.barcode] = barcode_location

        barcode_location.add_location(
            barcode_hits.blat_contig_ids[0],
            barcode_hits.blat_strands[0],
            barcode_hits.blat_starts[0],
            barcode_hits.blat_ends[0]
        )

        barcode_location.reads_count_inc()


class _SequenceIdGenerator:
    def __init__(self):
        self.__id = 0

    def next(self):
        self.__id += 1
        return self.__id


def process_fastq_files(fastq_source, tag, fname_prefix):
    '''
        As a result of the fastq files processeing (sepearately for the up and down tags)
        three types of info will be collected:

        1. Fasta file with extracted genomic regions (one combined file)
        2. File with barcode stat, including similar barcodes (one combined file)
        3. File with extracted barcodes (one per each fastq file)
    '''
    # # Hashtable to accumulate barcodes extracted from all fastq files
    # # it is a hashtable: {barcode sequnece => BarcodeStat}
    # barcode_stats = {}

    # to enerate a unique sequece id
    seq_id_generator = _SequenceIdGenerator()

    # Open fasta file to store extracted nucleotide sequence of genomic regions
    fasta_fp = open(Context.fasta_fname(fname_prefix), 'w')
    try:
        # proces all fastq files (extract genomic regions and colect barcodes)
        util.process_fastq_files(
            fastq_source,
            process_fastq_file,
            tag,
            seq_id_generator,
            fasta_fp)

    finally:
        fasta_fp.close()


def process_fastq_file(fastq_fname, tag, seq_id_generator, fasta_fp):
    '''
        It will:
        1. try to extract a barcode
        2. store the extracted barcode in barcodes file
        3. store the high quality barcode in barcode_stats for the downstream analysis
        4. extract a genomic sequnece for the high quality barcode and store in fasta file
    '''

    print("Doing file:%s" % fastq_fname)

    # To collect a general stats on the processed fastq file (will be logged)
    fastq_file_stat = FastqFileStat()

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
                barcode = tag.extract_barcode(
                    record, Context.primer_position_shifts)

                if barcode:
                    # count the reads with extracted barcodes
                    fastq_file_stat.barcode_extracted_reads_inc()

                    # store the extracted barcode
                    record_id = record.id.split(' ')[0]
                    barcodes_fp.write("%s\t%s\n" % (record_id, barcode))

                if barcode and barcode.min_quality >= Context.min_barcode_quality:

                    # Extract a genomic sequence
                    genomic_sequence = record.sequence[tag.tag_end:]

                    # if sequence is long enough, store in fasta file
                    if len(genomic_sequence) >= Context.min_genomic_sequence_length:
                        seq_id = seq_id_generator.next()
                        fasta_fp.write('>' + str(seq_id) + ":" +
                                       barcode.sequence + " " + record.id + '\n')
                        fasta_fp.write(genomic_sequence + '\n')

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
