import struct
from . import util
from .blat import BlatRecord


class Primer:
    __slots__ = ('__sequence', '__pos')

    def __init__(self, sequence, pos):
        self.__sequence = sequence
        self.__pos = pos

    @property
    def sequence(self):
        return self.__sequence

    @property
    def pos(self):
        return self.__pos

    @property
    def size(self):
        return len(self.__sequence)

    def check_primer(self, sequence, pos_shift, require_entire_primer=True):
        pos_from = self.__pos + pos_shift
        pos_to = pos_from + self.size

        seq_len = len(sequence)

        if pos_from >= seq_len:
            return False

        if pos_to > seq_len and not require_entire_primer:
            primer_to = len(self.sequence) - (pos_to - seq_len)
            pos_to = seq_len
            return self.sequence[0:primer_to] == sequence[pos_from: pos_to]

        return self.sequence == sequence[pos_from: pos_to]


class BarcodeTag:
    def __init__(self, primer_seq1, primer_pos1, primer_seq2, primer_pos2):
        self.__primer1 = Primer(primer_seq1, primer_pos1)
        self.__primer2 = Primer(primer_seq2, primer_pos2)

    @property
    def primer1(self):
        return self.__primer1

    @property
    def primer2(self):
        return self.__primer2

    def check_primers(self, sequence, pos_shift, require_entire_primer2=True):
        return self.primer1.check_primer(sequence, pos_shift) and \
            self.primer2.check_primer(
                sequence, pos_shift, require_entire_primer2)

    @property
    def tag_start(self):
        return self.primer1.pos

    @property
    def tag_end(self):
        return self.primer2.pos + self.primer2.size

    def extract_barcode(self, fastq_record, pos_shifts, require_entire_primer2=True):
        barcode = None
        for pos_shift in pos_shifts:
            if self.check_primers(fastq_record.sequence, pos_shift, require_entire_primer2):
                barcode = Barcode()
                barcode(fastq_record, self, pos_shift)
                break
        return barcode


class Barcode:
    __BARCODE_SIZE = 20
    __QUALITY_CHAR_BASE = 33
    __MAX_QUALITY = 100

    def __init__(self):
        self.__pos = 0
        self.__sequence = '-' * Barcode.__BARCODE_SIZE
        self.__quality_str = '-' * Barcode.__BARCODE_SIZE
        self.__min_quality = 0
        # self.__hex_code = '-' * int(Barcode.__BARCODE_SIZE / 2)

    def __call__(self, fastq_record, barcode_tag, pos_shift):
        pos_from = barcode_tag.primer1.pos + barcode_tag.primer1.size + pos_shift
        pos_to = barcode_tag.primer2.pos + pos_shift
        self.__pos = pos_from
        self.__sequence = fastq_record.sequence[pos_from: pos_to]
        self.__quality_str = fastq_record.quality[pos_from: pos_to]
        # self.__hex_code = util.to_hex_code(self.sequence)
        self.__min_quality = self._min_quality(self.quality_str)

    @staticmethod
    def header(prefix='', sep='\t'):
        return sep.join(prefix + x for x in
                        [
                            # 'hex_code',
                            'sequence', 'pos', 'quality_str', 'min_quality'])

    def __str__(self):
        return "\t".join([
            # self.hex_code,
            self.sequence, str(self.pos), self.quality_str, str(self.min_quality)])

    @property
    def pos(self):
        return self.__pos

    @property
    def sequence(self):
        return self.__sequence

    @property
    def quality_str(self):
        return self.__quality_str

    @property
    def min_quality(self):
        return self.__min_quality

    # @property
    # def hex_code(self):
    #     return self.__hex_code

    def _min_quality(self, quality_str):
        minQuality = Barcode.__MAX_QUALITY
        for ch in quality_str:
            quality = ord(ch) - Barcode.__QUALITY_CHAR_BASE
            if quality < minQuality:
                minQuality = quality
        return minQuality


EMPTY_BARCODE = Barcode()


class BarcodeStat:
    SIM_RATIO_THRESHOLD = 2

    __slots__ = ['__reads_count',
                 '__sim_reads_counts'
                 ]

    @staticmethod
    def header(prefix='', sep='\t'):
        return sep.join(prefix + x for x in [
            'reads_count',
            'sim_recommended',
            'sim_reads_count',
            'sim_reads_count_max',
            'sim_reads_counts'
        ])

    def __init__(self):
        self.__reads_count = 0
        self.__sim_reads_counts = []

    def __str__(self):
        return '\t'.join(str(x) for x in [
            self.reads_count,
            '+' if self.sim_recommended() else '-',
            self.sim_reads_count,
            self.sim_reads_count_max,
            ','.join(str(x) for x in self.__sim_reads_counts) + ','
        ])

    @property
    def reads_count(self):
        return self.__reads_count

    @property
    def sim_reads_count(self):
        return sum(self.__sim_reads_counts)

    @property
    def sim_reads_count_max(self):
        return max(self.__sim_reads_counts) if self.__sim_reads_counts else 0

    def reads_count_inc(self):
        self.__reads_count += 1

    def add_sim_reads_count(self, sim_reads_count):
        self.__sim_reads_counts.append(sim_reads_count)

    def sim_recommended(self):
        return self.__reads_count >= self.sim_reads_count_max * BarcodeStat.SIM_RATIO_THRESHOLD

    @staticmethod
    def find_similar_barcodes(barcodes):
        # update info about similar barcodes
        for barcode_sequence, barcode_stat in barcodes.items():
            barcode_sequence_chars = list(barcode_sequence)
            util.process_similar_sequences(
                barcode_sequence_chars,
                BarcodeStat._process_similar_barcode,
                barcode_stat,
                barcodes)

    @staticmethod
    def _process_similar_barcode(chars, barcode_stat, barcodes):
        similar_barcode_key = ''.join(chars)
        similar_barcode_stat = barcodes.get(similar_barcode_key)
        if similar_barcode_stat:
            barcode_stat.add_sim_reads_count(similar_barcode_stat.reads_count)
            # barcode_stat.sim_reads_count += similar_barcode_stat.reads_count

    @staticmethod
    def save_barcode_stats(file_name, barcodes):
        with open(file_name, 'w') as f:
            # write a header
            f.write("%s\t%s\n" % ('barcode', BarcodeStat.header()))
            for barcode_sequence, barcode_stat in barcodes.items():
                f.write("%s\t%s\n" % (barcode_sequence, str(barcode_stat)))


class PairedBarcodeStat(BarcodeStat):
    CHIM_RATIO_THRESHOLD = 2

    __slots__ = ['__pair_reads_counts']

    @staticmethod
    def header(prefix='', sep='\t'):
        return BarcodeStat.header(prefix, sep) + sep + sep.join(prefix + x for x in [
            'chim_recommended',
            'pair_reads_count_max',
            'pair_reads_count_submax',
            'pair_reads_counts'
        ])

    def __init__(self):
        BarcodeStat.__init__(self)
        self.__pair_reads_counts = []

    def __str__(self):
        return BarcodeStat.__str__(self) + '\t' + '\t'.join(str(x) for x in [
            '+' if self.chim_recommended() else '-',
            self.pair_reads_count_max,
            self.pair_reads_count_submax,
            ','.join(str(x) for x in self.__pair_reads_counts) + ','
        ])

    def chim_recommended(self):
        return self.pair_reads_count_max >= self.pair_reads_count_submax * PairedBarcodeStat.CHIM_RATIO_THRESHOLD

    def add_pair_reads_count(self, pair_reads_count):
        self.__pair_reads_counts.append(pair_reads_count)

    @property
    def pair_reads_count_max(self):
        # if self.__pair_reads_counts else 0
        return max(self.__pair_reads_counts)

    @property
    def pair_reads_count_submax(self):
        m1 = 0
        m2 = 0
        for count in self.__pair_reads_counts:
            if count > m1:
                m2 = m1
                m1 = count
            elif count > m2:
                m2 = count
        return m2

    @staticmethod
    def save_barcode_stats(file_name, barcodes):
        with open(file_name, 'w') as f:
            # write a header
            f.write("%s\t%s\n" % ('barcode', PairedBarcodeStat.header()))
            for barcode_sequence, barcode_stat in barcodes.items():
                f.write("%s\t%s\n" % (barcode_sequence, str(barcode_stat)))

        # class BarcodeStat:
        #     """ To store barcode statistics
        #         Optimizing memory usage by using __slots__ and packing into a single variable
        #     """

        #     __slots__ = ('__data')
        #     __counts = struct.Struct('ll')

        #     @staticmethod
        #     def header(prefix='', sep='\t'):
        #         return sep.join(prefix + x for x in [
        #             'reads_count',
        #             'sim_reads_count',
        #             'recommended'
        #         ])

        #     def __init__(self):
        #         self.__data = BarcodeStat.__counts.pack(0, 0)

        #     def __str__(self):
        #         return '\t'.join(str(x) for x in [
        #             self.reads_count,
        #             self.sim_reads_count,
        #             '+' if self.recommended() else '-'
        #         ])

        #     @property
        #     def reads_count(self):
        #         return self.__getter(0)

        #     @reads_count.setter
        #     def reads_count(self, value):
        #         self.__setter(0, value)

        #     @property
        #     def sim_reads_count(self):
        #         return self.__getter(1)

        #     @sim_reads_count.setter
        #     def sim_reads_count(self, value):
        #         self.__setter(1, value)

        #     def __getter(self, index):
        #         return BarcodeStat.__counts.unpack(self.__data)[index]

        #     def __setter(self, index, value):
        #         values = list(BarcodeStat.__counts.unpack(self.__data))
        #         values[index] = value
        #         self.__data = BarcodeStat.__counts.pack(*values)

        #     def recommended(self):
        #         return self.reads_count > self.sim_reads_count

        #     @staticmethod
        #     def find_similar_barcodes(barcodes):
        #         # update info about similar barcodes
        #         for barcode_sequence, barcode_stat in barcodes.items():
        #             barcode_sequence_chars = list(barcode_sequence)
        #             util.process_similar_sequences(
        #                 barcode_sequence_chars,
        #                 BarcodeStat._process_similar_barcode,
        #                 barcode_stat,
        #                 barcodes)

        #     @staticmethod
        #     def _process_similar_barcode(chars, barcode_stat, barcodes):
        #         similar_barcode_key = ''.join(chars)
        #         similar_barcode_stat = barcodes.get(similar_barcode_key)
        #         if similar_barcode_stat:
        #             barcode_stat.sim_reads_count += similar_barcode_stat.reads_count

        #     @staticmethod
        #     def save_barcode_stats(file_name, barcodes):
        #         with open(file_name, 'w') as f:
        #             # write a header
        #             f.write("%s\t%s\n" % ('barcode', BarcodeStat.header()))
        #             for barcode_sequence, barcode_stat in barcodes.items():
        #                 f.write("%s\t%s\n" % (barcode_sequence, str(barcode_stat)))


class BarcodeLocation(BarcodeStat):
    LOC_RATIO_THRESHOLD = 0.3

    __slots__ = ['__location_2_count']

    def __init__(self):
        BarcodeStat.__init__(self)
        self.__location_2_count = {}

    @staticmethod
    def header(sep='\t'):
        return BarcodeStat.header(sep) + sep + sep.join([
            'loc_recommended',
            'total_hits_count',
            'hit_count_max',
            'hit_count_submax',
            'top_hit_sequence_id',
            'top_hit_pos',
            'top_hit_strand',
            'hits_counts'
        ])

    def __str__(self):
        hit_count_max = 0
        hit_count_submax = 0
        top_hit_sequence_id = ''
        top_hit_pos = 0
        top_hit_strand = '-'
        hits_counts = []

        for location_id, hit_count in self.__location_2_count.items():

            hits_counts.append(hit_count)
            if hit_count > hit_count_max:
                hit_count_submax = hit_count_max

                hit_count_max = hit_count
                top_hit_sequence_id, top_hit_pos = location_id[1:].split(":")
                top_hit_strand = location_id[0]
            elif hit_count > hit_count_submax:
                hit_count_submax = hit_count

        loc_recommended = hit_count_max >= hit_count_submax * \
            BarcodeLocation.LOC_RATIO_THRESHOLD

        return BarcodeStat.__str__(self) + '\t' + '\t'.join(str(x) for x in
                                                            [
            '+' if loc_recommended else '-',
            self.total_hits_count,
            hit_count_max,
            hit_count_submax,
            top_hit_sequence_id,
            top_hit_pos,
            top_hit_strand,
            ','.join(str(x) for x in hits_counts)
        ])

    @property
    def total_hits_count(self):
        count = 0
        for hits_count in self.__location_2_count.values():
            count += hits_count
        return count

    def add_location(self, sequence_id, strand, start, end):
        location_id = "%s%s:%s" % (
            strand,
            sequence_id,
            str(start) if strand == '+' else str(end)
        )

        location_count = self.__location_2_count.get(location_id)
        self.__location_2_count[location_id] = location_count + \
            1 if location_count else 1

    @staticmethod
    def save_barcodes_locations(file_name, barcodes_locations):
        with open(file_name, 'w') as f:

            # write header
            f.write("%s\t%s\n" % ('barcode', BarcodeLocation.header()))

            # write all barcodes
            for barcode, barcode_location in barcodes_locations.items():
                f.write("%s\t%s\n" % (barcode, str(barcode_location)))


class BarcodeHits:

    def __init__(self):
        self.barcode = None
        self.barcode_index = None
        self.hits_count = None
        self.blat_contig_ids = None
        self.blat_matches = None
        self.blat_strands = None
        self.blat_starts = None
        self.blat_ends = None
        self.blat_block_counts = None
        self.blat_max_block_sizes = None

    def __call__(self, blat_record):
        self.barcode_index, self.barcode = blat_record.q_name.split(':')
        self.hits_count = 0
        self.blat_contig_ids = []
        self.blat_matches = []
        self.blat_strands = []
        self.blat_starts = []
        self.blat_ends = []
        self.blat_block_counts = []
        self.blat_max_block_sizes = []
        self.add_hit(blat_record)

    @staticmethod
    def header(sep='\t'):
        return sep.join([
            'barcode',
            'barcode_index',
            'hits_count',
            'blat_contig_ids',
            'blat_matches',
            'blat_strands',
            'blat_starts',
            'blat_ends',
            'blat_block_counts',
            'blat_max_block_sizes'
        ])

    def __str__(self):
        return '\t'.join(
            [
                self.barcode,
                str(self.barcode_index),
                str(self.hits_count),
                ','.join(str(x) for x in self.blat_contig_ids),
                ','.join(str(x) for x in self.blat_matches),
                ','.join(str(x) for x in self.blat_strands),
                ','.join(str(x) for x in self.blat_starts),
                ','.join(str(x) for x in self.blat_ends),
                ','.join(str(x) for x in self.blat_block_counts),
                ','.join(str(x) for x in self.blat_max_block_sizes)
            ]
        )

    def add_hit(self, blat_record):
        self.hits_count += 1
        self.blat_contig_ids.append(blat_record.t_name)
        self.blat_matches.append(blat_record.match)
        self.blat_strands.append(blat_record.strand)
        self.blat_starts.append(blat_record.t_start)
        self.blat_ends.append(blat_record.t_end)
        self.blat_block_counts.append(blat_record.block_count)
        self.blat_max_block_sizes.append(max(blat_record.block_sizes))
