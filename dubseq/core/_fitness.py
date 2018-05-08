import os
import pandas as pd


class BarseqLayoutItem:
    def __init__(self, itnum, item_type, experiment_condition):
        self.__itnum = itnum
        self.__item_type = item_type
        self.__experiment_condition = experiment_condition

    @property
    def itnum(self):
        return self.__itnum

    @property
    def item_type(self):
        return self.__item_type

    @property
    def experiment_condition(self):
        return self.__experiment_condition


class BarseqLayout:
    def __init__(self, layout_file_name):
        self.__layout_file_name = layout_file_name
        self.__df = None
        self.__load()

        # Check if the layout has time zero items
        if len(self.time_zero_items) == 0:
            raise ValueError(
                'No time zero experiments were found in the layout')

    def __load(self):
        self.__df = pd.read_csv(self.__layout_file_name, sep='\t')

    @property
    def layout_file_name(self):
        return self.__layout_file_name

    @property
    def time_zero_items(self):
        df = self.__df
        return self.__to_items(df[df.type == 'Time0'])

    @property
    def lb_items(self):
        df = self.__df
        return self.__to_items(df[df.type == 'LB'])

    @property
    def stress_items(self):
        df = self.__df
        return self.__to_items(df[df.type == 'stress'])

    @property
    def non_time_zero_items(self):
        df = self.__df
        return self.__to_items(df[df.type != 'Time0'])

    @property
    def all_items(self):
        return self.__to_items(self.__df)

    def __to_items(self, df):
        items = []
        for _, row in df.iterrows():
            items.append(BarseqLayoutItem(row.itnum, row.type, row.name))
        return items

    @property
    def experiment_types(self):
        return self.__df.type.unique()


class BpagItem:
    __slots__ = ['barcode_up', 'barcode_dn',
                 'bpair_read_count', 'up_read_count', 'dn_read_count',
                 'contig_id',  'pos_from', 'pos_to']

    def __init__(self, barcode_up, barcode_dn, bpair_read_count,
                 up_read_count, dn_read_count,
                 contig_id, pos_from, pos_to):
        self.barcode_up = barcode_up
        self.barcode_dn = barcode_dn
        self.bpair_read_count = bpair_read_count
        self.up_read_count = up_read_count
        self.dn_read_count = dn_read_count
        self.contig_id = contig_id
        self.pos_from = pos_from
        self.pos_to = pos_to


class BpagSet:
    def __init__(self, blag_file_name):
        self.__blag_file_name = blag_file_name
        self.__items = []
        self.__up_barcode_2_item = {}
        self.__load()

    @property
    def blag_file_name(self):
        return self.__blag_file_name

    @property
    def size(self):
        return len(self.__items)

    def get_item(self, index):
        return self.__items[index]

    def find_up_item(self, barcode_up):
        return self.__up_barcode_2_item.get(barcode_up)

    def __load(self):
        df = pd.read_csv(self.__blag_file_name, sep='\t')
        for _, row in df.iterrows():
            if row.recommended == '+':
                item = BpagItem(
                    row.barcode_up,
                    row.barcode_dn,
                    row.bpair_read_count,
                    row.up_read_count,
                    row.dn_read_count,
                    row.up_contig_id,
                    row.pos_from,
                    row.pos_to
                )
                self.__items.append(item)
                self.__up_barcode_2_item[item.barcode_up] = item


class TimeZeroItem:
    def __init__(self, barcode, time0_experiments_count):
        self.__barcode = barcode
        self.__read_counts = [0] * time0_experiments_count

    @property
    def barcode(self):
        return self.__barcode

    @property
    def total_read_count(self):
        return sum(self.__read_counts)

    @property
    def max_read_count(self):
        return max(self.__read_counts)

    def set_read_count(self, experiment_index, count):
        self.__read_counts[experiment_index] = count


class TimeZeroSet:
    def __init__(self, bpag_set, barseq_layout, barseq_dir):
        self.__time0_itnums = []
        for item in barseq_layout.time_zero_items:
            self.__time0_itnums.append(item.itnum)

        self.__barcode_2_item = {}
        self.__items = []
        self.__load(bpag_set, barseq_dir)

    def filter_items(self, good_item_method):
        for i in range(self.size)[::-1]:
            if not good_item_method(self.__items[i]):
                barcode = self.__items[i].barcode
                del self.__barcode_2_item[barcode]
                del self.__items[i]

    @property
    def size(self):
        return len(self.__items)

    def __load(self, bpag_set, barseq_dir):
        for experiment_index, itnum in enumerate(self.__time0_itnums):
            bstat_fname = self.__get_bstat_file(itnum, barseq_dir)
            if not bstat_fname:
                raise ValueError(
                    'Can not find bstat file for itnum %s in %s directory' % (itnum, barseq_dir))

            df = pd.read_csv(bstat_fname, sep='\t')
            for _, row in df.iterrows():

                if row.recommnended != '+':
                    continue

                if not bpag_set.find_up_item(row.barcode):
                    continue

                self.__register_read_count(
                    row.barcode, experiment_index, int(row.reads_count))

    def __get_bstat_file(self, itnum, barseq_dir):
        file_path = None
        for file_name in os.listdir(barseq_dir):
            if not file_name.endswith('.bstat.tsv'):
                continue
            if '_' + itnum + '_' in file_name:
                file_path = os.path.join(barseq_dir, file_name)
                break

        return file_path

    @property
    def experiment_count(self):
        return len(self.__time0_itnums)

    def __register_read_count(self, barcode, exp_index, read_count):
        t0_item = self.__barcode_2_item.get(barcode)
        if not t0_item:
            t0_item = TimeZeroItem(barcode, self.experiment_count)
            self.__barcode_2_item[barcode] = t0_item

        if t0_item:
            t0_item.set_read_count(exp_index, read_count)
