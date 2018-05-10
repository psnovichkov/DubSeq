import os
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class DubSeqViewer:

    def __init__(self, gscore_dir):
        self.__gscore_dir = gscore_dir
        self.__braseq_layout_df = pd.read_csv(
            os.path.join(gscore_dir, 'barseq_layout.tsv'), sep='\t')
        self.__fscore_base_df = pd.read_csv(
            os.path.join(gscore_dir, 'fscore_base.tsv'), sep='\t')
        self.__gscore_base_df = pd.read_csv(
            os.path.join(gscore_dir, 'gscore_base.tsv'), sep='\t')

        self.__window_size = 14000
        self.__min_fscore = -5
        self.__max_fscore = 20
        self.__gene_y = 18
        self.__gene_x_offset = 200

        self.__itnum = None
        self.__fscores = None
        self.__gscores = None
        self.__cur_gene_index = 0
        self.__score_type = 'score_cnnls'
        self.__fr_covered_color = '#00FF00'
        self.__fr_non_covered_color = '#AAAAAA'
        self.__cur_gene_color = '#FF0000'
        self.__gene_color = '#000000'
        self.__gene_score_color = '#FF0000'

    def set_score_type(self, score_type):
        self.__score_type = score_type

    def set_itnum(self, itnum):
        self.__itnum = itnum

        # Load fragment scores
        self.__fscores = pd.read_csv(os.path.join(
            self.__gscore_dir, '%s.fscore.tsv' % itnum), sep='\t')
        self.__fscores['pos_from'] = self.__fscore_base_df.pos_from
        self.__fscores['pos_to'] = self.__fscore_base_df.pos_to

        # Load gene scores
        self.__gscores = pd.read_csv(os.path.join(
            self.__gscore_dir, '%s.gscore.tsv' % itnum), sep='\t')
        self.__gscores['pos_from'] = self.__gscore_base_df.pos_from
        self.__gscores['pos_to'] = self.__gscore_base_df.pos_to
        self.__gscores['strand'] = self.__gscore_base_df.strand
        self.__gscores['name'] = self.__gscore_base_df['name']
        self.__gscores['locus_tag'] = self.__gscore_base_df.locus_tag
        self.__gscores['product'] = self.__gscore_base_df.product

        n = self.__gscore_base_df.shape[0]
        self.set_gene(index=self.__gscore_base_df.iloc[n // 2]['gene_index'])

    def set_gene(self, index=None, name=None, locus_tag=None):
        if index is not None:
            self.__cur_gene_index = index
        elif name is not None:
            genes = self.genes(name=name)
            if genes.shape[0] > 0:
                self.__cur_gene_index = genes.iloc[0]['gene_index']
        elif locus_tag is not None:
            genes = self.genes(locus_tag=locus_tag)
            if genes.shape[0] > 0:
                self.__cur_gene_index = genes.iloc[0]['gene_index']

    def __filter_range(self, d, pos_from, pos_to):
        if pos_from is not None:
            d = d[d.pos_from >= pos_from]
        if pos_to is not None:
            d = d[d.pos_to <= pos_to]
        return d

    def window(self):
        cur_gene = self.current_gene()
        gene_center = cur_gene.pos_from + \
            (cur_gene.pos_to - cur_gene.pos_from) / 2
        window_from = gene_center - self.__window_size / 2
        widnow_to = gene_center + self.__window_size / 2
        return (window_from, widnow_to)

    def set_window_size(self, window_size):
        self.__window_size = window_size

    def zoom_in(self):
        self.__window_size /= 1.2

    def zoom_out(self):
        self.__window_size *= 1.2

    def current_condition(self):
        d = self.__braseq_layout_df
        return d[d.itnum == self.__itnum].iloc[0]

    def current_gene(self):
        return self.__gscores.loc[self.__cur_gene_index]

    def next_gene(self):
        self.set_gene(index=self.__cur_gene_index + 1)

    def prev_gene(self):
        self.set_gene(index=self.__cur_gene_index - 1)

    def fscores(self, pos_from=None, pos_to=None):
        d = self.__fscores
        return self.__filter_range(d, pos_from, pos_to)

    def gscores(self, pos_from=None, pos_to=None):
        d = self.__gscores
        return self.__filter_range(d, pos_from, pos_to)

    def fragments(self, pos_from=None, pos_to=None):
        d = self.fscore_base
        return self.__filter_range(d, pos_from, pos_to)

    def genes(self, name=None, locus_tag=None, pos_from=None, pos_to=None):
        d = self.gscore_base
        if name is not None:
            d = d[d['name'].str.find(name) != -1]
        if locus_tag is not None:
            d = d[d['locus_tag'].str.find(locus_tag) != -1]
        return self.__filter_range(d, pos_from, pos_to)

    def conditions(self, name=None):
        d = self.braseq_layout
        if name is not None:
            d = d[d['name'].str.find(name) != -1]
        return d

    @property
    def fscore_base(self):
        return self.__fscore_base_df

    @property
    def gscore_base(self):
        return self.__gscore_base_df

    @property
    def braseq_layout(self):
        return self.__braseq_layout_df

    @property
    def gscore_dir(self):
        return self.__gscore_dir

    def show_next_gene(self):
        self.next_gene()
        self.show()

    def show_prev_gene(self):
        self.prev_gene()
        self.show()

    def show_gene(self, name=None, locus_tag=None):
        self.set_gene(name=name, locus_tag=locus_tag)
        self.show()

    def show_zoom_in(self):
        self.zoom_in()
        self.show()

    def show_zoom_out(self):
        self.zoom_out()
        self.show()

    def show(self):
        cur_gene = self.current_gene()

        (window_from, widnow_to) = self.window()
        genes = self.genes(pos_from=window_from, pos_to=widnow_to)
        fscores = self.fscores(pos_from=window_from, pos_to=widnow_to)

        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111)
#         ax.set_title('DubSeq Viewer: %s' % self.current_condition()['name'], fontsize=15)
        ax.set_title(self.current_condition()['name'], fontsize=15)
        ax.grid(True)

        # Do genes
        for _, gene in genes.iterrows():
            color = self.__cur_gene_color  \
                if gene['gene_index'] == cur_gene['index'] \
                else self.__gene_color

            arrowstyle = '->' if gene.strand == '+' else '<-'
            ax.annotate(
                gene['name'],
                xy=(gene.pos_from + (gene.pos_to - gene.pos_from) / 2, 18),
                fontsize=10,
                xytext=(-10, 20), textcoords='offset points', ha='left', va='top',
            )
            ax.annotate(
                '',
                xy=(gene.pos_to, 18),
                xytext=(gene.pos_from, 18),
                fontsize=20,
                arrowprops=dict(arrowstyle=arrowstyle, color=color)
            )

        # Do fragments
        for _, fscore in fscores.iterrows():
            color = self.__fr_covered_color \
                if fscore.pos_from <= cur_gene.pos_from and fscore.pos_to >= cur_gene.pos_to \
                else self.__fr_non_covered_color

            ax.annotate(
                '',
                xy=(fscore.pos_to, fscore.score),
                xytext=(fscore.pos_from, fscore.score),
                fontsize=20,
                arrowprops=dict(arrowstyle='-', color=color)
            )

        x_min = min(genes.pos_from.min(), fscores.pos_from.min()) - \
            self.__gene_x_offset
        x_max = max(genes.pos_to.max(),   fscores.pos_to.max()) + \
            self.__gene_x_offset
        y_min = min(self.__min_fscore, fscores.score.min())
        y_max = max(self.__max_fscore, fscores.score.max())

        # Do gene score
        gscore = cur_gene[self.__score_type]

        ax.annotate(
            '%s = %.2f' % (cur_gene['name'], gscore),
            xy=(x_max - self.__gene_x_offset, gscore),
            fontsize=10,
            xytext=(-70, 20), textcoords='offset points', ha='left', va='top', color=self.__gene_score_color
        )
        ax.annotate(
            '',
            xy=(x_min + self.__gene_x_offset, gscore),
            xytext=(x_max - self.__gene_x_offset, gscore),
            fontsize=20,
            arrowprops=dict(arrowstyle='-', color=self.__gene_score_color)
        )

        ax.axis([x_min, x_max, y_min, y_max])
        ax.get_xaxis().set_major_formatter(
            ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        plt.ylabel('fitness score')
        plt.show()
