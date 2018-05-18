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

        self.__browse_mode = 'gbrowse'
        self.__browse_mode_params = {
            'gbrowse': {
                'min_fscore': -5,
                'max_fscore': 20,
                'gene_y': 18,
                'plot_width': 15,
                'plot_height': 7,
                'plot_grid': True
            },
            'fbrowse': {
                'min_fscore': -5,
                'max_fscore': 20,
                'gene_y': 18,
                'plot_width': 15,
                'plot_height': 7,
                'plot_grid': True
            },
            'landscape': {
                'min_fscore': -1.5,
                'max_fscore': 3,
                'gene_y': 1.5,
                'plot_width': 15,
                'plot_height': 2,
                'plot_grid': False
            }
        }

        self.__color_model = 'gb'
        self.__color_model_params = {
            'gb': {
                'fr_covered_color': '#00FF00',
                'fr_non_covered_color': '#AAAAAA',
                'cur_gene_color': '#FF0000',
                'gene_color': '#000000',
                'gene_score_color': '#FF0000'
            }
        }

        self.__score_type = 'score_cnnls'
        self.__window_size = 14000
        self.__gene_x_offset = 200

        self.__itnum = None
        self.__fscores = None
        self.__gscores = None
        self.__cur_gene_index = 0

    def __getattr__(self, name):
        if name in self.__browse_mode_params[self.__browse_mode]:
            return self.__browse_mode_params[self.__browse_mode][name]
        elif name in self.__color_model_params[self.__color_model]:
            return self.__color_model_params[self.__color_model][name]

        raise AttributeError("DubSeq Viewer has no attribute '%s'" % name)

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

    @property
    def browse_mode(self):
        return self.__browse_mode

    def set_browse_mode(self, mode):
        if mode not in self.__browse_mode_params:
            raise AttributeError('%s - wrong browse mode. Available modes are: %s'
                                 % (mode, list(self.__browse_mode_params.keys())))
        self.__browse_mode = mode

    def set_color_model(self, cm):
        if cm not in self.__color_model_params:
            raise AttributeError('%s - wrong color model. Available models are: %s'
                                 % (cm, list(self.__color_model_params.keys())))
        self.__color_model = cm

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

    def show(self, fname=None):
        cur_gene = self.current_gene()

        (window_from, widnow_to) = self.window()

        fig = plt.figure(figsize=(self.plot_width, self.plot_height))
        ax = fig.add_subplot(111)

        # Do genes
        genes = self.genes(pos_from=window_from, pos_to=widnow_to)
        for _, gene in genes.iterrows():
            color = self.cur_gene_color  \
                if gene['gene_index'] == cur_gene['index'] \
                else self.gene_color

            arrowstyle = '->' if gene.strand == '+' else '<-'
            if self.browse_mode != 'landscape' or gene['gene_index'] == cur_gene['index']:
                ax.annotate(
                    gene['name'],
                    xy=(gene.pos_from + (gene.pos_to -
                                         gene.pos_from) / 2, self.gene_y),
                    fontsize=10,
                    xytext=(-10, 20), textcoords='offset points', ha='left', va='top',
                )
            ax.annotate(
                '',
                xy=(gene.pos_to, self.gene_y),
                xytext=(gene.pos_from, self.gene_y),
                fontsize=20,
                arrowprops=dict(arrowstyle=arrowstyle, color=color)
            )

        # Do fragments
        y_min = 0
        y_max = 0
        x_min = 0
        x_max = 0
        if self.browse_mode == 'landscape':
            fragments = self.fragments(pos_from=window_from, pos_to=widnow_to)
            y = self.min_fscore * 2 / 3
            y_range = np.abs(y * 2)
            y_min = y
            y_max = y + y_range
            x_min = fragments.pos_from.min()
            x_max = fragments.pos_to.max()

            for _, fragment in fragments.iterrows():

                color = self.fr_non_covered_color
                y += y_range / fragments.shape[0]
                ax.annotate(
                    '',
                    xy=(fragment.pos_to, y),
                    xytext=(fragment.pos_from, y),
                    fontsize=20,
                    arrowprops=dict(arrowstyle='-', color=color)
                )
        else:
            fscores = self.fscores(pos_from=window_from, pos_to=widnow_to)
            y_min = fscores.score.min()
            y_max = fscores.score.max()
            x_min = fscores.pos_from.min()
            x_max = fscores.pos_to.max()

            for _, fscore in fscores.iterrows():
                color = self.fr_covered_color \
                    if fscore.pos_from <= cur_gene.pos_from and fscore.pos_to >= cur_gene.pos_to \
                    else self.fr_non_covered_color
                ax.annotate(
                    '',
                    xy=(fscore.pos_to, fscore.score),
                    xytext=(fscore.pos_from, fscore.score),
                    fontsize=20,
                    arrowprops=dict(arrowstyle='-', color=color)
                )

        x_min = min(genes.pos_from.min(), x_min) - self.__gene_x_offset
        x_max = max(genes.pos_to.max(), x_max) + self.__gene_x_offset
        y_min = min(self.min_fscore, y_min)
        y_max = max(self.max_fscore, y_max)

        # Do gene score
        if self.browse_mode == 'gbrowse':
            gscore = cur_gene[self.__score_type]

            ax.annotate(
                '%s = %.2f' % (cur_gene['name'], gscore),
                xy=(x_max - self.__gene_x_offset, gscore),
                fontsize=10,
                xytext=(-70, 20), textcoords='offset points', ha='left', va='top', color=self.gene_score_color
            )
            ax.annotate(
                '',
                xy=(x_min + self.__gene_x_offset, gscore),
                xytext=(x_max - self.__gene_x_offset, gscore),
                fontsize=20,
                arrowprops=dict(arrowstyle='-', color=self.gene_score_color)
            )

        if self.browse_mode == 'landscape':
            plt.ylabel('')
            ax.set_yticklabels([])
        else:
            ax.set_title(self.current_condition()['name'], fontsize=15)
            plt.ylabel('fitness score')

        ax.grid(self.plot_grid)
        ax.axis([x_min, x_max, y_min, y_max])
        ax.get_xaxis().set_major_formatter(
            ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

        if fname is not None:
            plt.savefig(fname, type='pdf', bbox_inches='tight')
        else:
            plt.show()
