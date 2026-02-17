from srs import AnalyzeDNA
from bitarray import bitarray
import pandas as pnd
from matplotlib import pyplot as plt

from utils import _slicer, _nspace, _checkifbuilt

class AnalyzeORFs:

    ### Constructor

    def __init__(self, context):

        ## Private variables


        # Parent inheritance of AnalyzeDNA for amino and base bitarrays
        self._aminos5bit = None
        self.context = context
        self.indexer = _slicer
        self.orfs_indexes = None
        pass


    ### Public properties

    @property
    def orfs(self):

        return self.indexer(self._orfsdecode)

        pass

    ### Public builders

    def aminosbuild(self) -> None:

        inherit = AnalyzeDNA(self.context)
        inherit.aminosbuild()
        self._aminos5bit = inherit._aminos5bit
        del inherit


    def orfsbuild(self) -> None:

        start_codon = tuple(_nspace._aminobit5['MET'])
        stop_codon = tuple(_nspace._aminobit5['STP'])

        zipper = (
            self._aminos5bit[i::5].__iter__() for i in range(5)
        )

        overlap_indexes = []

        orfs_indexes = []

        running_orf = False

        index_overlap = 1

        running_pos = 0


        for i in zip(*zipper):

            if running_orf:

                if i == start_codon:
                    index_overlap += 1
                    orfs_indexes.append([running_pos, None])
                    overlap_indexes.append(
                        len(orfs_indexes) - 1
                    )

                elif i == stop_codon:

                    if index_overlap < 2:
                        orfs_indexes[-1] = orfs_indexes[-1] = (orfs_indexes[-1][0], running_pos + 1)

                    else:
                        for j in range(index_overlap):
                            orfs_indexes[-(j + 1)] = (orfs_indexes[-(j + 1)][0], running_pos + 1)

                        index_overlap = 1

                    running_orf = False

            elif i == start_codon:

                running_orf = True
                orfs_indexes.append([running_pos, None])


            running_pos += 1


        self.orfs_indexes = pnd.DataFrame(tuple(orfs_indexes), columns=['start', 'stop'])
        self.orfs_indexes['length'] = self.orfs_indexes['stop'] - self.orfs_indexes['start']
        self.orfs_indexes['overlap'] = False
        self.orfs_indexes.loc[overlap_indexes, 'overlap'] = True


    ### Band graph

    @_checkifbuilt
    def bandsgraph(self, start=0, stop=None, chart_height=10, chart_width=10, overlap=False):

        indexes = self.orfs_indexes[
            self.orfs_indexes['overlap'] == False
        ]

        per_row = round(indexes.__len__() / chart_height)

        slices = (slice(i*per_row, ((i+1)*per_row)+1)
                  for i in range(chart_height))


        fig, ax = plt.subplots()


        y_pos_start = -0.2
        some_argument_idk = 0.8
        prev = 0

        for y_pos in range(chart_height):

            line_bars = indexes.iloc[next(slices)]

            bar_data = line_bars.loc[:, ['start', 'length']]

            last_data = bar_data.iloc[-1]

            row_data_length = last_data['start'] + last_data['length']

            # I feel like this is not written optimally or most concisley
            bar_data['start'] = bar_data['start'].map(lambda x: float((x - prev) / (row_data_length - prev) * chart_width))
            bar_data['length'] = bar_data['length'].map(lambda x: float(x / (row_data_length - prev) * chart_width))

            bar_plot = bar_data.itertuples(index=False, name=None)

            ax.broken_barh(tuple([*bar_plot]), (y_pos_start, some_argument_idk), picker=True).set_gid(line_bars.to_dict(orient='records'))

            y_pos_start += 1
            prev = row_data_length

        ax.set_yticks(range(chart_height))
        ax.invert_yaxis()
        ax.set_xlim(0, chart_width)

        plt.show()










    ### Private methods
    @_checkifbuilt
    def _orfsdecode(self, index):

        if isinstance(index, int):

            pos = self.orfs_indexes.iloc[index]

            it = self._aminos5bit[pos['start']*5:pos['stop']*5].__iter__()

            return tuple(_nspace._bit5amino[i] for i in zip(it, it, it, it, it))

        elif isinstance(index, slice):

            newslice = slice(
                index.start if index.start != None else None,
                index.stop if index.stop != None else None,
                index.step if index.step != None else None
            )

            poses = self.orfs_indexes.iloc[newslice]

            orfs = []

            for index, row in poses.iterrows():

                it = self._aminos5bit[row['start']*5:row['stop']*5].__iter__()
                orfs.append(tuple(_nspace._bit5amino[i] for i in zip(it, it, it, it, it)))

            return tuple(orfs)

        pass
