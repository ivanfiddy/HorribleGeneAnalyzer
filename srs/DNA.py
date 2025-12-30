from io import TextIOWrapper
from bitarray import bitarray
import pandas as pnd

from utils import _slicer, _nspace, _checkifbuilt


class AnalyzeDNA:

    '''
        This class focuses on analysing the fundemental base sequence.
        It stores the base in a 2 bit array for optimized memory and parsing

        Builders methods
        Prerequisits for accessing public variables below:
        - basebuild() builds and stores fasta bases into a 2 bit array
        - transbuild() builds and stores translated bases into amino acids into a 5 bit array
        - reversebuild() builds and stores reverse template bases into a 2 bit array ( replaces storage of basebuild() )

        Public variables:
        - bases = index object that can be sliced to get bases / reverse complement bases quickly, e.g. bases[:] for all bases
        - codons = index object that can be sliced to get codons, e.g. codons[:200] for first 200 codons
        - trans = idnex object that can be sliecd to get amino acids, e.g. trans[::2] for every second amino acid

        Public methods:
        - baseratio(form=_) returns the ratio of bases. form is the fraction you want, that can be written in any way in
            the format of 'at/gc' or even by mixing bases into 'ag/ct'. Default is 'at/gc'
        - codonfreq() returns a pandas.Series() object with every codon and its frequency as found in the DNA sequence
        - transfrequi() returns a pandas.Series() object with every amino acid and its frequency

    '''


    ### Constructor
    def __init__(self, context=None) -> None:

        self.contextchecker(context)

        ## Public variables

        self.baselength = None
        self._aminos5bit = None

        ## Private variables

        # Not built private bitarray of bases
        self._base2bit = None
        # Not built private bitarray of amino acids (translated)
        self._aminos5bit = None
        # Callable indexer object
        self.indexer = _slicer
        # Store context
        self.context = context
        pass

    ### Public properties

    # Properties that can be sliced
    @property
    def codons(self) -> object:

        return self.indexer(self._codondecode)

    @property
    def bases(self) -> object:

        return self.indexer(self._basedecode)

    @property
    def aminos(self) -> object:

        return self.indexer(self._aminosdecode)

    ### Builders (public)

    def basebuild(self) -> None:

        self._bitgenebuilder()
        self.baselength = self._base2bit.__len__() / 2

        pass

    def aminosbuild(self) -> None:

        self._bitaminosbuilder()
        self.aminoslength = self._aminos5bit.__len__() / 5

        pass

    def reversebuild(self):

        self._reversebuilder()

        pass

    ### Public methods

    @_checkifbuilt
    def baseratio(self, form='at/gc') -> float:

        inp = form.replace(' ', '')

        frac = inp.split('/')

        if not len(inp) == 5 and not all([i in inp for i in ['a', 't', 'g', 'c', '/']]):
            raise ValueError('Incorrect formating for base ratio, use e.g. at/gc')

        numer = frac[0]
        denom = frac[1]

        g_c = self._base2bit.count(_nspace._base_encode[numer[0]]) \
              + self._base2bit.count(_nspace._base_encode[numer[0]])

        a_t = self._base2bit.count(_nspace._base_encode[denom[0]]) \
              + self._base2bit.count(_nspace._base_encode[denom[1]])

        return float(a_t / g_c)

    @_checkifbuilt
    def codonfreq(self) -> pnd.Series:

        codon_counter = {value: 0 for value in _nspace._bit6codon.keys()}

        it = self._base2bit.__iter__()

        for i in zip(it, it, it, it, it, it):
            codon_counter[i] += 1

        return pnd.Series({_nspace._bit6codon[key]: value for key, value in codon_counter.items()})

    @_checkifbuilt
    def aminofreq(self) -> pnd.Series:

        amino_counter = {keys: 0 for keys in _nspace._bit5amino.keys()}

        it = self._aminos5bit.__iter__()

        for i in zip(it, it, it, it, it):

            amino_counter[i] += 1

        return pnd.Series({_nspace._bit5amino[key]: value for key, value in amino_counter.items()})


    ### Script builders

    def _bitgenebuilder(self) -> None:

        temparray = bitarray()

        self.context.seek(0)

        for i in self.context.readlines()[1:]:
            temparray.encode(_nspace._base_encode, i.lower()[:-1])


        self._base2bit = temparray


    def _bitaminosbuilder(self) -> None:

        temparray = bitarray()

        self.context.seek(0)

        for i in self.context.readlines()[1:]:
            i = i.lower()[:-1]
            l = zip(i[0::3].__iter__(),
                    i[1::3].__iter__(),
                    i[2::3].__iter__())
            for j in l:
                temparray.extend(_nspace._basebit5amino[j])


        self._aminos5bit = temparray


    def _reversebuilder(self) -> None:

        temparray = bitarray()

        self.context.seek(0)

        for i in self.context.readlines()[:0:-1]:
            temparray.encode(_nspace._reverse_encode, i.lower()[:-1])

        self._base2bit = temparray


    ### Private methods

    @_checkifbuilt
    def _basedecode(self, index):

        if isinstance(index, int):

            return self._base2bit[
                        index*2:(index*2)+2
                        ].decode(_nspace._base_encode)

        elif isinstance(index, slice):

            newslice = slice(
                index.start*2 if index.start != None else None,
                index.stop*2 if index.stop != None else None,
                index.step*2 if index.step != None else None
            )

            return list(self._base2bit[newslice].decode(_nspace._base_encode))

        pass

    @_checkifbuilt
    def _codondecode(self, index):

        if isinstance(index, int):

            return ''.join(list(self._base2bit[
                        index*6:(index*6)+6
                        ].decode(_nspace._base_encode)))

        elif isinstance(index, slice):

            newslice = slice(
                index.start*6 if index.start != None else None,
                index.stop*6 if index.stop != None else None,
                index.step*6 if index.step != None else None
            )

            decoded_obj = list(self._base2bit[newslice].decode(_nspace._base_encode))

            wrapper = list(
                map(lambda x: ''.join(x),
                    zip(decoded_obj[0::3].__iter__(),
                        decoded_obj[1::3].__iter__(),
                        decoded_obj[2::3].__iter__())
                )
            )

            return wrapper

    @_checkifbuilt
    def _aminosdecode(self, index):

        if isinstance(index, int):

            return _nspace._bit5amino[tuple(
                self._aminos5bit[index * 5:(index * 5) + 5]
            )]

        elif isinstance(index, slice):

            newslice = slice(
                index.start * 5 if index.start != None else None,
                index.stop * 5 if index.stop != None else None,
                index.step * 5 if index.step != None else None
            )

            it = zip(self._aminos5bit[newslice].__iter__())

            return [_nspace._bit5amino[a + b + c + d + e] for a, b, c, d, e in zip(it, it, it, it, it)]

    def contextchecker(self, context):

        if context == None:
            raise ValueError("Analyze DNA takes in argument of file, open('__.fasta', 'r')")
        if not isinstance(context, TextIOWrapper):
            raise ValueError("The argument entered was not a type of file")


    ### Build in functions

    def __sizeof__(self) -> int:

        running_memory = 0

        for i in self.__dir__():

            try:
                running_memory += getattr(self, i).__sizeof__()

            except:
                pass


        return running_memory


    def __getitem__(self, index):

        return self.bases[index]