from bitarray import bitarray
from bitarray.util import int2ba
from types import SimpleNamespace

def _build_namespace():
    _nspace = SimpleNamespace()

    _nspace._base_encode = {
        'a': bitarray('00'),
        't': bitarray('01'),
        'g': bitarray('10'),
        'c': bitarray('11')
    }

    _nspace._reverse_encode = {
        'a': bitarray('01'),
        't': bitarray('00'),
        'g': bitarray('11'),
        'c': bitarray('10')
    }

    _nspace._codonbase = {
        'a': 't',
        'u': 'a',
        'c': 'g',
        'g': 'c'
    }

    _nspace._aminocodon_dict = {
        'phe': ['uuu', 'uuc'],
        'leu': ['uua', 'uug', 'cuu', 'cuc', 'cua', 'cug'],
        'ile': ['auu', 'auc', 'aua'],
        'MET': ['aug'],
        'val': ['guu', 'guc', 'gua', 'gug'],
        'ser': ['ucu', 'ucc', 'uca', 'ucg', 'agu', 'agc'],
        'pro': ['ccu', 'ccc', 'cca', 'ccg'],
        'thr': ['acu', 'acc', 'aca', 'acg'],
        'ala': ['gcu', 'gcc', 'gca', 'gcg'],
        'tur': ['uau', 'uac'],
        'STP': ['uaa', 'uag', 'uga'],
        'his': ['cau', 'cac'],
        'gln': ['caa', 'cag'],
        'asn': ['aau', 'aac'],
        'lys': ['aaa', 'aag'],
        'asp': ['gau', 'gac'],
        'glu': ['gaa', 'gag'],
        'cys': ['ugu', 'ugc'],
        'trp': ['ugg'],
        'arg': ['cgu', 'cgc', 'cga', 'cgg', 'aga', 'agg'],
        'gly': ['ggu', 'ggc', 'gga', 'ggg']
    }

    _nspace._codonamino_dict = {j: i[0] for i in _nspace._aminocodon_dict.items() for j in i[1]}

    amino_iterate = _nspace._aminocodon_dict.keys().__iter__()

    _nspace._aminobit5 = {next(amino_iterate): int2ba(i, 5) for i in range(21)}

    _nspace._bit5amino = {tuple(value): key for key, value in _nspace._aminobit5.items()}

    codon_iterate = _nspace._codonamino_dict.keys().__iter__()
    _nspace._codonbit6 = {next(codon_iterate): int2ba(i, 6) for i in range(64)}

    _nspace._bit6codon = {tuple(value): key for key, value in _nspace._codonbit6.items()}

    _nspace._basebit5amino = {tuple([_nspace._codonbase[i] for i in key]): _nspace._aminobit5[value] for key, value in _nspace._codonamino_dict.items()}

    _nspace._bit6codonbit5amino = {tuple(_nspace._codonbit6[key]): _nspace._aminobit5[value]
                                   for key, value in _nspace._codonamino_dict.items()}

    return _nspace