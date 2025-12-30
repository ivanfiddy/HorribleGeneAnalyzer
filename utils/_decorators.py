def _checkifbuilt(method) -> callable:

    def checker(*args, **kwargs):

        d = {
            '_basedecode': '_base2bit',
            '_codondecode': '_base2bit',
            '_aminosdecode': '_aminos5bit',
            'baseratio': '_base2bit',
            'codonfreq': '_base2bit',
            'aminofreq': '_aminos5bit',
            '__getitem__': '_base2bit',
            '_orfsdecode': ('_aminos5bit', 'orfs_indexes'),
            'bandsgraph': 'orfs_indexes'
        }

        helpdict = {'_base2bit': "Haven't built _base2bit variable yet, use basebuild()",
                    "_aminos5bit": "Haven't built _trans5bit variable yet, use aminosbuild()",
                    ('_aminos5bit', 'orfs_indexes'): "Haven't built _aminos5bit and orfs_indexes, need to use aminosbuild() and orfsbuild()",
                    'orfs_indexes': "Haven't built orfs_indexes, need to use orfsbuild()"}

        attrib = d[method.__name__]

        if not isinstance(args[0], object):
            raise ValueError('decorator used outside of class scope (first self argument not given)')


        if isinstance(attrib, str) and type(args[0].__dict__.get(attrib)) == type(None):
            raise ValueError(helpdict[attrib])
        elif isinstance(attrib, tuple) and any([type(args[0].__dict__.get(i)) == type(None) for i in attrib]):
            raise ValueError(helpdict[attrib])
        else:
            return method(*args, **kwargs)

    return checker