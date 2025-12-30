import builtins

def getsize(obj: object):

    if not '__sizeof__' in obj.__dir__():
        raise ValueError("Object has no attribute '__getsize__'")

    byte_size = obj.__sizeof__()


    index = 0
    units = ['bytes', 'kilobytes', 'megabytes', 'gigabytes']

    outp_size = float(byte_size)

    while outp_size // 1_000 > 0:
        outp_size /= 1024
        index += 1

    if index > len(units) - 1:
        raise ValueError('Object exceeds the maximum size measurements of gigabytes')

    return (outp_size, units[index])