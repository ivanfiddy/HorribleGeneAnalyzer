class _slicer:

    def __init__(self, wrapper):
        self.wrapper = wrapper


    def __getitem__(self, index):

        return self.wrapper(index)

