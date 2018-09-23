import os.path
import pickle


class Parabolics:
    type = 0

    # underlying weyl group
    weyl = 0

    # standard parabolics
    # - simple roots
    simple_roots = 0
    # - all roots
    roots = 0
    # - longest element
    w = 0
    # - w0 * longest element of para
    w0w = 0

    def __init__(self, t, with_roots=True, with_weyl=True):
        self.type = t

        self.load_parabolics()

        if with_roots:
            self.roots = t.rootsystem
        if with_weyl:
            self.weyl = t.weyl

    """
    Load data for standard parabolics from ... para.pickle
    """

    def load_parabolics(self):
        fname = "data/" + self.type.label + "/para.pickle"
        print("Loading:", fname)

        if os.path.isfile(fname):
            f = open(fname, 'rb')
            tmp = pickle.load(f)
            f.close()
            self.simple_roots = tmp[0]
            self.roots = tmp[1]
            self.w = tmp[2]
            self.w0w = tmp[3]
        else:
            print("File not found: " + fname)
