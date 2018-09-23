import os.path
import h5py
import pickle
import numpy as np


class Weyl:
    type = 0

    # underlying roots
    roots = 0

    # permutations of the Weyl group in the form
    # the image of the permuted roots
    perm = 0
    # the number of permutations, i.e. order of W
    nr_perm = 0
    # the order of the permutations
    perm_ord = 0
    # the inverse of the permutations
    perm_inv = 0

    # permutations of the Weyl group in the form
    # of words of simple permutations
    # --- here they are stored as strings,
    # --- because of hdf5 doesn't allow arbitrary length
    word_str = 0
    # the length of a string representation
    word_str_len = 0
    # the length of the individual words
    word_len = 0
    # the indices of the reflections
    # --- the first rang refs are simple refs
    # --- refs[r] is the index of the reflection in the root r
    # --- although it is enough to know them for r positive, we hold also -r
    refs = 0
    # the longest element:
    w0 = 0

    # tree for distinguished expressions
    expr_tree = 0
    # list of distinguished expressions
    # --- [[x,y,z],[dist expr]] such that
    # ------ expr_in_tree(x,y,z) is not []
    # ------ and non-distinguished expressions are filtered out
    dist_expr = 0
    # all subexpressions of x which multiplied with y give z
    expr = 0

    def __init__(self, t, with_expr_tree=False):
        self.type = t

        self.roots = t.rootsystem
        self.load_words()
        self.load_perms()

        if with_expr_tree:
            self.load_expr_tree()

    # def __del__(self):
    # file where the perms are stored
    # needs to be closed when object dies
    # self.perm.file.close()
    # file where the words are stored
    # needs to be closed when object dies
    # self.word_str.file.close()
    # file where the expr_tree is stored
    # needs to be closed (if opened) when object dies
    # if self.expr_tree !=0:
    #    self.expr_tree.file.close()

    def load_perms(self):
        fname = "data/" + self.type.label + "/w.perm.hdf5"
        print("Loading:", fname)

        if os.path.isfile(fname):
            f = h5py.File(fname, 'r')
            self.perm = f["perm"]
            self.nr_perm = f["nr_perm"].value
            self.perm_ord = f["perm_ord"]
            self.perm_inv = f["perm_inv"]
        else:
            print("File not found: " + fname)
            print("Generating permutations from words...")
            self.generate_perms_from_words()
            self.load_perms()

    def load_words(self):
        fname = "data/" + self.type.label + "/w.word.hdf5"
        print("Loading:", fname)

        if os.path.isfile(fname):
            f = h5py.File(fname, 'r')
            self.word_str = f["word_str"]
            self.word_str_len = f["word_str_len"].value
            self.word_len = f["word_len"]
            self.refs = f["refs"]
            # longest element w0
            self.w0 = f["w0"].value
        else:
            print("File not found: " + fname)

    """
    Loading the tree for distinguished expressions
    """

    def load_expr_tree(self):
        fname = "data/" + self.type.label + "/w.expr_tree.hdf5"
        print("Loading:", fname)

        if os.path.isfile(fname):
            f = h5py.File(fname, 'r')
            self.expr_tree = f["expr_tree"]
        else:
            print("File not found: " + fname)
            print("Generating expression tree...")
            self.generate_expr_tree()
            self.load_expr_tree()

    """
    Loading self.expr from file ...dist_expr.pickle
    """

    def load_expr(self):
        fname = "data/" + self.type.label + "/w.expr.pickle"
        print("Loading:", fname)

        if os.path.isfile(fname):
            f = open(fname, 'rb')
            self.expr = pickle.load(f)
            f.close()
        else:
            print("File not found: " + fname)
            print("Generating expressions...")
            self.generate_expr()
            self.load_expr()

    """
    Loading self.dist_expr (all distinguished expressions) from file ...dist_expr.pickle
    """

    def load_dist_expr(self):
        fname = "data/" + self.type.label + "/w.dist_expr.pickle"
        print("Loading:", fname)

        if os.path.isfile(fname):
            f = open(fname, 'rb')
            self.dist_expr = pickle.load(f)
            f.close()
        else:
            print("File not found: " + fname)
            print("Generating distinguished expressions...")
            self.generate_dist_expr()
            self.load_dist_expr()

    """
    Finds the position of a permutation given as list
    """

    def pindex(self, lis):
        f = self.perm.file
        lis = str(lis)
        if lis in f["index"].keys():
            return f["index"][lis].value
        return -1

    """
    Finds the position of a permutation given as word
    """

    def windex(self, w):
        f = self.word_str.file
        # form="{:<"+str(self.word_str_len)+"}"
        # w is given as a numpy string so we have to remove the prefix
        # w=w[1:len(w)]
        # w=form.format(w)
        # print(w)
        w = np.string_(str(w)).rjust(self.word_str_len)
        if w in f["index"].keys():
            return f["index"][w].value
        return -1

    """
    Return the i-th word, by converting it from string to list
    """

    def word(self, i):
        result = str(self.word_str[i])
        result = result[2:len(result) - 1]
        return eval(result)

    """
    Multiply two permutations together
    --- perm[i]*perm[j]: first apply j then i
    """

    def pprod(self, i, j):
        a = self.perm[i]
        b = self.perm[j]
        return self.pindex(a[b])

    """
    Multiply reflections in word together to obtain index of perm
    --- here we use self.refs
    --- normally(why?) only product of simple reflections are used,
    --- BUT any reflection may be in the list
    """

    def wprod(self, w):
        if not w:
            # the index of the identity
            return 0
        result = self.refs[w[0]]
        lung = len(w)
        s = 1
        while s < lung:
            result = self.pprod(result, self.refs[w[s]])
            s = s + 1
        return result

    """
    Determine set of positive roots which are made negative by perm i
    """

    def plen_set(self, i):
        a = self.perm[i]
        npr = self.roots.nr_pos_roots
        result = []
        # we only look at positive roots
        for r in range(npr):
            if a[r] >= npr:
                result.append(r)
        return result

    """
    Generate the entry for:
    --- expr_tree[x][y][z] = [[A,xx,zz],[-1,-1,-1]] or [[B,xx,z],[C,xx,zz]]
    """

    def expr_tree_entry(self, x, y, z):
        if x == 0:
            if y == z:
                # empty expression is distinguished
                return [[-1, -1, -1], [-1, -1, -1]]
            else:
                # there is no distinguished expression
                return [[-2, -2, -2], [-2, -2, -2]]

        lung = self.word_len
        refs = self.refs

        # word of x
        xw = self.word(x)
        # word of x'
        xxw = xw[1:len(xw)]
        # index of x'
        xx = self.windex(xxw)
        if xx == -1:
            print("expr_tree_entry: this should not be!")
        # index of s_x[0]*z
        sz = self.pprod(refs[xw[0]], z)

        if lung[z] < lung[sz]:
            result = [[xw[0], xx, sz], [-1, -1, -1]]
        else:
            result = [[-1, xx, z], [xw[0], xx, sz]]
        return result

    """
    Generates the tree representation for dist expr
    --- creates a hdf5 file
    --- expr_tree[x][y][z] = [[A,xx,zz],[-1,-1,-1]] or [[B,xx,z],[C,xx,zz]]
    """

    def generate_expr_tree(self):
        fname = "data/" + self.type.label + "/w.expr_tree.hdf5"

        f = h5py.File(fname, 'w')
        lis = self.nr_perm
        dt = np.zeros((lis, lis, lis, 2, 3), np.int32)
        for x in range(lis):
            for y in range(lis):
                for z in range(lis):
                    dt[x][y][z] = self.expr_tree_entry(x, y, z)
        f["expr_tree"] = dt
        f.close()
        # self.load_expr_tree()

    """
    Get expressions from expr_tree
    """

    def expr_in_tree(self, x, y, z):
        result = []
        if x == 0:
            if y == z:
                # empty expression is distinguished
                # and A, B, C are empty
                return [[[], [[], [], []]]]
            else:
                # there is no distinguished expression
                return []
        dt = self.expr_tree[x][y][z]
        # A -> dt[0]
        # the second term in dt is [-1,-1,-1]
        if dt[1][1] == -1:
            tmp = self.expr_in_tree(dt[0][1], y, dt[0][2])
            for d in tmp:
                d[0].insert(0, dt[0][0])
                d[1][0] = [i + 1 for i in d[1][0]] + [0]
                d[1][1] = [i + 1 for i in d[1][1]]
                d[1][2] = [i + 1 for i in d[1][2]]
                result.append(d)
        else:
            # B -> dt[0]
            tmp = self.expr_in_tree(dt[0][1], y, dt[0][2])
            for d in tmp:
                d[0].insert(0, dt[0][0])
                d[1][0] = [i + 1 for i in d[1][0]]
                d[1][1] = [i + 1 for i in d[1][1]] + [0]
                d[1][2] = [i + 1 for i in d[1][2]]
                result.append(d)
            # C -> dt[1]
            tmp = self.expr_in_tree(dt[1][1], y, dt[1][2])
            for d in tmp:
                d[0].insert(0, dt[1][0])
                d[1][0] = [i + 1 for i in d[1][0]]
                d[1][1] = [i + 1 for i in d[1][1]]
                d[1][2] = [i + 1 for i in d[1][2]] + [0]
                result.append(d)

        return result

    """
    Test to see if e is a distinguished expression in J(x,y)
    """

    def is_dist_expr(self, e, x, y):
        lung = self.word_len
        xw = self.word(x)
        if len(e) != len(xw):
            return False
        for i in range(len(xw)):
            if e[i] == -1:
                xx = [k for k in e[i + 1:len(e)] if k != -1]
                xxperm = self.wprod(xx)
                xxpermy = self.pprod(xxperm, y)
                sxxpermy = self.pprod(self.refs[xw[i]], xxpermy)
                if lung[sxxpermy] > lung[xxpermy]:
                    return False
        return True

    """
    Get expressions from expr_tree
    """

    def dist_expr_in_tree(self, x, y, z):
        expr = self.expr_in_tree(x, y, z)
        result = []
        for c in expr:
            if self.is_dist_expr(c[0], x, y):
                result.append(c)
        return result

    """
    Collect all expr_in_tree(x,y,z) which are not zero
    in the file ...dist_expr.hdf5
    """

    def generate_dist_expr(self):
        fname = "data/" + self.type.label + "/w.dist_expr.pickle"
        # ---
        result = []
        lis = self.nr_perm
        for x in range(lis):
            for y in range(lis):
                for z in range(lis):
                    tmp = self.dist_expr_in_tree(x, y, z)
                    if tmp:
                        result.append([[x, y, z], tmp])
        # ---
        f = open(fname, 'wb')
        pickle.dump(result, f)
        f.close()

    """
    Collect all expr_in_tree(x,y,z) which are not zero
    in the file ...expr.hdf5
    """

    def generate_expr(self):
        fname = "data/" + self.type.label + "/w.expr.pickle"
        # ---
        result = []
        lis = self.nr_perm
        for x in range(lis):
            for y in range(lis):
                for z in range(lis):
                    tmp = self.expr_in_tree(x, y, z)
                    if tmp:
                        result.append([[x, y, z], tmp])
        # ---
        f = open(fname, 'wb')
        pickle.dump(result, f)
        f.close()

    """
    Generate image of permutation on roots of simple reflection r
    r is an element in [1..rang], so 1+npr is not treated
    """

    def generate_sref_r(self, r):
        if r >= self.type.rank:
            print("generate_sref_r: r is not a root")
            return False

        nr = self.roots.nr_roots
        npr = self.roots.nr_pos_roots
        chain = self.roots.chain
        result = [-1] * nr
        for i in range(nr):
            if i == r:
                result[i] = r + npr
            elif i == r + npr:
                result[i] = r
            elif result[i] == -1:
                c = chain(r, i)
                lung = len(c)
                for j in range(int((lung + 1) / 2)):
                    result[c[j]] = c[lung - j - 1]
                    result[c[lung - j - 1]] = c[j]

        return result

    """
    Generate image of permutations perm from words
    """

    def generate_perms_from_words(self):
        self.nr_perm = len(self.word_len)
        sref = [self.generate_sref_r(i) for i in range(self.type.rank)]

        # add id-element as first element
        # result=[np.array([i for i in range(self.roots.nr_roots)])]
        result = [[i for i in range(self.roots.nr_roots)]]
        for i in range(1, self.nr_perm):
            w = self.word(i)

            tmp = np.array(sref[w[0]])
            j = 1
            while j < len(w):
                tmp = tmp[sref[w[j]]]
                j = j + 1

            result.append(list(tmp))

        po = []
        for i in result:
            po += [self.determine_perm_order(i)]

        pi = []
        for i in result:
            pi += [np.array(result.index(self.determine_inverse_perm(i)))]

        tmp = []
        for i in result:
            tmp += [np.array(i)]
        result = tmp

        # self.perm.file.close()
        fname = "data/" + self.type.label + "/w.perm.hdf5"
        fname = fname.lower()
        f = h5py.File(fname, 'w')
        f["perm"] = result
        f["nr_perm"] = self.nr_perm
        f["perm_ord"] = po
        f["perm_inv"] = pi
        g = f.create_group("index")
        for i in range(self.nr_perm):
            g[str(result[i])] = i
        f.close()

    @staticmethod
    def determine_perm_order(p):
        tmp = list(p)
        one = list(range(len(p)))
        o = 1
        while tmp != one:
            tmp = [tmp[i] for i in p]
            o += 1
        return o

    @staticmethod
    def determine_inverse_perm(p):
        result = list(p)
        for i in range(len(p)):
            result[p[i]] = i
        return result
