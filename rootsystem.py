import os.path
import json
import numpy as np


class RootSystem:
    type = 0

    label = "none"

    # list of roots as linear combination of simple roots:
    roots = []
    # number of roots
    nr_roots = 0
    # number of positive roots
    nr_pos_roots = 0
    # matrix (r,s) -> position of r+s
    rsum = []
    # array indicating if root is short or long
    short = []
    # square of ratio of long root length over short root length
    ls_square = 0
    # extra special pairs
    es_pairs = []

    def __init__(self, t):
        self.type = t

        self.load_roots()
        self.load_roots_sums()
        self.load_short_roots()

        if self.type.isSimpleLaced():
            self.ls_square = 1
        elif self.type.isDoubleLaced():
            self.ls_square = 2
        elif self.type.isTripleLaced():
            self.ls_square = 3
        else:
            print("Unknown Lie type!")

    def load_roots(self):
        fname = "data/" + self.type.label + "/roots.json"
        print("Loading:", fname)

        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                self.roots = json.load(f)

            self.nr_roots = len(self.roots)
            self.nr_pos_roots = int(self.nr_roots / 2)

            # convert roots to numpy.array
            for i in range(self.nr_roots):
                self.roots[i] = np.array(self.roots[i])
        else:
            print("File not found: " + fname)

    def load_roots_sums(self):
        fname = "data/" + self.type.label + "/roots_sum.json"
        print("Loading:", fname)

        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                self.rsum = json.load(f)
        else:
            print("File not found: " + fname)
            print("Generating root sums...")
            self.generate_root_sums()
            print("Exporting root sums to: " + fname)
            self.export_root_sums()

    def generate_root_sums(self):
        self.rsum = [[0] * self.nr_roots] * self.nr_roots
        for r in range(self.nr_roots):
            for s in range(self.nr_roots):
                self.rsum[r][s] = self.index(self.roots[r] + self.roots[s])

    def export_root_sums(self, fname=0):
        if fname == 0:
            fname = "data/" + self.type.label + "/roots_sum.json"
        with open(fname, 'w') as f:
            json.dump(self.rsum, f)

    def load_short_roots(self):
        fname = "data/" + self.type.label + "/short_roots.json"
        print("Loading:", fname)

        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                self.short = json.load(f)
        else:
            print("File not found: " + fname)
            print("Determining short/long roots...")
            self.short = self.determine_short_roots()
            print("Exporting short/long roots to: " + fname)
            self.export_short_roots()

    def export_short_roots(self, fname=0):
        if fname == 0:
            fname = "data/" + self.type.label + "/short_roots.json"
        with open(fname, 'w') as f:
            json.dump(self.short, f)

    """
    Get index of a root expressed as np.array r=[0,1,0..] in roots
    """

    def index(self, r):
        result = [np.array_equal(r, x) for x in self.roots]
        if True in result:
            return result.index(True)
        return -1

    """
    Determine the length of an r-chain through s
    """

    def chain_len(self, r, s):
        rr = self.roots
        p = 0
        while self.index(-p * rr[r] + rr[s]) != -1:
            p = p + 1
        q = 0
        while self.index(q * rr[r] + rr[s]) != -1:
            q = q + 1
        return q + p - 2

    """
    Determine the length of the left part of the r-chain through s
    p in -pr+s, ...,-r+s,s
    """

    def chain_len_left(self, r, s):
        rr = self.roots
        p = 0
        while self.index(-p * rr[r] + rr[s]) != -1:
            p = p + 1
        return p - 1

    """
    Determine the left most root in the r-chain through s
    [-pr+s], ...,-r+s,s
    """

    def chain_left_root(self, r, s):
        rr = self.roots
        p = 0
        while self.index(-p * rr[r] + rr[s]) != -1:
            p = p + 1
        p = p - 1
        return self.index(-p * rr[r] + rr[s])

    """
    Determine the r-chain through s (in order from left to right)
    -pr+s, ...,-r+s,s
    """

    def chain(self, r, s):
        rsum = self.rsum
        prs = self.chain_left_root(r, s)
        result = [prs]
        while rsum[r][prs] != -1:
            prs = rsum[r][prs]
            result.append(prs)
        return result

    """
    Determine the length of the right part of the r-chain through s
    q in s, s+r, s+2r, ..., s+qr
    """

    def chain_len_right(self, r, s):
        rr = self.roots
        q = 0
        while self.index(q * rr[r] + rr[s]) != -1:
            q = q + 1
        return q - 1

    """
    Determine the position i of simple roots in roots sum(roots[i])=1
    """

    def simple_roots(self):
        result = []
        for i in range(self.nr_pos_roots):
            if sum(self.roots[i]) == 1:
                result.append(i)
        return result

    """
    Method for determining which simple roots are short
    from the decomposition in terms of simple roots
    """

    def determine_short_simple_roots(self):
        if self.type.isSimpleLaced():
            result = [True] * self.type.rank
        else:
            sr = self.simple_roots()
            result = [-1] * self.type.rank
            while -1 in result:
                for r in sr:
                    for s in sr:
                        v = self.rsum[r][s]
                        if v != -1:
                            # if I can append another r,
                            # then r is short and s is long
                            if self.rsum[v][r] != -1:
                                result[r] = True
                                result[s] = False
                            # if I can append another s,
                            # then r is long and s is short
                            elif self.rsum[v][s] != -1:
                                result[r] = False
                                result[s] = True
                            # if I cannot append r or s
                            # but I know one, then I know the other
                            elif result[r] != -1:
                                result[s] = result[r]
                            elif result[s] != -1:
                                result[r] = result[s]
        return result

    """
    Method for determining which roots are short and which roots are long
    """

    def determine_short_roots(self):
        if self.type.isSimpleLaced():
            result = [True] * self.nr_roots
        elif self.type.isTripleLaced():
            result = self.determine_short_roots_g2()
        else:
            # we are in doubly laced type
            result = [-1] * self.nr_roots
            ssr = self.determine_short_simple_roots()
            for i in range(len(ssr)):
                result[i] = ssr[i]
            sr = self.simple_roots()
            modif = True
            while modif:
                modif = False
                for s in sr:
                    for r in range(self.nr_pos_roots):
                        if result[r] != -1:
                            v = self.rsum[s][r]
                            if v != -1 and result[v] == -1:
                                # long + long = long
                                if not result[s] and not result[r]:
                                    result[v] = False
                                # short + long = short
                                elif result[s] and not result[r]:
                                    result[v] = True
                                # long + short = short
                                elif not result[s] and result[r]:
                                    result[v] = True
                                # short + short
                                else:
                                    # = long if the s-chain > 1
                                    if self.chain_len(s, r) == 2:
                                        result[v] = False
                                    else:
                                        result[v] = True
                                modif = True
            for r in range(self.nr_pos_roots):
                result[r + self.nr_pos_roots] = result[r]

        return result

    """
    For G2, the method for determining which roots are 
    short and which roots are long is hardcoded
    """

    def determine_short_roots_g2(self):
        result = [-1] * self.nr_roots

        # determine the shortness of simple roots
        ssr = self.determine_short_simple_roots()
        s = ssr.index(True)
        r = ssr.index(False)
        result[s] = True
        result[r] = False

        # determine the shortness of the other roots
        rr = self.roots
        result[self.rsum[s][r]] = True
        result[self.index(2 * rr[s] + rr[r])] = True
        result[self.index(3 * rr[s] + rr[r])] = False
        result[self.index(3 * rr[s] + 2 * rr[r])] = False

        # determine negative roots
        result[s + self.nr_pos_roots] = True
        result[r + self.nr_pos_roots] = False
        result[self.index(-(rr[s] + rr[r]))] = True
        result[self.index(-(2 * rr[s] + rr[r]))] = True
        result[self.index(-(3 * rr[s] + rr[r]))] = False
        result[self.index(-(3 * rr[s] + 2 * rr[r]))] = False

        return result

    """
    Determine extra special pairs:
       - Carter p58. for extra special pairs
       - Carter p13. for total order of vector space
    """

    def extra_special_pairs(self):
        result = []
        for i in range(self.type.rank, self.nr_pos_roots):
            cauta = True
            j = 0
            while cauta and j < self.type.rank:
                v = self.rsum[i][j + self.nr_pos_roots]
                # find first non-zero coefficient
                if v != -1:
                    k = 0
                    while k < self.type.rank and self.roots[v][k] == 0:
                        k = k + 1
                    if k < self.type.rank and self.roots[v][k] > 0:
                        result.append([j, v])
                        cauta = False
                j = j + 1
        return result

    """
    Method for checking that the four roots are pair-wise
    not opposite
    """

    def no_pair_opp(self, r1, r2, r3, r4):
        pnr = self.nr_pos_roots
        opp = [r1 + pnr, r1 - pnr, r2 + pnr, r2 - pnr, r3 + pnr, r3 - pnr, r4 + pnr, r4 - pnr]
        if r1 in opp or r2 in opp or r3 in opp or r4 in opp:
            return False
        return True

    """
    Method for checking that the four roots are pair-wise
    distinct
    """

    @staticmethod
    def pair_distinct(r1, r2, r3, r4):
        if r1 in [r2, r3, r4] or r2 in [r3, r4] or r3 == r4:
            return False
        return True

    """
    Generate roots which are linear combinations
    with positive coefficients of given roots
    """

    def generate_roots(self, rr):
        rsum = self.rsum
        modif = True
        while modif:
            modif = False
            i = 0
            while i < len(rr):
                j = i + 1
                while j < len(rr):
                    poz = rsum[rr[i]][rr[j]]
                    if poz != -1 and poz not in rr:
                        rr.append(poz)
                        modif = True
                    j = j + 1
                i = i + 1
        return rr

    """
    Determine if all roots are in a half-plane
    """

    def in_half_plane(self, rr):
        npr = self.nr_pos_roots
        rr = self.generate_roots(rr)
        for i in rr:
            if i < npr and i + npr in rr:
                print("roots.in_half_plane: not in one half_plane!")
                return False
            if i >= npr and i - npr in rr:
                print("roots.in_half_plane: not in one half_plane!")
                return False
        return True

    """
    Determine if all roots in the list rr are positive
    """

    def all_positive(self, rr):
        npr = self.nr_pos_roots
        test = True
        i = 1
        while test and i < len(rr):
            if rr[i] >= npr:
                return False
            i = i + 1
        return True

    """
    Given a root r return the index of -r
    """

    def minus_r(self, r):
        npr = self.nr_pos_roots
        if r < npr:
            return r + npr
        else:
            return r - npr
