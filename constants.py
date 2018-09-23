import os.path
import json


class Constants:
    type = 0

    # a link to the underlying roots
    roots = 0

    #
    # structure constants Nrs, Ars, Mrsi, Cijrs, etars
    #
    N = []
    A = []
    M = []
    C = []
    eta = []

    def __init__(self, t):  # tip,rang,with_roots=True):
        self.type = t

        self.roots = t.rootsystem

        self.load_N()
        self.load_A()
        self.load_M()
        self.load_C()
        self.load_eta()

    """
    Function for loading Nrs, Ars, Mrsi, Cijrs, etars (as in Carter)
    """

    def load_N(self):
        fname = "data/" + self.type.label + "/n.json"
        print("Loading:", fname)

        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                self.N = json.load(f)
        else:
            print("File not found: " + fname)
            print("Generating N...")
            self.N = self.Nrs()
            print("Exporting N to: " + fname)
            self.export_N()

    def load_A(self):
        fname = "data/" + self.type.label + "/a.json"
        print("Loading:", fname)

        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                self.A = json.load(f)
        else:
            print("File not found: " + fname)
            print("Generating A...")
            self.A = self.Ars()
            print("Exporting A to: " + fname)
            self.export_A()

    def load_M(self):
        fname = "data/" + self.type.label + "/m.json"
        print("Loading:", fname)

        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                self.M = json.load(f)
        else:
            print("File not found: " + fname)
            print("Generating M...")
            self.M = self.Mrsi()
            print("Exporting M to: " + fname)
            self.export_M()

    def load_C(self):
        fname = "data/" + self.type.label + "/c.json"
        print("Loading:", fname)

        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                self.C = json.load(f)
        else:
            print("File not found: " + fname)
            print("Generating C...")
            self.C = self.Cijrs()
            print("Exporting C to: " + fname)
            self.export_C()

    def load_eta(self):
        fname = "data/" + self.type.label + "/eta.json"
        print("Loading:", fname)

        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                self.eta = json.load(f)
        else:
            print("File not found: " + fname)
            print("Generating eta...")
            self.eta = self.etars()
            print("Exporting eta to: " + fname)
            self.export_eta()

    def export_N(self, fname=0):
        if fname == 0:
            fname = "data/" + self.type.label + "/n.json"
        with open(fname, 'w') as f:
            json.dump(self.N, f)

    def export_A(self, fname=0):
        if fname == 0:
            fname = "data/" + self.type.label + "/a.json"
        with open(fname, 'w') as f:
            json.dump(self.A, f)

    def export_M(self, fname=0):
        if fname == 0:
            fname = "data/" + self.type.label + "/m.json"
        with open(fname, 'w') as f:
            json.dump(self.M, f)

    def export_C(self, fname=0):
        if fname == 0:
            fname = "data/" + self.type.label + "/c.json"
        with open(fname, 'w') as f:
            json.dump(self.C, f)

    def export_eta(self, fname=0):
        if fname == 0:
            fname = "data/" + self.type.label + "/eta.json"
        with open(fname, 'w') as f:
            json.dump(self.eta, f)

    """
    Determine absolute values of Nrs=+/-(p+1)
      - Carter p.55
    """

    def absNrs(self):
        result = []
        for r in range(self.roots.nr_roots):
            result.append([])
            for s in range(self.roots.nr_roots):
                if r == s:
                    result[r].append(False)
                elif self.roots.rsum[r][s] != -1:
                    result[r].append(self.roots.chain_len_left(r, s) + 1)
                else:
                    result[r].append(0)
        return result

    """
    Determine the signs of Nrs
      - Carter p.55 Theorem 4.1.2
    """

    def sgnNrs(self):
        nr = self.roots.nr_roots
        S = [[0] * nr] * nr

        # set sign on extra special pairs to +1
        esp = self.roots.extra_special_pairs()
        for e in esp:
            S[e[0]][e[1]] = 1

        # use the properties in Carter p.55 to deduce other signs
        absN = self.absNrs()
        modif = True
        while modif:
            [modif_i, S] = self.sgnNrs_fill_i(S)
            [modif_ii, S] = self.sgnNrs_fill_ii(S)
            [modif_iii, S] = self.sgnNrs_fill_iii(S)
            # modif_iv=False
            [modif_iv, S] = self.sgnNrs_fill_iv(S, absN)
            modif = modif_i or modif_ii or modif_iii or modif_iv

        return S

    """
    Carter p.55 Theorem 4.1.2 (i)
    """

    def sgnNrs_fill_i(self, S):
        rsum = self.roots.rsum
        nr = self.roots.nr_roots
        modif = False
        for r in range(nr):
            for s in range(nr):
                if rsum[r][s] != -1 and S[r][s] == 0 and S[s][r] != 0:
                    S[r][s] = -S[s][r]
                    modif = True
        return [modif, S]

    """
    Carter p.55 Theorem 4.1.2 (ii)
    """

    def sgnNrs_fill_ii(self, S):
        rsum = self.roots.rsum
        nr = self.roots.nr_roots
        pnr = self.roots.nr_pos_roots
        modif = False
        for r1 in range(nr):
            for r2 in range(nr):
                if rsum[r1][r2] != -1:
                    # determine -(r1+r2)
                    if rsum[r1][r2] >= pnr:
                        r3 = rsum[r1][r2] - pnr
                    else:
                        r3 = rsum[r1][r2] + pnr
                    # now r1+r2+r3=0

                    # we make a change only if one of the
                    # corresponding Srs is non-zero
                    tmp = S[r1][r2] * S[r2][r3] * S[r3][r1]
                    if tmp == 0:
                        if S[r1][r2] != 0:
                            S[r2][r3] = S[r1][r2]
                            S[r3][r1] = S[r1][r2]
                            modif = True
                        elif S[r2][r3] != 0:
                            S[r1][r2] = S[r2][r3]
                            S[r3][r1] = S[r2][r3]
                            modif = True
                        elif S[r3][r1] != 0:
                            S[r1][r2] = S[r3][r1]
                            S[r2][r3] = S[r3][r1]
                            modif = True
        return [modif, S]

    """
    Carter p.55 Theorem 4.1.2 (iii)
    """

    def sgnNrs_fill_iii(self, S):
        pnr = self.roots.nr_pos_roots
        modif = False
        for r in range(pnr):
            for s in range(pnr):
                # known r,s
                if S[r + pnr][s + pnr] == 0 and S[r][s] != 0:
                    S[r + pnr][s + pnr] = -S[r][s]
                # known r,-s
                if S[r + pnr][s] == 0 and S[r][s + pnr] != 0:
                    S[r + pnr][s] = -S[r][s + pnr]
                # known -r,s
                if S[r][s + pnr] == 0 and S[r + pnr][s] != 0:
                    S[r][s + pnr] = -S[r + pnr][s]
                # known -r,-s
                if S[r][s] == 0 and S[r + pnr][s + pnr] != 0:
                    S[r][s] = -S[r + pnr][s + pnr]

        return [modif, S]

    """
    Carter p.55 Theorem 4.1.2 (iv)
    """

    def sgnNrs_fill_iv(self, S, absN):
        rsum = self.roots.rsum
        nr = self.roots.nr_roots
        index = self.roots.index
        rr = self.roots.roots
        lss = self.roots.ls_square
        short = self.roots.short
        npo = self.roots.no_pair_opp
        pdi = self.roots.pair_distinct
        modif = False
        for r1 in range(nr):
            for r2 in range(nr):
                for r3 in range(nr):
                    r4 = index(-(rr[r1] + rr[r2] + rr[r3]))
                    if r4 != -1 and npo(r1, r2, r3, r3) and pdi(r1, r2, r3, r4):
                        # now r1+r2+r3+r4=0
                        cond = absN[r1][r2] != 0 and S[r1][r2] != 0
                        if cond:
                            cond = cond and absN[r3][r4] != 0 and S[r3][r4] != 0
                        if cond:
                            cond = cond and absN[r2][r3] != 0 and S[r2][r3] != 0
                        if cond:
                            cond = cond and absN[r1][r4] != 0 and S[r1][r4] != 0
                        if cond:
                            cond = cond and absN[r3][r1] != 0 and S[r3][r1] != 0
                        if cond:
                            cond = cond and absN[r2][r4] != 0 and S[r2][r4] == 0
                        # now we can deduce S[r2][r4]:
                        if cond:
                            if not short[rsum[r1][r2]]:
                                num1 = lss
                            else:
                                num1 = 1
                            if not short[rsum[r2][r3]]:
                                num2 = lss
                            else:
                                num2 = 1
                            if not short[rsum[r3][r1]]:
                                num3 = lss
                            else:
                                num3 = 1

                            s = 0
                            s = s + S[r1][r2] * absN[r1][r2] * S[r3][r4] * absN[r3][r4] / num1
                            s = s + S[r2][r3] * absN[r2][r3] * S[r1][r4] * absN[r1][r4] / num2
                            s = s + S[r3][r1] * absN[r3][r1] * absN[r2][r4] / num3
                            if s == 0:
                                S[r2][r4] = 1
                            else:
                                S[r2][r4] = -1

                            modif = True

        return [modif, S]

    """
    Method for calculating Nrs from absNrs() and sgnNrs()
    """

    def Nrs(self):
        absN = self.absNrs()
        S = self.sgnNrs()
        nr = self.roots.nr_roots
        result = []
        for r in range(nr):
            result.append([])
            for s in range(nr):
                result[r].append(S[r][s] * absN[r][s])
        return result

    """
    Method for calculating Ars
      - Carter p.38
      - pay attention to r=+/-s
    """

    def Ars(self):
        result = []
        pnr = self.roots.nr_pos_roots
        p = self.roots.chain_len_left
        q = self.roots.chain_len_right
        nr = self.roots.nr_roots
        for r in range(nr):
            result.append([])
            for s in range(nr):
                if r == s:
                    result[r].append(2)
                elif r in [s + pnr, s - pnr]:
                    result[r].append(-2)
                else:
                    result[r].append(p(r, s) - q(r, s))

        return result

    """
    Method for calculating Mrsi
      - Carter p.61
      - i=1,2,..
      - Mrs1=Nrs
    """

    def Mrsi(self):
        nr = self.roots.nr_roots
        index = self.roots.index
        rr = self.roots.roots
        result = []
        for r in range(nr):
            result.append([])
            for s in range(nr):
                result[r].append([])
                i = 1
                tmp1 = 1
                tmp2 = 1
                while index((i - 1) * rr[r] + rr[s]) != -1:
                    v = index((i - 1) * rr[r] + rr[s])
                    if self.N[r][v] == 0:
                        break
                    tmp1 = tmp1 * self.N[r][v]
                    tmp2 = tmp2 * i
                    result[r][s].append(int(tmp1 / tmp2))
                    i = i + 1
        return result

    """
    Method for calculating Cijrs
    - Carter p.77 Theorem 5.2.2
    - the order of the terms must be such that i+j is increasing
    - i,j=1,2,3
    """

    def Cijrs(self):
        rsum = self.roots.rsum
        M = self.M
        nr = self.roots.nr_roots
        result = []
        for r in range(nr):
            result.append([])
            for s in range(nr):
                result[r].append([])
                l1 = len(M[r][s])
                l2 = len(M[s][r])
                i = 0
                # we don't want Mrs11 twice
                j = 1
                while i < l1 or j < l2:
                    if i < l1:
                        # careful with the indexing
                        result[r][s].append([i + 1, 1, M[r][s][i]])
                        i = i + 1
                    if j < l2:
                        # careful with the exponent of -1
                        result[r][s].append([1, j + 1, ((-1) ** (j + 1)) * M[s][r][j]])
                        j = j + 1
                # for G2 there are extra terms
                if self.type.isTripleLaced():
                    v = rsum[r][s]
                    if v != -1:
                        if len(M[v][r]) > 1:
                            result[r][s].append([3, 2, int(M[v][r][1] / 3)])
                        if len(M[v][s]) > 1:
                            result[r][s].append([2, 3, int(-2 * M[v][s][1] / 3)])

        return result

    """
    Method for calculating eta_rs for a pair of roots
    - Carter pp.94-95 for eta
    - Carter p.84 for epsilon
    """

    def etars_pair(self, r, s):
        p = self.roots.chain_len_left(r, s)
        q = self.roots.chain_len_right(r, s)
        rr = self.roots.roots
        # the starting root for the r-chain through s
        ss = self.roots.index(-p * rr[r] + rr[s])
        epsilon = []
        for i in range(p + q):
            v = self.roots.index(i * rr[r] + rr[ss])
            if self.N[r][v] > 0:
                epsilon.append(1)
            elif self.N[r][v] < 0:
                epsilon.append(-1)
            else:
                print("etars_pair: Nrv==0 - This should not be!")
        tmp1 = 1
        for i in epsilon[0:p]:
            tmp1 = tmp1 * i
        tmp2 = 1
        # q on p.94 (Carter) is the length p+q of the chain
        for i in epsilon[0:q]:  # q-p-1]:
            tmp2 = tmp2 * i
        return int(((-1) ** p) * tmp1 / tmp2)

    """
    Method for calculating eta_rs
    - Carter pp.94-95 for eta
    - Carter p.84 for epsilon
    """

    def etars(self):
        nr = self.roots.nr_roots
        pnr = self.roots.nr_pos_roots
        result = []
        for r in range(nr):
            result.append([])
            for s in range(nr):
                if r in [s, s + pnr, s - pnr]:
                    result[r].append(-1)
                else:
                    result[r].append(self.etars_pair(r, s))
        return result
