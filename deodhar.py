class Deodhar:
    type = 0
    # the Chevalley group
    group = 0
    # the Weyl group
    weyl = 0
    # bruhat forms
    bruhat = 0

    def __init__(self, t):
        self.type = t

        self.weyl = t.weyl
        self.group = t.group
        self.bruhat = t.bruhat
        if self.type.weyl.dist_expr == 0:
            self.type.loadDistinguishedExpressions()

    """
    Construct cell given an expression

    eg. cell_Di(word(14),dist_expr[14][15][9][0][1])
    de=tmp.group.weyl.dist_expr[54]
    tmp.cell_D(tmp.group.weyl.word(de[0][0]),de[1][0][1],xx)
    """

    # def cell_D(self,w,ABC,xx):
    def cell_D(self, i, j, xx=0):
        de = self.weyl.dist_expr[i]
        w = self.weyl.word(de[0][0])
        ABC = de[1][j][1]
        if xx == 0:
            xx = self.group.var.x
        npr = self.group.roots.nr_pos_roots
        [A, B, C] = ABC
        result = []
        ll = len(w)
        for r in range(ll):
            if r in A:
                result.append(["u", w[r], xx[r]])
                result.append(["n", w[r], 1])
            elif r in B:
                result.append(["u", w[r] + npr, xx[r]])
                # result.append(["u",w[r]+npr,xx[r]**-1])
                # result.append(["u",w[r],xx[r]])
                # result.append(["n",w[r],1])
            elif r in C:
                result.append(["n", w[r], 1])
            else:
                print(">>>>>>", w, A, B, C)
                print("cell_Di: this should not be!")
                return False
        return result

    """
    Construct cell given an expression
    --- HERE I WRITE U_{-a}(t^-1) as U_{a}(-t)N_{a}(-t)U_{a}(-t)

    eg. cell_D2i(word(14),dist_expr[14][15][9][0][1])
    de=tmp.group.weyl.dist_expr[54]
    tmp.cell_DI(tmp.group.weyl.word(de[0][0]),de[1][0][1],xx)
    """

    def cell_DI(self, i, j, xx=0):
        de = self.weyl.dist_expr[i]
        w = self.weyl.word(de[0][0])
        ABC = de[1][j][1]
        if xx == 0:
            xx = self.group.var.x
        [A, B, C] = ABC
        result = []
        ll = len(w)
        for r in range(ll):
            if r in A:
                result.append(["u", w[r], xx[r]])
                result.append(["n", w[r], 1])
            elif r in B:
                result.append(["u", w[r], xx[r] ** -1])
                result.append(["n", w[r], -xx[r] ** -1])
                result.append(["u", w[r], xx[r] ** -1])
            elif r in C:
                result.append(["n", w[r], 1])
            else:
                print(">>>>>>", w, A, B, C)
                print("cell_Di: this should not be!")
                return False
        return result

    """
    Move the u atoms on the left or on the right such that DD=UxU
    !!! if possible
    """

    def move_u_left_right(self, DD):
        # print("deodhar.move_u_left_right")
        perm = self.group.weyl.perm
        npr = self.group.roots.nr_pos_roots
        modif = True
        # print(">>>.>>> DD=",DD)
        while modif:
            modif = False
            i = 0
            while i < len(DD) and DD[i][0] == "u" and DD[i][1] < npr:
                i = i + 1
            # print("deodhar.move_u_left_right: before splitting in Bruhat form!")
            DDD = self.bruhat.split_Bruhat(DD)
            # print(">>>...>>> DDD=",i,DDD)
            # print("deodhar.move_u_left_right: splitted in Bruhat form!")
            if not self.group.all_thn(DDD[1]) or \
                    not self.group.all_up(DDD[0]) or \
                    not self.group.all_up(DDD[2]):
                while i < len(DD) and not modif:
                    if DD[i][0] == "u":
                        # print(">>>.>>>",i,DD)
                        # print(">>>.>>>",i,self.group.latex(DD))
                        w_right = self.group.e_to_w(DD[i + 1:len(DD)])
                        w_right = self.group.weyl.wprod(w_right)
                        w_right = self.group.weyl.perm_inv[w_right]
                        # print("w_right:",w_right)
                        if perm[w_right][DD[i][1]] < npr:
                            DD = DD[0:i] + self.group.move_atom_right(DD[i], DD[i + 1:len(DD)])
                            DD = self.group.canonic_pos_u_part(DD)
                            # DD=self.group.simplify_u(DD)
                            DD = self.group.simplify_params(DD)
                            modif = True
                        else:
                            w_left = self.group.e_to_w(DD[0:i])
                            w_left = self.group.weyl.wprod(w_left)
                            # print("w_left:",w_left)
                            if perm[w_left][DD[i][1]] < npr:
                                DD = self.group.move_atom_left(DD[0:i], DD[i]) + DD[i + 1:len(DD)]
                                DD = self.group.canonic_pos_u_part(DD)
                                # DD=self.group.simplify_u(DD)
                                DD = self.group.simplify_params(DD)
                                modif = True
                            else:
                                # print(">>>..>>>",i,DD)
                                # print(">>>>>de",de)
                                print("..deodhar.move_u_left_right: This should not be (root gorup not movable)!")
                                return False

                    i = i + 1
                # print(">>>.end",i,self.group.latex(DD))

        DDD = self.bruhat.split_Bruhat(DD)
        if not self.group.all_thn(DDD[1]):
            print("deodhar.move_u_left_right: This should not be (thn)!")

        return DD

    """
    Bring cel D to zUy^-1-form, i.e. determine the U(y^-1)T
    """

    def cell_UyiU(self, i, j):
        DD = self.cell_D(i, j)
        de = self.weyl.dist_expr[i]
        #
        # x = de[0][0]
        # y = de[0][1]
        z = de[0][2]

        # z^-1
        zz = self.group.w_to_ni(self.group.weyl.word(z))
        # if I use self.group.invert() it will convert n(1) to n(-1)
        zz.reverse()
        DD = zz + DD

        DD = self.move_u_left_right(DD)
        DDD = self.bruhat.split_Bruhat(DD)

        DDD = self.group.canonic_u(DDD[0]) + self.group.canonic_nih(DDD[1]) + self.group.canonic_u(DDD[2])

        return DDD

    """
    Bring cel D to Ux-form, i.e. move the u on the left
    """

    def cell_Ux(self, i, j):
        DD = self.cell_DI(i, j)
        DD = self.move_u_left_right(DD)
        DDD = self.bruhat.split_Bruhat(DD)

        DDD = self.group.canonic_u(DDD[0]) + self.group.canonic_nh(DDD[1]) + self.group.canonic_u(DDD[2])

        return DDD
