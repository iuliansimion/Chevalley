class Bruhat:
    type = 0
    # underlying roots
    roots = 0
    # underlying weyl group
    weyl = 0

    def __init__(self, t, with_weyl_group=True):
        self.type = t

        self.roots = t.rootsystem
        if with_weyl_group:
            self.weyl = t.weyl

    """
    Test to see that an element is in Bruhat form
    --- U_w*w*T*U
    --- the expression of w needs to be in self.weyl.words
    """

    def in_Bruhat_form(self, E, n_coef=1, check_Uw=True):
        if not E:
            return True
        i = 0
        while i < len(E) and E[i][0] == "u":
            i = i + 1
        j = i
        while j < len(E) and E[j][0] == "n":
            j = j + 1
        k = j
        while k < len(E) and E[k][0] in ["t", "h"]:
            k = k + 1
        l = k
        while l < len(E) and E[l][0] == "u":
            l = l + 1

        EUw = E[0:i]
        EN = E[i:j]
        EH = E[j:k]
        EU = E[k:l]

        if not EN:
            if EH:
                if EUw:
                    print("group.in_Bruhat_form: (1) elements are not grouped!")
                    return False
            else:
                EU = EUw
                EUw = []

        # print(">> >> >>",EUw,EN,EH,EU)

        if len(EUw) + len(EN) + len(EH) + len(EU) != len(E):
            print("group.in_Bruhat_form: (2) elements are not grouped!", len(EUw), len(EN), len(EH), len(EU), len(E))
            return False

        w = []
        for n in EN:
            if n[2] != n_coef:
                print("test_Bruhat_form: not n_a(", n_coef, ")!")
                return False
            else:
                w.append(n[1])

        w = self.weyl.wprod(w)
        if w == -1:
            print("group.in_Bruhat_form: w not a known word!")
            return False
        #
        # U_w aici conjugam cu inversul
        # da, dar fiind pe dreapta, este deja in ordine inversa
        #
        if check_Uw:
            pwi = self.weyl.perm[self.weyl.perm_inv[w]]
            npr = self.roots.nr_pos_roots
            for u in EUw:
                if pwi[u[1]] < npr:
                    print("group.in_Bruhat_form: Uw not alright! w=", w, " r=", u[1])
                    return False

        return True

    """
    Separates L into three parts:
    --- the first consecutive u
    --- middle
    --- the last consecutive u

    --- IF L is unipotent, then this will be considered the last consecutive u
    """

    def split_Bruhat(self, L):
        # if L==[]:
        #    return [[],[],[]]

        L = list(L)
        npr = self.roots.nr_pos_roots
        i = len(L) - 1
        while i >= 0 and L[i][0] == "u" and L[i][1] < npr:  # move on positive roots only
            i = i - 1
        if i == -1:
            return [[], [], L]

        j = 0
        while j < len(L) and L[j][0] == "u" and L[j][1] < npr:  # move on positive roots only
            j = j + 1

        return [L[0:j], L[j:i + 1], L[i + 1:len(L)]]

    """
    Separates L into Bruhat parts:
    --- the first consecutive u
    --- n
    --- t
    --- the last consecutive u

    --- IF L is unipotent, then this will be considered the last consecutive u
    """

    def split_strict_Bruhat(self, L, n_coef=1):

        if not self.in_Bruhat_form(L, n_coef):
            print("group.split_strict_Bruhat: This should not be!")
            return False

        L = list(L)

        i = len(L) - 1
        while i >= 0 and L[i][0] == "u":
            i = i - 1
        if i == -1:
            return [[], [], [], L]

        j = 0
        while j < len(L) and L[j][0] == "u":
            j = j + 1

        k = j
        while k < len(L) and L[k][0] == "n":
            k = k + 1

        return [L[0:j], L[j:k], L[k:i + 1], L[i + 1:len(L)]]
