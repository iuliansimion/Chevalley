import sympy


class Variables:
    # lists of sympy variables
    x = []
    y = []
    z = []

    def __init__(self, nr_x, nr_y, nr_z):
        self.load_variables("x", nr_x)
        self.load_variables("y", nr_y)
        self.load_variables("z", nr_z)

    """
    This generates <code> to be used with exec(<code>)
    in order to define global indeterminates
    --- number a prefix
    """

    @staticmethod
    def com_variables(pref, sym_nr):
        sym_name = [pref + str(i) + ", " for i in range(1, sym_nr)]
        sym_name.append(pref + str(sym_nr))
        sym_name = "".join(sym_name)

        com1 = "(" + sym_name + ")=sympy.symbols(\"" + sym_name + "\")"
        com2 = "self." + pref + "=[" + sym_name + "]"

        return [com1, com2]

    """
    This generates <code> to be used with exec(<code>)
    in order to define global indeterminates
    --- uses a list of names
    """

    @staticmethod
    def com_variables_l(pref, names):
        sym_name = [names[i] + ", " for i in range(len(names) - 1)]
        sym_name.append(names[len(names) - 1])
        sym_name = "".join(sym_name)

        com1 = "(" + sym_name + ")=sympy.symbols(\"" + sym_name + "\")"
        com2 = "self." + pref + "=[" + sym_name + "]"

        return [com1, com2]

    """
    Define global indeterminates
    """

    def load_variables(self, pref, sym_nr):
        [com1, com2] = self.com_variables(pref, sym_nr)

        exec(com1)
        exec(com2)

    """
    Give Latex form of variable
    """

    def latex(self, v):
        xx = self.x
        yy = self.y
        zz = self.z
        result = str(v)
        for i in range(len(xx)):
            result = result.replace(str(xx[i]), "x_{" + str(i + 1) + "}")
        for i in range(len(yy)):
            result = result.replace(str(yy[i]), "y_{" + str(i + 1) + "}")
        for i in range(len(zz)):
            result = result.replace(str(zz[i]), "z_{" + str(i + 1) + "}")
        for i in range(3):
            for s in ["a", "b", "c", "d"]:
                result = result.replace(s + str(i), s + "_{" + str(i + 1) + "}")

        result = result.replace("**", "^")
        result = result.replace("*", "")
        return result
