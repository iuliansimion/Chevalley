import sympy
import curtis


class CurtisG2(curtis.Curtis):
    # variables for the torus
    v = []
    # "primed" variables
    vv = []
    # "double primed" variables
    vvv = []

    def __init__(self, t):
        super().__init__(t)
        cmd = self.group.var.com_variables_l("v", ["a0", "b0", "c0", "d0"])
        exec(cmd[0])
        exec(cmd[1])

        v = self.v
        self.tori = [
            [["t", 0, v[0]], ["t", 1, v[1]]],
            [["t", 1, v[2]]],
            [],
            [["t", 0, v[3]]]
        ]

        cmd = self.group.var.com_variables_l("vv", ["a1", "b1", "c1", "d1"])
        exec(cmd[0])
        exec(cmd[1])

        vv = self.vv
        self.tori2 = [
            [["t", 0, vv[0]], ["t", 1, vv[1]]],
            [["t", 1, vv[2]]],
            [],
            [["t", 0, vv[3]]]
        ]

        cmd = self.group.var.com_variables_l("vvv", ["a2", "b2", "c2", "d2"])
        exec(cmd[0])
        exec(cmd[1])

        vvv = self.vvv
        self.tori3 = [
            [["t", 0, vvv[0]], ["t", 1, vvv[1]]],
            [["t", 1, vvv[2]]],
            [],
            [["t", 0, vvv[3]]]
        ]
