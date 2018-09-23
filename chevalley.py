import rootsystem as rs
import weylgroup as wg
import constants as cs
import parabolics as pa
import variables as vr
import group as gr
import deodhar as dr
import liealgebra as la
import adaction as ad
import bruhat as br
from curtisa2 import CurtisA2
from curtisb2 import CurtisB2
from curtisg2 import CurtisG2


class LieTypes(object):
    existing = {}

    # Create based on class name:
    @staticmethod
    def get_lie_type(family, rank):
        lie_type: str = family + str(rank)
        if lie_type in LieTypes.existing.keys():
            return LieTypes.existing[lie_type]
        else:
            LieTypes.existing[lie_type] = LieType(family, rank)
            return LieTypes.existing[lie_type]

    factory = staticmethod(get_lie_type)


class LieType(object):
    # Type and rank
    lie_type = ""
    rank = 0
    label = ""
    dim = 0
    # the corresponding root system
    rootsystem = 0
    # the corresponding weyl group
    weyl = 0
    # parabolic subgroups
    parabolics = 0
    # structure constants of lie algebra and group
    constants = 0
    # indeterminates/variables to calculate with
    var = 0
    # the chevalley group (fixed points of adjoint group)
    group = 0
    # the lie algebra
    lie = 0
    # adjoint action
    ad_action = 0

    # specific
    bruhat = 0
    deodhar = 0
    curtis = 0

    def __init__(self, lie_type, rank):
        self.lie_type = lie_type
        self.rank = rank
        self.label = self.lie_type + str(self.rank)

        self.rootsystem = rs.RootSystem(self)
        self.dim = self.rootsystem.nr_roots + self.rank

        self.weyl = wg.Weyl(self)
        self.parabolics = pa.Parabolics(self)
        self.constants = cs.Constants(self)
        nr = self.rootsystem.nr_roots
        self.var = vr.Variables(2 * nr, 2 * nr, 2 * self.rank)
        self.group = gr.Group(self)
        self.lie = la.LieAlgebra(self)
        self.ad_action = ad.AdAction(self)

        self.bruhat = br.Bruhat(self)
        self.deodhar = dr.Deodhar(self)
        if self.label == "a2":
            self.curtis = CurtisA2(self)
        if self.label == "b2":
            self.curtis = CurtisB2(self)
        if self.label == "g2":
            self.curtis = CurtisG2(self)

    def isSimpleLaced(self):
        return self.lie_type in ["a", "d", "e"]

    def isDoubleLaced(self):
        return self.lie_type in ["b", "c", "f"]

    def isTripleLaced(self):
        return self.lie_type in ["g"]

    def loadDistinguishedExpressions(self):
        self.weyl.load_expr_tree()
        self.weyl.load_expr()
        self.weyl.load_dist_expr()
