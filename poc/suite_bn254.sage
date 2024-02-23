#!/usr/bin/sage
# vim: syntax=python

import hashlib
import sys
from hash_to_field import XMDExpander
try:
    from sagelib.common import test_dst
    from sagelib.h2c_suite import BasicH2CSuiteDef, BasicH2CSuite, IsoH2CSuiteDef, IsoH2CSuite
    from sagelib.svdw_generic import GenericSvdW
    from sagelib.sswu_generic import GenericSSWU
    from sagelib.sswu_optimized import OptimizedSSWU
    from sagelib.ft import FouqueTibouchi
    from sagelib.iso_values import iso_bn254
except ImportError:
    sys.exit("Error loading preprocessed sage files. Try running `make clean pyfiles`")

p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
F = GF(p)
A = F(0)
B = F(3)
E = EllipticCurve(GF(p), [A, B])
# curve isogenous to bn254 with degree 59
Ap = F(9087994317191712533568698403530528306233527979934880849865820425505218365052)
Bp = F(3059101143800926337153883959975852125336293569895750485959800095292563537400)
iso_map = iso_bn254()

if __name__ == "__main__":
    u = 7105195380181880595384217009108718366423089053558315283835256316808390512725
    print("u: ", u)
    # ===
    svdw = GenericSvdW(F, A, B)
    x1, y1 = svdw.straight_line(u)
    print("svdw:\n", x1, y1, "\n")
    # ===
    ft = FouqueTibouchi(F, A, B)
    x1, y1 = ft.map_to_point_ft(u)
    print("ft:\n", x1, y1, "\n")
    # ===
    sswu = GenericSSWU(F, Ap, Bp) # on isogenous curve E'
    x1, y1 = sswu.straight_line(u)
    x2, y2 = sswu.not_straight_line(u)
    assert x1 == x2 and y1 == y2, "straight line/non straight line mismatch"
    print("generic sswu:\n", iso_map((x1, y1)), "\n")
    # ===
    sswu = OptimizedSSWU(F, Ap, Bp) # on isogenous curve E'
    x1, y1 = sswu.straight_line(u)
    x2, y2 = sswu.not_straight_line(u)
    assert x1 == x2 and y1 == y2, "straight line/non straight line mismatch"
    print("optimised sswu:\n", iso_map((x1, y1)), "\n")
