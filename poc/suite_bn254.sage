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
    from sagelib.sswu_opt_3mod4 import OptimizedSSWU_3mod4
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

def generate_svdw_test_vecs(svdw):
    print("[")
    for i in range(0, 1000):
        print("{")
        u = svdw.F.random_element()
        print(f"\"u\": \"{u}\",")
        x, y = svdw.not_straight_line(u)
        print(f"\"p\": [\"{x}\",\"{y}\"]")
        print("},")
    print("]")

if __name__ == "__main__":
    u = 3540903031681319421922757684101610645767707797048988415875375111724680581685
    svdw = GenericSvdW(F, A, B)
    x1, y1 = svdw.not_straight_line(u)
    generate_svdw_test_vecs(svdw)
