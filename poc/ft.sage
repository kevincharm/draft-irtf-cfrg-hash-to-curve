#!/usr/bin/sage
# vim: syntax=python

import sys
try:
    from sagelib.common import CMOV
    from sagelib.generic_map import GenericMap
except ImportError:
    sys.exit("Error loading preprocessed sage files. Try running `make clean pyfiles`")

class FouqueTibouchi(GenericMap):
    def __init__(self, F, A, B):
        self.name = "FT"
        self.F = F
        self.A = F(A)
        self.B = F(B)
        self.E = EllipticCurve(F, [F(A), F(B)])

    def sqrt_mod_p_checked(self, a):
        p = self.F.order()
        ret = pow(a, (p+1)/4, p)
        has_root = ret*ret == a
        return ret , has_root
    
    def sqrt_mod_p(self, a):
        p = self.F.order()
        return pow(a, (p+1)/4, p)

    def legendre(self, a):
        p = self.F.order()
        x = pow(a, (p-1)//2, p)
        if x == 0 or x == 1:
            return x
        if x == p-1:
            return -1
        assert False

    def map_to_point_ft(self, u):
        sqrt = self.sqrt_mod_p_checked
        inv0 = self.inv0
        F = self.F
        x = F(u)
        Z0 = F(0x0000000000000000b3c4d79d41a91759a9e4c7e359b6b89eaec68e62effffffd)
        Z1 = F(0x000000000000000059e26bcea0d48bacd4f263f1acdb5c4f5763473177fffffe)

        _, decision = sqrt(x)

        a0 = x * x
        a0 = a0 + 4
        a1 = x * Z0
        a2 = a1 * a0
        a2 = inv0(a2)
        a1 = a1 * a1
        a1 = a1 * a2

        # x1
        a1 = x * a1
        x = Z1 - a1
        # check curve
        a1 = x * x
        a1 = a1 * x
        a1 = a1 + 3
        a1, found = sqrt(a1)
        if found:
            if not decision:
                a1 = -a1
            return x, a1

        # x2
        x = -(x+1)
        # check curve
        a1 = x * x
        a1 = a1 * x
        a1 = a1 + 3
        a1, found = sqrt(a1)
        if found:
            if not decision:
                a1 = -a1
            return x, a1

        # x3
        x = a0 * a0
        x = x * x
        x = x * a2
        x = x * a2
        x = x + 1
        # must be on curve
        a1 = x * x
        a1 = a1 * x
        a1 = a1 + 3
        a1, found = sqrt(a1)
        if not found:
            raise Exception("MapToPointFailed")
        if not decision:
            a1 = -a1
        return x, a1

    # Another impl of FT to test against
    # From "Indifferentiable Hashing to Barreto-Naehrig Curves"
    # https://www.di.ens.fr/~fouque/pub/latincrypt12.pdf
    # yoinked from: https://github.com/randombit/pairings.py/blob/master/bn256.py#L1097
    def map_to_point_ft2(self, u):
        inv0 = self.inv0
        legendre = self.legendre
        sqrt = self.sqrt_mod_p
        t = self.F(u)
        A = self.A
        B = self.B
        p = self.F.order()

        def g(x):
            return x*x*x + B

        # constants
        sqrt_neg_3 = sqrt(p-3)
        inv_2 = inv0(2)
        b = B

        t2 = (t*t) % p

        chi_t = legendre(t)

        w = sqrt_neg_3 * t * inv0(1 + b + t2)

        x1 = ((sqrt_neg_3 - 1) * inv_2 - t*w) % p
        g_x1 = g(x1)
        if legendre(g_x1) == 1:
            x1_sqrt = sqrt(g_x1)
            return x1, (chi_t * x1_sqrt) % p

        x2 = (-1 - x1) % p
        g_x2 = g(x2)

        if legendre(g_x2) == 1:
            x2_sqrt = sqrt(g_x2)
            return x2, (chi_t * x2_sqrt) % p

        x3 = 1 + inv0(w*w)
        g_x3 = g(x3)

        assert legendre(g_x3) == 1
        x3_sqrt = sqrt(g_x3)
        return x3, (chi_t * x3_sqrt) % p

if __name__ == "__main__":
    u = 7105195380181880595384217009108718366423089053558315283835256316808390512725
    F = GF(21888242871839275222246405745257275088696311157297823662689037894645226208583)
    ft = FouqueTibouchi(F, 0, 3)
    x0, y0 = ft.map_to_point_ft(u)
    x1, y1 = ft.map_to_point_ft2(u)
    assert x0 == x1 and y0 == y1, "FT implementations don't match"
    print(x0, y0)
    assert x0 == 19485131671658523517646027848906165907640971588430452127920614621547697012573, "x doesn't match ref"
    assert y0 == 7252485661626054658053752721536032361940074412825453078837989033903251969412, "y doesn't match ref"
