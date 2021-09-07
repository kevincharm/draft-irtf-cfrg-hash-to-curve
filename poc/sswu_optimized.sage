#!/usr/bin/sage
# vim: syntax=python

import sys
try:
    from sagelib.common import CMOV
    from sagelib.generic_map import GenericMap
    from sagelib.z_selection import find_z_sswu
    from sagelib.sqrt import sqrt_checked, sqrt_ratio_straightline
except ImportError:
    sys.exit("Error loading preprocessed sage files. Try running `make clean pyfiles`")

class OptimizedSSWU(GenericMap):
    def __init__(self, F, A, B):
        self.name = "SSWU"
        self.F = F
        self.A = F(A)
        self.B = F(B)

        if self.A == 0:
            raise ValueError("S-SWU requires A != 0")
        if self.B == 0:
            raise ValueError("S-SWU requires B != 0")
        self.Z = find_z_sswu(F, F(A), F(B))
        self.E = EllipticCurve(F, [F(A), F(B)])

        # values at which the map is undefined
        # i.e., when Z^2 * u^4 + Z * u^2 = 0
        # which is at u = 0 and when Z * u^2 = -1
        c = -F(1) / self.Z
        self.undefs = [F(0)]
        if c.is_square():
            ex = c.sqrt()
            self.undefs += [ex, -ex]

    def not_straight_line(self, u):
        inv0 = self.inv0
        is_square = self.is_square
        sgn0 = self.sgn0
        sqrt = self.sqrt
        u = self.F(u)
        A = self.A
        B = self.B
        Z = self.Z

        tv1 = inv0(Z^2 * u^4 + Z * u^2)
        x1 = (-B / A) * (1 + tv1)
        if tv1 == 0:
            x1 = B / (Z * A)
        gx1 = x1^3 + A * x1 + B
        x2 = Z * u^2 * x1
        gx2 = x2^3 + A * x2 + B
        if is_square(gx1):
            x = x1
            y = sqrt(gx1)
        else:
            x = x2
            y = sqrt(gx2)
        if sgn0(u) != sgn0(y):
            y = -y
        return (x, y)

    def sqrt_ratio(self, u, v):
        x = self.F(u) / self.F(v)
        r1 = sqrt_checked(self.F, x, self.Z)
        r2 = sqrt_ratio_straightline(self.F, u, v, self.Z)
        assert r1 == r2
        return r2

    def straight_line(self, u):
        A = self.A
        B = self.B
        Z = self.Z
        u = self.F(u)
        sqrt_ratio = self.sqrt_ratio
        sgn0 = self.sgn0

        tv1 = u^2
        tv1 = Z * tv1
        tv2 = tv1^2
        tv2 = tv2 + tv1
        tv3 = tv2 + 1
        tv3 = B * tv3
        tv4 = CMOV(Z, -tv2, tv2 != 0)
        tv4 = A * tv4
        tv2 = tv3^2
        tv6 = tv4^2
        tv5 = A * tv6
        tv2 = tv2 + tv5
        tv2 = tv2 * tv3
        tv6 = tv6 * tv4
        tv5 = B * tv6
        tv2 = tv2 + tv5
        x = tv1 * tv3
        (is_gx1_square, y1) = sqrt_ratio(tv2, tv6)
        y = tv1 * u
        y = y * y1
        x = CMOV(x, tv3, is_gx1_square)
        y = CMOV(y, y1, is_gx1_square)
        e1 = sgn0(u) == sgn0(y)
        y = CMOV(-y, y, e1)
        x = x / tv4

        return (x, y)

if __name__ == "__main__":
    for _ in range(0, 32):
        OptimizedSSWU.test_random()
