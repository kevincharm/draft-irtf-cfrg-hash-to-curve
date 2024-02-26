#!/usr/bin/env sage

"""
    From Wahby & Boneh: https://eprint.iacr.org/2019/403.pdf
"""

import sage.schemes.elliptic_curves.isogeny_small_degree as isd
# look for isogenous curves having j-invariant not in {0, 1728}
# Caution: this can take a while!
def find_iso(E):
    for p_test in primes(100):
        print(p_test)
        isos = [ i for i in isd.isogenies_prime_degree(E, p_test)
            if i.codomain().j_invariant() not in (0, 1728) ]
        if len(isos) > 0:
            return isos[0].dual()
    return None

if __name__ == "__main__":
    # BN254 parameters
    p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    assert is_prime(p)
    # E
    F = GF(p)
    E = EllipticCurve(F, [0, 3]) # y^2 = x^3 + 0*x + 3
    # Find isogenous curve E'
    iso_G1 = find_iso(E) # an isogeny from E’ to E,
    E_prime = iso_G1.domain() # where this is E’
    assert iso_G1(E_prime.random_point()).curve() == E
    print("E'", E_prime)
