from khovanov import KhovanovComplex
from config import CONVENTION

class Jones:

    def __init__(self, kh):
        """
        kh is a KhovanovComplex object.
        Assumes:
            chi_q(Kh) = (q + q^-1) V(q^2)
        """
        self.kh = kh

    def khovanov_euler_poly(self):
        """
        Returns:
            dict exponent -> coefficient

        Represents:
            sum_{i,j} (-1)^i rank Kh^{i,j} q^j
        """

        poly = {}

        for (i, j) in self.kh.chain_basis:
            rank = self.kh.homology_rank(i, j)

            if rank == 0:
                continue

            coeff = ((-1) ** i) * rank
            poly[j] = poly.get(j, 0) + coeff

            if poly[j] == 0:
                del poly[j]

        return poly

    def laurent_divide(self, numerator, denominator):
        """
        Laurent polynomial division.

        Polynomials are dicts:
            exponent -> coefficient

        Assumes exact division.
        """

        num = dict(numerator)
        den = dict(denominator)

        quotient = {}

        while num:
            lead_num_exp = max(num.keys())
            lead_den_exp = max(den.keys())

            lead_num_coeff = num[lead_num_exp]
            lead_den_coeff = den[lead_den_exp]

            exp_shift = lead_num_exp - lead_den_exp

            if lead_num_coeff % lead_den_coeff != 0:
                raise ValueError("Non-integer coefficient during Laurent division.")

            coeff = lead_num_coeff // lead_den_coeff

            quotient[exp_shift] = quotient.get(exp_shift, 0) + coeff
            if quotient[exp_shift] == 0:
                del quotient[exp_shift]

            for e, c in den.items():
                new_e = e + exp_shift
                num[new_e] = num.get(new_e, 0) - coeff * c

                if num[new_e] == 0:
                    del num[new_e]

        return quotient

    def jones_q2_poly(self):
        """
        Returns V(q^2), using:
            chi_q(Kh) = (q + q^-1) V(q^2)
        """

        euler = self.khovanov_euler_poly()
        divisor = {1: 1, -1: 1}

        return self.laurent_divide(euler, divisor)

    def jones_poly(self):
        """
        Returns V(t), where t = q^2.

        Output:
            dict exponent -> coefficient
        """

        q2_poly = self.jones_q2_poly()
        t_poly = {}

        for q_exp, coeff in q2_poly.items():
            if q_exp % 2 != 0:
                raise ValueError(f"Expected even q-exponent, got q^{q_exp}")

            t_poly[q_exp // 2] = coeff

        return t_poly

    def format_laurent_poly(self, poly, var="q"):
        """
        Pretty-prints dict exponent -> coefficient.
        """

        if not poly:
            return "0"

        pieces = []

        for exp in sorted(poly.keys(), reverse=True):
            coeff = poly[exp]

            if coeff == 0:
                continue

            abs_coeff = abs(coeff)

            if exp == 0:
                mono = str(abs_coeff)
            else:
                if abs_coeff == 1:
                    coeff_part = ""
                else:
                    coeff_part = str(abs_coeff)

                if exp == 1:
                    mono = coeff_part + var
                else:
                    mono = coeff_part + f"{var}^{exp}"

            if not pieces:
                if coeff < 0:
                    pieces.append("-" + mono)
                else:
                    pieces.append(mono)
            else:
                if coeff < 0:
                    pieces.append(" - " + mono)
                else:
                    pieces.append(" + " + mono)

        return "".join(pieces)
    
    def invert_poly(self, poly):
        """
        Given P(q), returns P(q^-1).

        If:
            P(q) = q^2 - q + 1

        then:
            P(q^-1) = q^-2 - q^-1 + 1
        """

        inverted = {}

        for exp, coeff in poly.items():
            inverted[-exp] = inverted.get(-exp, 0) + coeff

        return inverted

    def print_euler(self):
        poly = self.khovanov_euler_poly()
        print(self.format_laurent_poly(poly, "q"))

    def print_jones_q2(self):
        poly = self.jones_q2_poly()
        print(self.format_laurent_poly(poly, "q"))

    def print_jones_original(self):
        poly = self.jones_poly()
        print(self.format_laurent_poly(poly, "t"))

    def jones_poly_inverse(self):
        """
        Returns V(t^-1), where V is the Jones polynomial in variable t.
        """

        return self.invert_poly(self.jones_poly())

    def print_jones_inverse(self):
        poly = self.jones_poly_inverse()
        print(self.format_laurent_poly(poly, "q"))

    def print_jones(self):
        if CONVENTION == "regina":
            self.print_jones_original()
        elif CONVENTION == "katlas":
            self.print_jones_inverse()

# Usage
# K = KhovanovComplex("dabcabcv-")
# J = Jones(K)
# J.print_jones()