from RichFiniteFieldClass import *
from RichPolynomialClass import *

q = 4

# Create a RichFiniteField instance of size q, with generator displayed as "a" - F_q = F_p(a)- and a polynomial ring in the variable displayed as "X":
RFF = RichFiniteField(q, "a", "X")
a = RFF.gen

# Get all information on this RichFiniteField
print(RFF)

# The SageMath finite field GF(q) hidden inside of the RichFiniteField can be accessed with:
print(RFF.F)

# The RichFiniteField already contains a Polynomial Ring in the variable "X"
# Its generator can be accessed with
x = RFF.x

# Create a RichPolynomial f over RFF:
# f = X^2+a
f = RichPolynomial([1,1,"a"], RFF)

# Get all information on f:
print(f.get_full_info())

# or every property separately:
print(f.get_weight())
print(f.get_order())
print(f.get_knormal())
print(f.get_primitive())
print(f.is_irreducible())

# Factorization:
print(f.get_factorization())
print(f.get_factors_asRP())

# Get roots of f over F_q:
print(f.roots())

# Create another RichPolynomial g over RFF:
g = RichPolynomial(x^3+a^2*x, RFF)
print(g.get_full_info())
print(g.roots())

# Basic arithmetic operations are implemented for RichPolynomials:
print(f+g)
print(f*g)
print(f==g)

# Lets create an extension of RFF:
k = 2
REF = RichExtensionField(RFF, k, "b", "Y")
print(REF)

# Cast f to the extension field and consider it as a polynomial over F_{q^k}
f_inREF = REF.pol_to_extf(f)
print(f_inREF)

# Get roots over F_{16}:
print(f_inREF.roots())

# It is not necessary to build the splitting field of f in general, it is already implemented in f:
print(f.get_all_roots())

# If you wish to access the splitting field:
print(f.get_splitting_field())






