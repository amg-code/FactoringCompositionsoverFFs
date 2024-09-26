# FactoringCompositionsoverFFs

This repository is licensed under the terms of the GNU General Public License v3.0 (GPL-3.0).

In this repository you can find an implementation of the factorization algorithm induced by my new closed formula for the factorization of every composition of the form f(X^n) for arbitrary monic irreducible polynomials f over a finite field F_q and for arbitrary positive integers n satisfying gcd(n,q)=1 in SageMath/Python.  

The implementation makes use of the two new SageMath/Python classes RichFiniteField and RichPolynomial which are useful for working with univariate polynomials over different finite fields and their extensions.

### Table of contents:
- [The new factorization algorithm](https://github.com/amg-code/FactoringCompositionsoverFFs#the-new-xn-factorization-algorithm)
- [How rich are the two RichClasses?](https://github.com/amg-code/FactoringCompositionsoverFFs#how-rich-are-the-two-richclasses)
- [Author](https://github.com/amg-code/FactoringCompositionsoverFFs#author)

## The new factorization algorithm




## How rich are the two RichClasses?
__RichFiniteField__ 

This class stores a SageMath FiniteField together with its PolynomialRing so that these two can be treated as a unit and used together. The finite field will always be of a primitive modulus so that primitive roots of unity in this finite field can easily be constructed by taking the generator to the respective exponent. 

The class has a kid called __RichExtensionField__. It stores a RichFiniteField which is an extension field of another RichFiniteField-instance - called the basefield. This class can cast elements and RichPolynomials from one of the two fields to the other. Furthermore, it can compute the minimal polynomial and the characteristic polynomial of elements of the extension field over the basefield. 

__RichPolynomial__ 

This class stores a polynomial over a given finite field (stored as a RichFiniteField) as a list and as a polynomial. This makes working with the polynomial as a list and as a polynomial in parallel easy.

It is enRICHed with many functions, some of which are redirections to existing SageMath-functions, many others are new implementations. 
All functions are split into private functions doing the computations and storing the result in a class attribute and a public function returning this attribute. This has the big advantage that all computations are done exactly once and only if needed. 



#### Author
Anna-Maurin Graner

