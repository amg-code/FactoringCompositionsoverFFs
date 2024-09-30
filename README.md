# FactoringCompositionsoverFFs

This repository is licensed under the terms of the GNU General Public License v3.0 (GPL-3.0).

I present a SageMath/Python implementation of the factorization algorithm induced by the new closed formula for 
the factorization of every composition of the form f(X^n) over a finite field F_q for arbitrary monic irreducible polynomials f and for arbitrary positive integers n satisfying gcd(n,q)=1 as presented in the ArXiV preprint 

__"Closed formulas for the factorization of X^n−1, the n-th cyclotomic polynomial, X^n−a and f(X^n) over a finite field for arbitrary positive integers n" by Anna-Maurin Graner__ [https://arxiv.org/abs/2306.11183]. 

The implementation makes use of the two new SageMath/Python classes  __RichFiniteField__ and __RichPolynomial__  which are useful for working with univariate polynomials over different finite fields and their extensions. You can find a short tutorial on how to use the two classes in the file __RichTutorial.sage__. 

### Table of contents:
- [How to preparse .sage files](https://github.com/amg-code/FactoringCompositionsoverFFs#how-to-preparse-sage-files)
- [How to use the new factorization algorithm](https://github.com/amg-code/FactoringCompositionsoverFFs#how-to-use-the-new-factorization-algorithm)
- [Measuring the new factorization algorithm's computation time](https://github.com/amg-code/FactoringCompositionsoverFFs#measuring-the-new-factorization-algorithms-computation-time)
- [How rich are the two RichClasses?](https://github.com/amg-code/FactoringCompositionsoverFFs#how-rich-are-the-two-richclasses)

## How to preparse .sage files

Several .sage files in this repository have to be preparsed so that they can be loaded as Python modules into other .sage files. We include the preparsed versions of the files in this repository. However, if you wish to preparse the files yourself, please navigate to the folder that the file (i.e. fXnAlgorithm.sage) is in and use the following command:

`sage --preparse fXnAlgorithm.sage`

The resulting file __fXnAlgorithm.sage.py__ will appear in the current working directory and can be renamed to __fXnAlgorithm.py__. 



## How to use the new factorization algorithm

The new factorization algorithm is given in the file __fXnAlggorithm.sage__. This file can be preparsed and loaded as a Python module. The algorithm itself is called with the function

`factor_fXn(f, n, printing)`

where __f__ is a RichPolynomial instance, __n__ a positive integer and __printing__ either `True` or `False`, depending on whether you wish to print all steps of the algorithm or not. Note that printing makes the program much slower. 

The function can be used as follows:

    from fXnAlgorithm import * 
    from UsefulFunctions import *      # for Magma-style printing
    
    q=2
    
    RFF=RichFiniteField(q,False,"X")
    
    f=RichPolynomial([1,1,1,0,1,0,1],RFF)
    
    n=2^3*7*3^3*11^2
    
    print(factor_fXn(f,n,1)     # with printing all algorithm steps
    
    print(factor_fXn(f,n,0)     # without printing all algorithm steps 

    print_Magma_style(factor_fXn(f,n,0))    # print in the same style as Magma

## Measuring the new factorization algorithm's computation time

If you wish to measure the execution time of the new algorithm and/or compare it with the existing SageMath `factor()` function, you can use the file __TestNewAlgorithm.sage__.  Look for `USER INPUT` and specify the examples you wish to compute as tuples `(q,n,f)` (look for USER INPUT). 

Note that the program supports parallel computations for multiple polynomials at the same time. For this the Python module __multiprocessing__ is used. However, the PARI implementation of SageMath's `factor()` function is not compatible with the Python module __signal__ which is used for the computation timeout. We recommend that you set the variable `parallel_computations` to `False` if you run into problems or out of RAM. 

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

