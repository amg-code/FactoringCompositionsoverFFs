# *************************************************************************
#       Copyright (C) 2023 Anna-Maurin Graner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *************************************************************************

# This file contains the following Python Classes for computations with polynomials over finite fields and their extension fields 
# - RichPolynomial

# The class RichPolynomial is useful for working with a polynomial over a given RichFiniteField-instance as a list and as a polynomial in parallel. 
# Furthermore, it is enRICHed with many useful functions returning properties of the polynomial such as its weight, its order, its k-normality, its primitivity. 

# This class contains many functions some of which are redirections to existing SageMath-functions, many others are new implementations.
# All functions are split into private functions doing the computations and storing the result in a class attribute and a public function returning this attribute. This has the big advantage that all computations are done exactly once and only if needed.

from sage.arith.functions import LCM_list
import RichFiniteFieldClass

sep = "\t"

class RichPolynomial:
    def __init__(self,pol,RFF):
        self._RFF = RFF

        if isinstance(pol, list):
            self.list = [RFF.F(pol[i]) for i in range(len(pol))]
            self.deg = len(self.list)-1
        else:
            self._polynomial = pol
            self.deg = pol.degree()
            self.list = pol.list()
            (self.list).reverse()
        return

    def __str__(self):
        return str(self.list)+" = "+str(self.get_polynomial())

    def __mul__(self,pol):
        return RichPolynomial(self.get_polynomial()*pol.get_polynomial(), self._RFF)

    def __add__(self,pol):
        return RichPolynomial(self.get_polynomial()+pol.get_polynomial(), self._RFF)

    def __eq__(self,pol):
        return self.list == pol.list 

    def __neq__(self,pol):
        return not self.list == pol.list

    # If self is a binomial with nonzero constant term, this functions returns a string representing this polynomial of the form X^t-a
    def str_binomial(self):
        if self.get_weight() >2 or (self.list)[-1] == (self._RFF).zero:
            raise ValueError("The polynomial "+str(self)+" is no binomial!")
        return str((self._RFF).charvar)+"^"+str(self.deg)+" - ("+str(-(self.list)[-1])+")"

    # This function returns a string detailing all information that can be obtained from a RichPolynomial-instance
    def get_full_info(self):
        s = str(self.list)+" = "+str(self.get_polynomial())+"\n"+sep+"irreducible: "+sep+str(self.is_irreducible())+"\n"+sep+"weight: "+sep+str(self.get_weight())
        if self.is_irreducible():
            s= s+"\n"+sep+"order: "+sep+str(self.get_order())+" = "+str((self.get_order()).factor())+"\n"+sep+"primitive: "+sep+str(self.get_primitive())+"\n"+sep+"k-normality: "+sep+str(self.get_knormal())
        return s

    # This function returns the finite field over which self is defined
    def get_RFF(self):
        return self._RFF

    #public: returns the polynomial in F[X]
    def get_polynomial(self):
        if hasattr(self,"_polynomial"):
            return self._polynomial
        else:
            self._add_polynomial()
            return self._polynomial

    #private: adds polynomial in F[X] to attributes
    def _add_polynomial(self):
        self._polynomial = sum([self.list[i]*(self._RFF).x^(self.deg-i) for i in range(self.deg+1)])
        return

    # public: get splitting field of polynomial as RichExtensionField-instance:
    def get_splitting_field(self):
        if not hasattr(self, "_RSplF"):
            self._add_splitting_field()
        return self._RSplF

    # private: if polynomial is irreducible, this function adds the splitting field as a RichExtensionField-instance
    def _add_splitting_field(self):
        if not hasattr(self, "_RSplF"):
            if self.is_irreducible():
                self._RSplF = RichFiniteFieldClass.RichExtensionField(self._RFF, self.deg, "z", "X")
            else: 
                max_deg = LCM_list([fac[0].deg for fac in self.get_factors_asRP()])
                self._RSplF = RichFiniteFieldClass.RichExtensionField(self._RFF,max_deg, "z", "X")            
        return

    #public: returns self as RichPolynomial-instance over the splitting field
    def get_RP_in_splitting_field(self):
        if not hasattr(self, "_in_RSplF"):
            self._add_pol_in_RSplF()
        return self._in_RSplF

    #private: adds self as RichPolynomial-instance over the splitting field
    def _add_pol_in_RSplF(self):
        self._in_RSplF = RichPolynomial([(self.get_splitting_field()).to_extf(e) for e in self.list], self.get_splitting_field())
        return

    #public: returns whether polynomial is irreducible
    def is_irreducible(self):
        if not hasattr(self,"_irreducible"):
            self._add_irreducible()
        return self._irreducible

    #private: determines whether polynomial is irreducible
    def _add_irreducible(self):
        self._irreducible = (self.get_polynomial()).is_irreducible()
        return
    
    # public: returns the factorization of the polynomial over the finite field (as SageMath Factorization object)
    def get_factorization(self):
        if not hasattr(self, "_factorization"):
            self._add_factorization()
        return self._factorization
    
    # private: adds the factorization to the attributes of the polynomial
    def _add_factorization(self):
        self._factorization = (self.get_polynomial()).factor()
        return

    #public: returns the list of factors as RichPolynomial-instances
    def get_factors_asRP(self):
        if not hasattr(self, "_factors"):
            self._add_factors()
        return self._factors
    
    # private: adds a list of factors as RichPolynomial-instances
    def _add_factors(self):
        self._factors = [(RichPolynomial( g[0], self._RFF), g[1]) for g in self.get_factorization()]
        return

    # Returns one root of the polynomial in _RFF
    def root(self):
        if len(self.roots())>0:
            return ((self.roots())[0])[0]
        else: 
            raise ValueError("The polynomial "+str(self)+" does not have a root in the finite field of size "+str((self._RFF).q)+".")
        return 

    # Returns all roots of the polynomial in _RFF
    def roots(self):
        if not hasattr(self, "_roots"):
            self._add_roots()
        return self._roots

    # private: 
    def _add_roots(self):
        self._roots = (self.get_polynomial()).roots()
        return

    # public: Returns a list of the roots of all irreducible factors of the polynomial over the respective extension fields
    def get_all_roots(self):
        if not hasattr(self, "_all_roots"):
            self._add_all_roots()
        return self._all_roots  

    # private: computes the roots of all irreducible factors over the respective extension field
    # Returns a list (sorted in the same order as the list self._factors) of the form 
    # [(factor, q^k, multiplicity, [root1, root2, root3]), ... ] 
    def _add_all_roots(self):    
        self._all_roots = [(fac[0].list, (fac[0].get_splitting_field()).q, fac[1], (fac[0].get_RP_in_splitting_field()).roots()) for fac in self.get_factors_asRP()]
        return     
      
    #private: adds the variable order to polynomial  
    def _add_order(self):
        if not self.is_irreducible():
            print("The polynomial "+str(self)+" is not irreducible. Therefore has no order.")
            self._order = 0
        else:
            self._order = ((self.get_RP_in_splitting_field()).root()).multiplicative_order()
        return 

    #public: Returns the order of the polynomial
    def get_order(self):
        if not hasattr(self,"_order"):
            self._add_order()
        return self._order

    #private: Computes and sets the variable primitive to True or False
    def _add_primitive(self):        
        self._primitive = (self.get_order()==((self._RFF).q)^self.deg-1)
        return

    #public: Return True/False whether polynomial is primitive
    def get_primitive(self):
        if not hasattr(self, "_primitive"):
            self._add_primitive()
        return self._primitive

    # public: Returns the weight of the polynomial
    def get_weight(self):
        if not hasattr(self,"_weight"):
            self._add_weight()
        return self._weight

    #private: compute weight of polynomial -> package multiset needed!
    def _add_weight(self):
        #Version without multiset dependence since this package is missing from the SageMath version of Python
        self._weight = sum([1 for a in self.list if not a==0])
        return 

    # public: Returns the positive integer k such that the polynomial is k-normal
    def get_knormal(self):
        if hasattr(self,"_knormal"):
            return self._knormal
        else:
            self._add_knormal()
            return self._knormal
    
    #private: computes the k-normality of the polynomial and stores it in self._knormal
    def _add_knormal(self):
        alpha = (self.get_RP_in_splitting_field()).root()
        f=((self.get_splitting_field()).x)^(self.deg)-1
        g=sum([alpha^((self._RFF).q^i)*((self.get_splitting_field()).x)^(self.deg-1-i) for i in range(self.deg)])
        self._knormal = (f.gcd(g)).degree()
        return

    # Checks if the polynomial is a composition of the form g(X^m)
    def is_composition(self, m):
        for i in [j for j in range(len(self.list)) if j%m !=0]:
            if self.list[i] !=0:
                return False
        return True

    # Extracts the polynomial g from the composition f = g(X^m)
    def extract(self, m):
        if self.is_composition(m):
            return RichPolynomial([self.list[i] for i in range(len(self.list)) if i%m == 0], self._RFF) 
        else:
            raise ValueError("The polynomial "+str(self)+" is no composition of the form g(X^"+str(m)+").")
        return

    # This function constructs m_{beta^p} from self, where beta is a root of self and where p is a prime dividing (q-1)
    # It is an implementation of Corollary 8 from the paper <<Constructing irreducible polynomials recursively with a reverse composition method>> by Anna-Maurin Graner and Gohar M. Kyureghyan 
    def construction_mbetap(self, p):
        # checks if p divides q-1:
        if p not in (self._RFF).get_primefactors():
            raise ValueError(str(p)+" is not a prime factor of "+str((self._RFF).q)+"-1.")

        if self.is_composition(p):
            # if f = g(X^p), then m_{beta^p}(X^p)=f
            m_betap_Xp = self 
        else:
            # if f is no composition with X^p, then m_{beta^p}(X^p) is computed with formula:
            m_beta = self.get_polynomial()
            zeta_p = (self._RFF).get_unityroot(p)

            m_betap_Xp = (-1)^(self.deg*(p+1)) * prod([m_beta(zeta_p^j*(self._RFF).x) for j in range(p)]) 
            m_betap_Xp = RichPolynomial(m_betap_Xp, self._RFF)

        return m_betap_Xp.extract(p)

    
