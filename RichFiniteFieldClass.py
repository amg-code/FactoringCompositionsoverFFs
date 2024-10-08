

# This file was *autogenerated* from the file RichFiniteFieldClass.sage
from sage.all_cmdline import *   # import sage library

_sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_10 = Integer(10)# ****************************************************************************
#       Copyright (C) 2023 Anna-Maurin Graner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# This file contains the following Python Classes for computations with polynomials over finite fields and their extension fields 
# - RichFiniteField
# - RichExtensionField (kid of RichFiniteField)

# The class stores a finite field together with its polynomial ring so that these two structures can be used together easily


from sage.coding.relative_finite_field_extension import *
import RichPolynomialClass

sep = "\t"

# CLASSES

# The Class RichFiniteField stores a Sagemath finite field instance _FF and a polynomial ring over this finite field 

class RichFiniteField:
    def __init__(self, qq, charroot, charvar):
        self.q = qq
        self.F = GF(self.q, charroot)#, modulus="primitive")
        self.gen = (self.F).gen() 
        self.zero = (self.F)(_sage_const_0 )
        self.one = (self.F)(_sage_const_1 )
        self._PRF = PolynomialRing(self.F,charvar)
        self.charvar = charvar
        self.x = (self._PRF).gen()
        return        
    
    def __str__(self):
        return "RichFiniteField of size q = "+str(self.q)+" = "+str(factor(self.q))+"\n"+sep+"q-1 = "+str(self.q-_sage_const_1 )+" = "+str(factor(self.q-_sage_const_1 ))+"\n"+sep+"modulus: "+sep+str((self.F).modulus())+"\n"+sep+"generator: "+sep+str(self.gen)+"\n"+sep+"polynomial ring in: "+sep+str(self.charvar)

    # public (redirection): returns the modulus of the finite field
    def get_modulus(self):
        return (self.F).modulus()

    # public: returns a list of the prime factors of q-1
    def get_primefactors(self):
        if not hasattr(self, "_primefactors"):
            self._add_primefactors()
        return self._primefactors

    # private: adds a list of the prime factors of q-1
    def _add_primefactors(self):
        self._primefactors = [fact[_sage_const_0 ] for fact in factor(self.q-_sage_const_1 )]
        #self._primefactors.reverse()
        return 
    
    def get_primitive_element(self):
        if not hasattr(self, "_primitive_element"):
            self._add_primitive_element()
        return self._primitive_element
    
    def _add_primitive_element(self):
        self._primitive_element = (self.F).multiplicative_generator()

    # public: returns a primitive n-th root of unity if n divides q-1
    def get_unityroot(self, n):
        if (self.q-_sage_const_1 )%n == _sage_const_0 :
            return (self.get_primitive_element())**(int((self.q-_sage_const_1 )/n))
        else:
            raise ValueError(str(n)+" does not divide "+str(self.q)+"-1. Thus, there is no primitive "+str(n)+"-th root of unity in F_{"+str(self.q)+"}.")
        return 

    # public: returns a dictionary of one n-th primitive root of unity for all divisors n of q-1
    def get_unityroots(self):
        self._add_unityroots()
        return self._unityroots

    # private: add the dictionary of n-th primitive roots of unity for all divisors n of q-1
    def _add_unityroots(self):
        if not hasattr(self, "_unityroots"):
            self._unityroots = dict()

            for d in (self.q-_sage_const_1 ).divisors():
                self._unityroots[d] = self.gen**(int((self.q-_sage_const_1 )/d))
        return 


# The class RichExtensionField which inherits from the class RichFiniteField is mainly an extension of another RichFiniteField-instance such that the elements of the base field are embedded in the extension field.
class RichExtensionField(RichFiniteField):
    def __init__(self, RFF, nn, charroot, charvar):
        super().__init__(RFF.q**nn, charroot, charvar)
        self._basefield = RFF 
        self._RFFE = RelativeFiniteFieldExtension(self.F,RFF.F)
        self._base_embedding = (self._RFFE).embedding()
        self._extension_degree = nn
        return

    def __str__(self):
        return "RichExtension"+(super().__str__())[_sage_const_10 :]+"\n"+sep+"extension of degree : "+sep+str(self._extension_degree)+"\n"+sep+"over : "+sep+((str(self._basefield)).split("\n"))[_sage_const_0 ]
    
    #public: function returns basefield as RichFiniteField-instance
    def get_basefield(self):
        return self._basefield

    #public: function returns the extension degree over basefield
    def extension_degree(self):
        return self._extension_degree

    # public: function embeds an element el from the base field in extension field self
    def to_extf(self, el):
        return self._base_embedding(el)

    def to_basef(self, el):
        return (self._RFFE).cast_into_relative_field(el)

    def in_basis_over_basef(self,el):
        return (self._RFFE).relative_field_representation(el)

    # Converts a RichPolynomial-instance over the basefield to a RichPolynomial-instance over the extension field
    def pol_to_extf(self,pol):
        return RichPolynomialClass.RichPolynomial([self.to_extf(el) for el in pol.list], self)

    # Converts a RichPolynomial-instance over the extension field to a RichPolynomial instance over the basefield
    def pol_to_basef(self, pol):
        return RichPolynomialClass.RichPolynomial([self.to_basef(el) for el in pol.list], self._basefield)


    # public: compute the minimal polynomial over the base field of an element of extension field
    def minimal_polynomial_over_basefield(self, el):
        mp  = self.x - el 
        for i in range(_sage_const_1 , self._extension_degree):
            eltoqi = el**(((self._basefield).q)**i)
            if eltoqi == el:
                break
            else:
                mp = mp * (self.x-eltoqi)
        mp = [self.to_basef(e) for e in mp.list()]
        mp.reverse()
        return RichPolynomialClass.RichPolynomial(mp,self._basefield) 

    # public: compute the characteristic polynomial over the base field of an element of the extension field
    def characteristic_polynomial_over_basefield(self, el):
        mp = self.minimal_polynomial_over_basefield(el)
        return RichPolynomialClass.RichPolynomial((mp.get_polynomial())**(int(self._extension_degree/mp.deg)), self._basefield)
        

