# *************************************************************************
#       Copyright (C) 2023 Anna-Maurin Graner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *************************************************************************

import itertools

def pause():
    input("\n(Press ENTER to continue and CTR+C to stop)")
    return

#___________________________________________________
# Functions in number theory:

# function splits the positive integer n into two positive integers n' and m such that gcd(n',q)=1 and rad(m)|q
# where q=p^l is a prime power
def make_gcd_one(q,n):
    p=((q.factor())[0])[0]
    l=n.valuation(p)
    return (int(n/p^l), p^l)

# function checks whether a given integer is a prime power
def is_prime_power(q):
    if q==0:
        return False
    else:
        if len(factor(q))==1:
            return True
        else: 
            return False
    return 

# function returns the radical of a positive integer n:
def rad(n):
    return  prod([fac[0] for fac in factor(n)])

# function returns all k-subsets of the set s:
def k_subsets(s,k):
    return set(itertools.combinations(s,k))

# This function returns the cyclotomic coset of i mod d over Fq
def cycl_coset(q,d,i):
    return [i*q^j % d for j in range((Mod(q,d/(gcd(d,i))).multiplicative_order()))]

# This function returns a representative system for the q-cyclotomic classes modulo d
def cycl_coset_repsystem(q,d):
    if d==1:
        coset_reps = {0}
    elif d>1:            
        coset_reps = set()
        tbd = list(range(0,d))
        while len(tbd)>0:
            i = tbd[0]
            coset_reps.add(i)
            tbd.remove(i)
            cyclel = i*q % d 

            while not cyclel == i:
                tbd.remove(cyclel)
                cyclel = cyclel*q % d 
    return coset_reps

#___________________________________________________
# Functions in finite fields:

# inverse function on a finite field
def inverse(a,F):
    return a^(F.cardinality()-2)

#___________________________________________________
# Functions for polynomials:

# f(X^n) for f given as a list:
def composition_Xn(f,x, n):
    k=len(f)
    return sum([f[i]*x^(n*(k-1-i)) for i in range(k)])

# function computes the q-spin of a given polynomial pol in the variable x
# coefficient degree is already given as cd 
def spin(q,x,pol,cd):
    poll = pol.list()
    spin = pol
    for j in range(1,cd):
        spin = spin * sum([poll[i]^(q^j)*x^i for i in range(len(poll))])
    return spin 

# function computes the q-spin of a binomial X^t-b given as (x,t,b)
# coefficient degree is already given as cd 
def spin_binomial(q,x,t,b,cd):
    spin =0
    for i in range(cd+1):
        bproduct = sum([prod([b^(q^j) for j in J]) for J in k_subsets(range(cd),cd-i)])
        spin = spin+x^(t*i)*(-1)^(cd-i)*bproduct
    return spin

# function casts polynomial in X^d to the basefield of REF
def pol_inXd_to_basefield(pol, REF, x, d):
    pol = pol.list()
    pol = [pol[l] for l in range(len(pol)) if l%d ==0]
    pol = sum(REF.to_basef(pol[l])*x^(d*l) for l in range(len(pol)))
    return pol  


# function computes the q-spin of the linear factor X-b given as x,b
# with the explicit formula based on k-subsets
# coefficient degree is already given as cd
def spin_linear_factor_explicit(q,x,b,cd):
    spin =0
    for i in range(cd+1):
        bproduct = sum([prod([b^(q^j) for j in J]) for J in k_subsets(range(cd),cd-i)])
        spin = spin+x^(i)*(-1)^(cd-i)*bproduct
    return spin

# functions prints the factorization as Magma does
def print_magma_style(factrztn):
    print(str_magma_style(factrztn))
    return

# function returns a string of a Factorization object in Magma-style
def str_magma_style(factrztn):
    return "\t"+"\n\t".join(["<"+str(fact[0])+","+str(fact[1])+">" for fact in factrztn])

# function returns a string of a polynomial for use with LaTeX:
def latex_pol(f,x):
    f=f.list()
    s=""
    for i in range(len(f)):
        if not f[i]==0:
            if not f[i]==1:
                sn = str(f[i])
            else:
                sn =""

            if not i in [0,1]:
                sn=sn+x+"^{"+str(i)+"}"
            if i==0:
                if f[i]==1:
                    sn = str(f[i])
            if i ==1:
                sn = sn+x
                
            if len(s)>0:
                s=sn+"+"+s
            else:
                s=sn
    return s
