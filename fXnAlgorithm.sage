# *************************************************************************
#       Copyright (C) 2023 Anna-Maurin Graner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *************************************************************************

# This module contains the function factor_fXn which computes the factorization of f(X^n) for any positive integer n and any irreducible polynomial f over Fq
# with the new algorithm by Anna-Maurin Graner

from RichPolynomialClass import *
from RichFiniteFieldClass import *
from UsefulFunctions import *

# function computes the q-spin of a the linear factor x-b
# coefficient degree is already given as cd 
def spin_linear_factor(q,x,b,cd):
    spin = x-b
    for u in range(1,cd):
        spin = spin * (x-b^(q^u))
    return spin 

# Cast list in X^t to basefield of tower of two extension fields:
# mode:
#   3 - q_ks-> q_k -> q
#   2 - q_ks -> q_k=q
#   1 - q_k -> q
#   0 - none 
def cast_from_qks_to_q(mode, pol, RF_qks, RF_qk, x, d): 
    if mode==0:
        return pol(x^d)
    else:
        pol = pol.list()
        if mode == 3:
            return sum(RF_qk.to_basef(RF_qks.to_basef(pol[l]))*x^(l*d) for l in range(len(pol)))
        elif mode == 2:
            return sum(RF_qks.to_basef(pol[l])*x^(l*d) for l in range(len(pol)))
        elif mode == 1:
            return sum(RF_qk.to_basef(pol[l])*x^(l*d) for l in range(len(pol)))
    return 

# Factorization of f(X^n) with the new algorithm by Anna-Maurin Graner
# Input: RichPolynomial f, positive integer n
# Output: Factorization of f(X^n)

def factor_fXn(f,n, printing):
    if not f.is_irreducible():
        raise ValueError("The given polynomial f is not irreducible.")
    else:
        RF_q = f.get_RFF()
        q=RF_q.q
        k = f.deg
                    
        # Determine root of f and F_{q^k}
        if k==1:
            RF_qk = RF_q
            alpha = f.root()
        else:
            RF_qk = f.get_splitting_field()
            alpha = (f.get_RP_in_splitting_field()).root()
        
        q_k = RF_qk.q       

        # Determine nu_p(n) and tilde-alpha
        p_to_l = 1
        # Case gcd(q,n)>1: Using remark of Chapter 2 we consider b instead of alpha where b^(p^l)=alpha
        if gcd(q_k, n)>1:
            (n, p_to_l) = make_gcd_one(q_k, n)
            alpha = alpha.nth_root(p_to_l)

        # Determine ord(alpha)
        e = alpha.multiplicative_order()
            
        # Determine F_{q^s}, the extension field such that rad(n)|(q^s-1) and 4 !| n or q^s = 1 mod 4
        w = Mod(q_k,rad(n)).multiplicative_order()
        if (not n%4==0) or q_k^w%4 ==1:
            s = w
        else:
            s = 2*w

        if s>1:
            RF_qks = RichExtensionField(RF_qk, s , "y", "T")        
            alpha = RF_qks.to_extf(alpha)
        else:
            RF_qks = RF_qk

        # Determine n_1, n_2, d_1, d_2
        q_ks = RF_qks.q
        n_1 = 1
        n_2 = 1
        for fac in factor(n):
            if e%fac[0] ==0:
                n_1 = n_1 * (fac[0])^(fac[1])
            else:
                n_2 = n_2 * (fac[0])^(fac[1])
        
        d_1s = gcd(n_1, int((q_ks-1)/e))
        d_11 = gcd(n_1, int((q_k-1)/e))        

        # Determine s_1:
        if (not d_1s%4==0) or q_k%4 ==1:
            s_1 = int(d_1s/d_11)
        else:
            s_1= int(d_1s/(d_11*2^(min([(Integer(n/4)).valuation(2), (Integer((q_k+1)/2)).valuation(2)]))))        
        
        #Determine d_2^s:
        d_2s = gcd(n_2, q_ks-1)
        s_2 = Mod(q_k, d_2s).multiplicative_order()      

        # Determine T and d_2t s:
        ts = s_2.divisors()
        ts.sort()
        d2sd2ts = dict()
        for l in range(len(ts)):
            d2sd2t = int(d_2s/gcd(d_2s, q_k^(ts[l])-1))
            if d2sd2t not in d2sd2ts.values():
                d2sd2ts[ts[l]] =  d2sd2t
        ts = list(d2sd2ts.keys())
        
        # Determine r:
        if e==1:
            r=1
        else:
            r= n_2.inverse_mod(e*d_1s)

        # Determine gamma:
        #gamma = (RF_qks.F).multiplicative_generator()
        #print("gamma="+str(gamma))
        
        # Find primitive unity roots 
        zeta_d_1s = RF_qks.get_unityroot(d_1s)
        zeta_d_2s = RF_qks.get_unityroot(d_2s)
        zeta_d_1s_e = RF_qks.get_unityroot(d_1s*e)

        # Find beta:
        for l in range(d_1s*e):
            if gcd(l,d_1s*e)==1:
                if zeta_d_1s_e^(l*d_1s) == alpha:
                    beta = zeta_d_1s_e^l
        
        # Compute beta_j:
        bjs = dict()
        js = dict()
        for j in range(d_1s):
            #bj = (zeta_d_1s^j*beta)^r
            bj = (zeta_d_1s^j*beta)^r
            bjs[j] = bj
            js[bj] = j
        
        #--------------------------------------------------------------------
        # Determine a representative system of the equivalence classes of ~:
        if s_1 ==1:
            j_repsystem = set(range(d_1s))
        else:
            j_repsystem = set()
            J = set(range(d_1s))
            while len(J)>0:
                j = J.pop()
                j_repsystem.add(j)
                J = J.difference(set([js[bjs[j]^(q_k^l)] for l in range(s_1)]))
        
        # Cyclotomic representatives and t_is:
        if s_2 ==1:
            i_repsystem = set(range(d_2s))
            tis = dict([(i,1) for i in i_repsystem])
            cis = dict([(i,s_1) for i in i_repsystem])
        else:
            i_repsystem = set()
            tis = dict()
            cis = dict()
            I = set(range(0,d_2s))
            while len(I)>0:
                i= I.pop()                        
                i_repsystem.add(i)
                for l in range(len(ts)):
                    if i%d2sd2ts[ts[l]]==0:
                        tis[i] = ts[l]
                        cis[i] = lcm(s_1, ts[l])
                        break
                I= I.difference(set([(i*q^(k*l))%d_2s for l in range(tis[i])]))

        n1d1 = int(n_1/d_1s)
        irred_factors = [] 

        if s>1:
            if k>1:
                casting = 3
            else:
                casting =  2
        else:
            if k>1:
                casting =1
            else:
                casting =0 

        if printing:
            print("n="+str(n)+"="+str(factor(n)))
            if gcd(q_k, n)>1:
                print("n = "+str(n)+"*"+str(p_to_l))
                print("beta = "+str(alpha)+" satisfies beta^"+str(p_to_l)+"="+str(alpha^p_to_l))
            print("k="+str(k))
            print("q^k = "+str(q_k)+"="+str(factor(q_k)))
            print("F_q^k: "+str(RF_qk))
            print("alpha="+str(alpha))
            print("ord(alpha) = "+str(e))
            print("rad(n)="+str(factor(rad(n))))
            print("rad(n)|q^w-1 for w="+str(w))
            print("s = "+str(s))
            print("q^s-1 = "+str(RF_qk.q^s-1)+" = "+str(factor(RF_qk.q^s-1)))        
            print("n_1 = "+str(n_1)+"="+str(factor(n_1)))
            print("n_2 = "+str(n_2)+"="+str(factor(n_2)))
            print("d_1s = "+str(d_1s)+"="+str(factor(d_1s)))
            print("d_11 = "+str(d_11)+"="+str(factor(d_11)))
            print("s_1="+str(s_1))
            print("d_2s="+str(d_2s)+"="+str(factor(d_2s)))
            print("s_2 = "+str(s_2))  
            print("d_2s/d_2t : "+str(d2sd2ts))
            print("r="+str(r))
            print(RF_qks)
            print("zeta_{d_1s}="+str(zeta_d_1s))
            print("zeta_{d_2s}="+str(zeta_d_2s))
            print("zeta_{d_1s_e}="+str(zeta_d_1s_e))
            print("beta="+str(beta))
            print("bjs = "+str(bjs))
            print("j_reps : "+str(j_repsystem))
            print("i_reps : "+str(i_repsystem))
            for i in i_repsystem:
                print("i="+str(i))
                print("t_i="+str(tis[i]))
                print("c_i="+str(cis[i]))
                print("C_q,d_2s(i)="+str(set([(i*q^(k*l))%d_2s for l in range(tis[i])])))
            print("casting method: "+str(casting))
        
        # Product over j, v, i, m
        for j in j_repsystem:
            for v in divisors(int(n_2/d_2s)):
                for i in i_repsystem:
                    if gcd(i,v)==1:
                        for m in range(gcd(s_1,tis[i])): 
                                                  
                            fact = spin_linear_factor(q,RF_qks.x,zeta_d_2s^(i*q_k^m)*bjs[j]^v,cis[i]*k)                            

                            fact = cast_from_qks_to_q(casting, fact, RF_qks, RF_qk, RF_q.x, n1d1*v)
                            
                            if printing:
                                print("(j,v,i,m)=("+str(j)+","+str(v)+","+str(i)+","+str(m)+")\n\t"+str(fact)) 
                                #print("("+str(j)+","+str(v)+","+str(i)+","+str(m)+")&: "+latex_pol(fact,"X"))
                            irred_factors.append((fact,p_to_l))

    return Factorization(irred_factors)