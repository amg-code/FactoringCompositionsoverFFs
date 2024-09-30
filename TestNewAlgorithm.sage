# *************************************************************************
#       Copyright (C) 2023 Anna-Maurin Graner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *************************************************************************

# With this program the implementation fXnAlgorithm of the new factorization algorithm for polynomials of the form f(X^n) by Anna-Maurin Graner can be tested and compared with the existing SageMath function factor

from datetime import datetime
import time 
import multiprocessing
import signal
from RichPolynomialClass import *
from RichFiniteFieldClass import *
from fXnAlgorithm import *

def main():
    print("This program tests the implementation of the new factorization algorithm for polynomials of the form f(X^n) by Anna-Maurin Graner and compares it with the existing SageMath function factor (if you wish).")

    # USER INPUT:
    filepath = "/put/your/filepath/here/"
    fileID = "your_file_id"
    sep = "\t"

    timeout = 30 # timeout of the computation after timeout seconds
    comparison = True # comparison with the SageMath algorithm (PARI)?
    print_factorization = True # Shall the factorization be printed?

    parallel_computations = False  # Do you wish to do the computations simultaneously or not? 
    # IMPORTANT: If the computation seems to run endlessly or your computer is out of RAM, then set parallel_computations=False. Otherwise, set parallel_computations=True

    # Polynomials given as follows: (q,n,f), where  f is the coefficient vector of an irreducible polynomial over F_q = F_p(c)
    # f = X^2+c <-> f=[1,0,"c"]

    qnfs = [(16,11^3,[1,"-c^3"]), (2,19^2,[1,1,1,0,1,0,1])]

    #[(2, 2^3*7*3^3*11^2, [1,1,1,0,1,0,1]), (2,19^2,[1,1,1,0,1,0,1]),(2,19^3,[1,1,1,0,1,0,1]), (2,19^4,[1,1,1,0,1,0,1]),(2,19^5,[1,1,1,0,1,0,1]), (2,19^6,[1,1,1,0,1,0,1])]

    #[(4, 2^3*7*3^3*11^2, [1,0,0,1,0,0,"c"]),(4,19^2,[1,0,0,1,0,0,"c"]),(4,19^3,[1,0,0,1,0,0,"c"]), (4,19^4,[1,0,0,1,0,0,"c"]),(4,19^5,[1,0,0,1,0,0,"c"]),(4,19^6,[1,0,0,1,0,0,"c"])] 
 
    #[(16,11^3,[1,"-c^3"]),(16,1321,[1,"-c^3"]),(16,4513,[1,"-c^3"]),(16, 4177, [1,"-c^3"]),(16,3^4*5^4,[1,"-c^3"]),(16,3^5*5^5,[1,"-c^3"]),(16,3^6*5^6,[1,"-c^3"]),(16,3^7*5^7,[1,"-c^3"])]

    # END of USER INPUT. Do not change anything after here.

    print("\ntimeout after: \t\t\t"+str(timeout)+" seconds\nparallel computations: \t\t"+str(parallel_computations)+"\ncomparison with SageMath: \t"+str(True)+"\nprinting factorization: \t"+str(print_factorization)+"\n\nPlease be patient - measurements are in progress.\n")

    if comparison:
        print("Info: Due to the comparison with the SageMath function, the polynomial f(X^n) needs to be computed explicitly and this can take a long time. Don't worry if you see an error thrown by the SageMath function (PARI). Unfortunately, we cannot suppress all error messages but the program will terminate properly.\n")

    measurements(qnfs, timeout, comparison, sep, filepath, fileID, parallel_computations,print_factorization)

    return  

def computation_time(new, n, f):
    if not new:
        f = composition_Xn(f.list,(f.get_RFF()).x,n)
    wst1 = time.time()
    cst1 = time.process_time() 
    if new:
        factorization = factor_fXn(f,n,0)
    else:
        try:
            factorization = f.factor()
        except NotImplementedError:
            factorization = False
    wet1 = (time.time()-wst1)  # in seconds
    cet1 = time.process_time()-cst1
    if factorization==False:
        wet1="stackovf"
        cet1="stackovf"
    return (cet1, wet1, factorization)
    
def computation_time_to_Queue(new, n, f, que):
    que.put(computation_time(new, n, f))
    return

def computation_time_max_mprocessing(new, n, f, seconds):
    q = multiprocessing.Queue()
    p = multiprocessing.Process(target=computation_time_to_Queue, name="computation_time_max_fXn_mprocessing", args=(new, n, f, q,))
    p.start()
    p.join(timeout=seconds)
    if not q.empty():
        result = q.get()
    else:
        result = ("infty","infty",False)
    p.terminate()    
    return result

class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException

def computation_time_max_signal(new, n, f, seconds):
    # Change the behavior of SIGALRM
    signal.signal(signal.SIGALRM, timeout_handler)

    signal.alarm(seconds)
    try:
        result = computation_time(new, n, f)
    except TimeoutException:
        result = ("infty", "infty", False)
    else:
        signal.alarm(0)
    return result

def measure_tupl(q,n,f,seconds, comparison, multiprocess, printing):
    RFF = RichFiniteField(q,"c","X")
    c = RFF.gen
    f = RichPolynomial(f, RFF)
    
    n_tilde = make_gcd_one(q,n)[0]
    w=Mod(q^f.deg,rad(n_tilde)).multiplicative_order()
    if not n_tilde%4==0 or q^(f.deg*w)==1:
        s=w
    else:
        s=2*w

    string = "\n_________________________________\nq="+str(q)+", n="+str(factor(n))+", f="+str(f.get_polynomial())+",  w="+str(w)+", s="+str(s)+"\n"

    if multiprocess:
        (ctn, wtn, fac) = computation_time_max_signal(True, n, f, seconds)
    else:
        (ctn, wtn, fac) = computation_time_max_mprocessing(True, n, f, seconds)
    
    if fac==False:
            print(string+"The new AMG algorithm did not yield a result after "+str(seconds)+" seconds.")
    else:
        if printing:
            string2= string+"f(X^n)=\n"+str_magma_style(fac)+"\n\n"
        else: 
            string2=string
        print(string2+"New AMG Algorithm: \n\tCPU time: \t"+str(ctn)+" seconds\n\tWall time: \t"+str(wtn)+" seconds")
    
    if comparison:
        if multiprocess:
            (cto, wto, faco) = computation_time_max_signal(False, n, f, seconds)
            string2 = string
        else:
            (cto, wto, faco) = computation_time_max_mprocessing(False, n, f, seconds)
            string2 = ""
        
        if cto == "stackovf":
            print(string2+"The SageMath function factor() exited due to a \n<<PariError: the PARI stack overflows>>.")
        elif faco==False:
            print(string2+"The SageMath function factor() did not yield a result after "+str(seconds)+" seconds.")

        if not fac==False and not faco==False:
            if not fac==faco:
                raise SystemError("\nThe factorizations of the new and the SageMath Algorithm are not equal. Something must have gone wrong!")
            else:
                rno = int(round(ctn/cto,0))
                ron = int(round(cto/ctn,0))
        elif faco==False and not fac==False:
            rno=0
            ron="infty"
        elif fac==False and faco==False:
            rno=False
            ron = False
        elif fac==False and not faco==False:
            rno="infty"
            ron =0 
            fac=faco
            if printing:
                string2 = string2+"\nf(X^n)=\n"+str_magma_style(fac)+"\n"
        
        if not faco==False:
            print(string2+"SageMath Algorithm: \n\tCPU time: \t"+str(cto)+" seconds\n\tWall time: \t"+str(wto)+" seconds\nratio NA:SM = "+str(rno)+":1\nratio SM:NA = "+str(ron)+":1\n")
            
    if not fac==False:
        max_deg = ((fac[-1])[0]).degree()
    else:
        max_deg = False

    result_str = "\n"+str(q)+sep+str(n)+sep+str(factor(n))+sep+str(f.list)+sep+str(f.get_order())+sep+str(w)+sep+str(s)+sep+str(max_deg)+sep+str(ctn)

    if comparison:
        result_str = result_str+sep+str(cto)+sep+str(rno)+" : 1"+sep+str(ron)+" : 1"
    if printing:
        result_str=result_str+sep+str(fac)

    return result_str

def measurements(qnfs, seconds, comparison, sep, filepath, filename, multiprocess, printing):
    filehead = "NewfXnAlg_measurements_"
    file = open(filepath+filehead+filename+".csv", "w")
    file.write("Computation time measurements \nfor the new f(X^n) factorization algorithm \nby Anna-Maurin Graner\n\n"+str(datetime.now())+"\n\nq"+sep+"n"+sep+"factor(n)"+sep+"f"+sep+"ord(f)"+sep+"w"+sep+"s"+sep+"max deg"+sep+"New Alg (CPU time)") 
    if comparison: 
        file.write(sep+"SageMath Alg (CPU time)"+sep+"NA:SM"+sep+"SM:NA")
    if printing:
        file.write(sep+"Factorization f(X^n)")

    if multiprocess:
        pool = multiprocessing.Pool()
        
        for st in pool.starmap(measure_tupl,[(q,n,f,seconds,comparison,multiprocess,printing) for (q,n,f) in qnfs]):
            file.write(st)

        pool.close()
    else:
        for (q,n,f) in qnfs:
            file.write(measure_tupl(q,n,f,seconds,comparison,multiprocess,printing))

    file.close()
    print("\nThe results of the measurements can be found in "+filepath+filehead+filename+".csv.")
    return  

if __name__ == "__main__":
    main()