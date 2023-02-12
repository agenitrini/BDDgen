# -*- coding: utf-8 -*-
# simplest recurrence for spines

from math import sqrt, ceil
from poly_func import prod_opt, prod_trunc, add, deg, prod_naive, prod_val, prod_val_naive
import sys as sys
import functools
from timeit import default_timer as timer
#import unrank_bdd as check
TRUNC = ()
CUMUL_TRUNC = ()
stats = {}
##############################
@functools.lru_cache(maxsize=None)
def binomial(n, k):
    # since C(n, k) = C(n, n - k)
    if(k > n - k):
        k = n - k
    # initialize result
    res = 1
    # Calculate value of
    # [n * (n-1) *---* (n-k + 1)] / [k * (k-1) *----* 1]
    for i in range(k):
        res = res * (n - i)
        res = res // (i + 1)
    return res
##############################
@functools.lru_cache(maxsize=None)
def stirling2(n, k):
    """
    S2(n, 0) and S2(0, k) = 0 # for n, k > 0
    S2(n, n) = 1
    S2(n + 1, k) = k * S2(n, k) + S2(n, k - 1)
    """
    if n == k:
        return 1
    if n == 0 or k == 0:
        return 0
    return k*stirling2(n-1, k)+stirling2(n-1, k-1)
#####################
@functools.lru_cache(maxsize=None)
def _P(n):
    """
    Product((x^2-x-i), i=0..n-1)
    """
    p = (1,)
    for i in range(n):
        p = prod_opt(p, (-i, -1, 1))
    return p
#####################
@functools.lru_cache(maxsize=None)
def _R(n, k):
    return prod_opt(tuple(binomial(k, j)*stirling2(k-j, n) for j in range(k-n+1)), _P(n))
##########################################
def init(n, k):
    global TRUNC
    global CUMUL
    p = max_profile(max_size(k))
    L = tuple(min(x, n) for x in p) 
    TRUNC= tuple(L[::-1])
    print("# TRUNC={}".format(TRUNC))
    L = []
    v = 0
    for i in range(k+1):
        L.append(min(n, v + p[k-i]))
        v = L[-1]
    CUMUL = tuple(L)
    print("# CUMUL={}".format(CUMUL))
##########################################
def max_profile(n):
    top = 0
    bottom = n
    L = []
    while top <n-2:
        m = min(top+1, bottom-ceil(sqrt(bottom)))
        L.append(m)
        top += m
        bottom -= m
    L.append(2)
    return tuple(L)
################################
def max_size(k):
    """
    N.B.
    A327461: Maximal size of a Binary Decision Diagram (or BDD) of index k.
    def A327461(k):
      return 2**(k-(k-k.bit_length()+1).bit_length()+1)+2**2**((k-k.bit_length()+1).bit_length()-1)-1
    # Pontus von Brömssen, Apr 08 2020
    Also:
    n = 2
    for i in range(1,k+1):
        t = min(2**(2**(i-1))*(2**(2**(i-1))-1), 2**(k-i))
        n += t
    return n
    """
    return 2**(k-(k-k.bit_length()+1).bit_length()+1)+2**2**((k-k.bit_length()+1).bit_length()-1)-1

#############################
@functools.lru_cache(maxsize=None)
def _D(n, l):
    """
    varphi(X^l).collect(X) mod X^n
    """
    max_degree = min(2*l, n-1) # max number of nodes that can be attached to l anchors
    P = [(),]*(max_degree+1)
    for x in range(0, min(l+1, n-l)):
        R_X = _R(x, l)
        dR = len(R_X)-1
        for d in range(min(dR, max_degree)+1):
            c = R_X[d]
            if c != 0:
                P[d] = shift_add(P[d], c, x)
    return tuple(P)
###########################
def shift_add(a, c, x):
    """
    return add(a, (0,)*x+(c,))   
    N.B: return sa + c*u^c if a seen as a polynomial in Z[u]
    """
    da = len(a)-1
    if x <= da:
        return a[:x]+(c+a[x],)+a[x+1:]
    return  a + (0,)*(x-da-1)+(c,)
###############################
@functools.lru_cache(maxsize=None)
def expand_monomial(n, l, k):
    """
    return varphi^k(X^l) subs(X=2) mod u^(n)
    """
    global CUMUL
    global TRUNC
    if k == 0:
        return (0, 0, 2**l) # u^2 2**l
    Q =()
    Ru = _D(n, l)
    # print("# n={} l={} k ={}".format(n, l, k))
    # for j, Pu  in enumerate(Ru):
    #     print("#\t j={}\t{}".format(j, Pu))
    #     print("# R(n={}, l={})={} (len={})".format(n, l, R, len(R)))
    #for j, Pu  in enumerate(Ru): # Pu in Z[u]
    deg_P = TRUNC[k]
    deg_M = CUMUL[k-1]
    deg_Q = CUMUL[k]
    for j in range(len(Ru)):
        Pu = Ru[j] #[:deg_P+1]
        Mu = expand_monomial(n, j, k-1) #[:deg_M+1] # Mu in Z[u]
        Q = add(Q, prod_trunc(Pu, min(deg(Pu), deg_P), Mu, min(deg_M, deg(Mu)), deg_Q)) #deg_trunc_u))
    pos = 0
    while pos < len(Q) and Q[pos]>= 0:
        pos += 1
    Q = Q[:pos]
    return Q

###############################################################
if __name__ == '__main__':
    n = -1
    if len(sys.argv) == 3:
        n = int(sys.argv[2])
        k = int(sys.argv[1])
    elif len(sys.argv) == 2:
        k = int(sys.argv[1])
        n = max_size(k)
    if n > max_size(k) or n == -1:
        print("# max size is {} for {} vars!".format(max_size(k), k))
        n = max_size(k)
    print("# #nodes: {}\t#vars:{}".format(n,k))
    print("# max profile for {} nodes: {}".format(n, max_profile(n)))
    init(n, k)
    print("# truncation profile: {}".format(TRUNC))
    start = timer()
    L = enumerate(expand_monomial(n, 1,k))
    end = timer()
    for x, count in L:
        print("{}\t{}".format(x, count))
    print("# done in {}".format(end-start))
    print("# expand_monomial: {}".format(expand_monomial.cache_info()))
    print("# D: {}".format(_D.cache_info()))
    # for (a, b) in sorted(stats):
    #     print("# n={} X^{}".format(a, b))
    #     for j, Pu  in enumerate(stats[a, b]):
    #         print("#\tj={}\t{}".format(j, Pu))
    # print("# maximum at layer {}:\n#\tMax={}\n#\tproportion={}%".format(res.index(M)+2, M, M/2**(2**k)))
