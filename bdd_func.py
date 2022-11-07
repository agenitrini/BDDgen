# -*- coding: utf-8 -*-
# simplest recurrence for spines

from math import *
import sys as sys
import time
import functools
from timeit import default_timer as timer
#import unrank_bdd as check
from itertools import zip_longest

########################################
def calc(p, m):
    """
    caution: profile p is incomplete (without the last value 2 for constants
    a is the profile of the multiset entries
    """
    P = (1,)
    for x, y in zip_longest(p, m, fillvalue=0): # padding if needed
        P = iter(x, (0,)*y + P)
    return P
########################################
@functools.lru_cache(maxsize=None)
def count_multi_BDDs(p, m):
    k = len(p)
    s = calc(p, m)
    return eval(s)

###############################
#@functools.lru_cache(maxsize=None)
def eval(P):
    i = 1
    s = 0
    for x in P:
        s += i*x
        i *= 2
    return s
##############################
def add(x, y):
    return tuple(a+b for a, b in zip_longest(x, y, fillvalue=0))
##############################
def prod(x, y):
    dx = len(x)-1
    dy = len(y)-1
    L = [0 for _ in range(dx+dy+1)]
    for i in range(dx+1):
        for j in range(dy+1):
            L[i+j] += x[i]*y[j]
    return tuple(L)
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
##############################
@functools.lru_cache(maxsize=None)
def R(n, k):
    #    Product((x^2-x-i), i=0..n-1)
    p = (1,)
    for i in range(n):
        p = prod(p, (-i, -1, 1))
    q = tuple(binomial(k, j)*stirling2(k-j, n) for j in range(k-n+1))
    return prod(p,q)
###############################
#@functools.lru_cache(maxsize=None)
def iter(n, P):
    Q = tuple()
    for k in range(n, len(P)):
        coeff = P[k]
        #        if coeff != 0: # does not happen often, not worth the test
        Q = add(Q,  tuple(coeff*x for x in R(n, k)))
    return Q
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
######################################
def unrank_profile(rank,  n, nb_vars):
    #print("n={} nb_vars={} rank={}".format(n, nb_vars, rank))
    r = rank
    p =()
    for layer in range(len(p), nb_vars):
        top = sum(p)
        bottom = n - top
        """
        dichotomic search
        search mid such that count(..mid-1) =< rank and count(..mid) > rank
        """
        low = 0
        high = 1 +min(top+1, bottom-ceil(sqrt(bottom)))
        mid = 0
        P = 0
        while low + 1 < high:
            mid = (high + low) // 2
            res = layers_build(n, nb_vars, p, constraint=(low, mid))
            if len(res)<= n-2:
                nb = 0
            else:
                nb = eval(res[n-2])
            if P + nb <= r:
                P += nb
                low = mid
            elif P + nb > r:
                high = mid
        p = p + (low,)
        r -= P
    return r, p
###########################################
def layers_init(p):
    P = (0,1)
    for x in p:
        # if x > 0: # does not happen often, not worth the test
        P = iter(x, P)
    return ((),)*sum(p)+(P,)
######################################
def layers_add(n, nb_vars, current, constraint=None):
    max_top = 2*(len(current)-1)+1
    m = min(max_top+1, n-2)
    next = [()]*(m+1)
    low = 0
    high = len(current)+1
    if constraint is not None:
        low, high = constraint
    for top in range(0, len(current)):
        bottom = n-top
        m = min(top+1, bottom-ceil(sqrt(bottom)))
        for x in range(low, min(m+1, high)):
            if x + top <= n-2:
                # P = current[top]
                # # if x > 0: # does not happen often, not worth the test
                # P = iter(x, P)
                next[top+x] = add(next[top+x], iter(x, current[top]))
    return next
#####################################
def layers_build(n, nb_vars, p, constraint=None):
    """
    add layer from start to nb_vars-1
    """
    current = layers_init(p)
    len_prefix = len(p)
    if constraint is not None:
        current = layers_add(n, nb_vars, current, constraint)
        len_prefix += 1
    for layer in range(len_prefix, nb_vars):
        current = layers_add(n, nb_vars, current)
    return current
###################
def count(n, nb_vars, p=()):
    res = layers_build(n, nb_vars, p)
    if len(res) <= n-2:
        return 0
    return eval(res[n-2])
######################################
def count_all_size(n, nb_vars, p=()):
    start_total = timer()
    print("#n={}, nb_vars={}, initial p={}".format(n, nb_vars, p))
    T = sum(p)
    previous = ()
    current = layers_init(p)
    for layer in range(len(p), nb_vars):
        S = min(2*T+1, n-2) # max size possible with new layer T +(T+1)
        next = [()]*(S+1)
        print("#layer={}".format(layer),end="", flush=True)
        start = timer()
        for top in range(0, T+1): #len(current)):
            print(".", end="", flush=True)
            bottom = n-top
            m = min(top+1, bottom-ceil(sqrt(bottom))) # max #nodes to add
            for x in range(0, m+1): # starting at 0 is important
                next[top+x] = add(next[top+x], iter(x, current[top]))
        end = timer()
        print("#\n#time for layer {}: {} s".format(layer, end - start), flush=True)
        current = next
        T = S
        #print("POST: len(current) = {} T={} S={}".format(len(current), T, S))
        #print("current={}".format(current))
    end = timer()
    print("#\n#total time: {} s".format(end - start_total), flush=True)
    return tuple(eval(y) for y in current)
###############################################################
if __name__ == '__main__':

    if len(sys.argv)== 1:
        print("arguments are missing")
        exit()
    k = int(sys.argv[1])
    n = max_size(k)
    if len(sys.argv) > 2:
        n = int(sys.argv[2])
        if n > max_size(k):
            print("max size is {}".format(max_size(k)))
            n = max_size(k)
    p=() #p=(1,) # if we want ROBDDs with a root in layer k
    if len(sys.argv) > 3:
        p =  tuple(int(x) for x in sys.argv[3].split())
    print("#  {} vars, max size = {}, initial profile={}\n#nodes\t#bdds".format(k, max_size(k), p))
    res = count_all_size(n, k, p)
    for i in range(len(res)):
        print("{}\t{}".format(i+2, res[i]))
    print("# {}".format(R.cache_info()))
    print("# total = {} ROBDDs\t({} boolean functions)".format(sum(res), 2**(2**k)))
    M = max(res)
    print("# maximum at layer {}:\n#\tMax={}\n#\tproportion={}%".format(res.index(M)+2, M, M/2**(2**k)))
