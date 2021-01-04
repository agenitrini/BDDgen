#!/usr/bin/env python
# -*- coding: utf-8 -*-
import random, time, math
import sys as sys
from Spine import *
from Utils import *
memo = dict()
#######################################################################
def too_big(n, k, p, s):
    S = sum(p)+s+n
    return n > MAX_SIZE[k] or S > CUMUL_MAX[k]

#########################################################################
def count(n, p, s):
    k = len(p)+1
    L = []
    if too_big(n, k, p, s):
        return []
    if n == 0: # useful only for decompose
        L = [(tuple([0]*k), sum(p)+s+2)]
    elif n == 1: #external
        S = sum(p)+2
        S = S*(S-1)-s
        if S > 0:
            L = [(tuple([0]*(k-1)+[1]),  S)]
    elif (n,p) in memo: # memoization
        L = memo[(n, p)]
    else:
        d = {}
        for i in range(n): # Left and right subtree of size $i$ and $n-i-1$
            d0 = []
            if i == 0:
                d0 = [(tuple(), sum(p)+2)]
            else:
                for k0 in range(0, k-1): # left node has index $k_0$
                    d0.extend(count(i, tuple(p[:k0]), p[k0]))
            if len(d0) == 0: # no need to look at greater values of i
                break
            for l, w0 in d0:
                d1 = []
                pp = list(p)
                add(pp, l)
                if (n-1-i == 0): # right subtree is empty
                    d1  = [(tuple(),sum(pp)+2-1)]
                else:
                    for k1 in range(0, k-1): # right node has index $k_0$
                        d1.extend(count(n-1-i, tuple(pp[:k1]), pp[k1]))
                for r, w1 in d1:
                    t = [0 for _ in range(k)]
                    add(t, l)
                    add(t, r)
                    t[-1]= 1
                    w = w0 * w1
                    if w> 0:
                        update_dict(d, tuple(t), w)
        L = [(key, d[key]) for key in d]
        if len(L) >0:
            memo[(n, p)] = L
    return L
##############################################################################
def decompose(rank, n, spine, target_profile):
    """
    rank: rank within structures with profile target_profile
    n: size of the bdd
    spine: current spine
    target_profile: target spine profile  
    returns a tuple (lo_size, r_left, lo_profile, hi_size, r_right, hi_profile) where
    * lo_size is the size of the left bdd
    * r_left is the rank of the left bdd 
    * lo_profile is the profile of the left bdd 
    * hi_size is the size of the right bdd
    * r_right is the rank of the right bdd 
    * hi_profile is the profile of the right bdd 
    """
    r = rank
    k = len(target_profile)
    total_profile = spine.get_spine_profile()
    p = total_profile[:k-1]
    #    rank_sibling = total_profile[k-1] # not useful
    i = n - 1
    while i >= 0:
        if i == 0:
            d0 = [(tuple([0]*(k-1)),sum(p)+2)]
        else:
            d0 = []
            for k0 in range(0, k-1): # left node has index $k_0$
                d0.extend(count(i, tuple(p[:k0]), p[k0]))
        for lo_profile, w0 in d0:
            pp = list(p)
            add(pp, lo_profile)
            if (n-1-i == 0): # right subtree is empty
                d1 = [(tuple([0]*(k-1)), sum(pp) +2 -1)]
            else:
                d1 = []
                for k1 in range(0, k-1): # right node has index $k_0$
                    d1.extend(count(n-1-i, tuple(pp[:k1]), pp[k1]))
            for hi_profile, w1 in d1:
                tree_profile = [0 for _ in range(k)]
                add(tree_profile, lo_profile)
                add(tree_profile, hi_profile)
                tree_profile[-1] = 1
                tree_profile=tuple(tree_profile)
                w = w0 * w1
                if tree_profile == target_profile:
                    if w > 0:
                        r -= w
                    if r < 0:
                        r += w
                        r_left = r % w0
                        r_right = r // w0
                        lo_size = i
                        hi_size = n-1-i
                        return (lo_size, r_left, lo_profile, hi_size, r_right, hi_profile)
        i -= 1
########################################################################
def generate(rank, n, spine, target_profile):
    """
    rank: rank within BDDs with spine profile target_profile
    n: size of the spine
    spine: current spine
    target_profile : spine profile to achieve
    """
    k = len(target_profile)
    if n ==0:
        i_node = spine.unrank_singleton(rank)
    elif n == 1:
        lo, hi = spine.unrank_pair(rank, k)
        i_node = spine.add_node(lo, hi, k)
    else:
        lo_size, r_left, lo_profile, hi_size, r_right, hi_profile = decompose(rank, n, spine, target_profile)
        # bdd with left bdd of size i>0 and right bdd of size n-1-i>0
        i_lo = generate(r_left, lo_size,  spine, lo_profile)
        if hi_size == 0 and spine.get_rank(i_lo) <= r_right:
            r_right += 1
            print("Problem : i_lo={} rank(i_lo)={}, r_right={}\nspine={}".format(i_lo,spine.get_rank(i_lo), r_right, spine))
        i_hi = generate(r_right, hi_size, spine, hi_profile)
        i_node = spine.add_node(i_lo, i_hi, k)
    return  i_node
#################################################################################
def rand_bdd(n, k, rank = None):
    sibling_rank = 0
    d = count(n,tuple([0]*(k-1)), 0)
    total = sum([x for _, x in d])
    r = rank
    if r == None:
        r = random.randint(0,total-1)
    elif r >= total:
        r = rank % total
    for tree_profile, w in d:
        r -= w
        if r < 0:
            r += w
            spine = Spine(tree_profile)
            generate(r, n, spine, tree_profile)
            return spine
    return None
#################################################################################
def gen_bdd(size, k):
    sibling_rank = 0
    d = count(n,tuple([0]*(k-1)), 0)
    for tree_profile, w in d:
        for rank in range(w):
            spine = Spine(tree_profile) 
            generate(rank, size, spine, tree_profile)
            yield spine
            
####################################
def printAll(n, k):
    d = count(n,tuple([0]*(k-1)), 0)
    total = sum([x for _,x in d])
    r =0
    for b in gen_bdd(n, k):
        print("*** {}:\n{}".format(r, b.to_list()))
        r += 1
####################################
def pgcd(a,b):
    """pgcd(a,b): calcul du 'Plus Grand Commun Diviseur' entre les 2 nombres entiers a et b"""
    while b != 0:
        r= a % b
        a,b=b,r
    return a
####################################
def pgcd_l(L):
    if len(L) == 0:
        return 1
    gcd = L[0]
    for x in L[1:]:
        gcd = pgcd(gcd, x)
    return gcd
###################################################
n = int(sys.argv[1])
k = int(sys.argv[2])

MAX_SIZE = [sum(max_profile(i)) for i in range(k+1)]
MAX_PROFILE = max_profile(k)
CUMUL_MAX = [ sum(MAX_PROFILE[:i]) for i in range(len(MAX_SIZE)+1)]

print("# MAX_SIZE={}\n# CUMUL_MAX={}\n".format(MAX_SIZE, CUMUL_MAX))

if n > MAX_SIZE[k]:
    print("n={} is too big for k={} variables!!! setting n={}".format(n, k, MAX_SIZE[k]))
    n = MAX_SIZE[k]
d = count(n, tuple([0]*(k-1)), 0)
for t,w in sorted(d):
    print("# {} : {}".format(t,w))

#print("# =========================")
#printAll(n, k)
#print("#BDD = {} (#profiles={})".format(sum([w for _, w in d]), len(d)))
     
# sample uniformly at random a BDD of size n and k variables
b = rand_bdd(n, k)
print(b.to_list())
# save the dot file corresponding to the sampled BDD
dot_save(b.to_list(), n, k, 'rand')

