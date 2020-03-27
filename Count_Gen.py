#!/usr/bin/env python
# -*- coding: utf-8 -*-
import random, time, math
import sys as sys
from Pool import *
from Utils import *

memo = dict()

##############################################################################
def add(t, f):
    """
    add two lists t <- t + f 
    len(t) MUST BE >= len(f)
    """
    lf = len(f)
    l = len(t)
    if l < lf:
        print("pb")
        print("{} + {}".format(t,f))
    for i in range(lf):
        t[i] += f[i]
########################################################################
def update_dict(R, key, value):
    if key in R:
        R[key] += value
    else:
        R[key] = value
#########################################################################
def count(n, p, s):
    """
    n: size of the enumerated pool-BDDs
    p: profile of the pool (the number of variables k is len(p)-1+1)
    s: sibling rank of root of the pool-BDD
    """
    
    if (n, p, s) in memo: # memoization
        return memo[(n, p, s)]
    k = len(p)
    d = {}

    for i in range(n): # Left and right subtree of size $i$ and $n-i-1$
        d0 = {}
        if i == 0:
            d0[tuple()] = sum(p)
        else:
            for k0 in range(1, k): # left node has index $k_0$
                d0.update(count(i, tuple(p[:k0]), p[k0]))
        for l in d0:
            w0 = d0[l]
            d1 = {}
            pp = list(p)
            add(pp, l)
            if (n-1-i == 0): # right subtree is empty
                d1[tuple()] = sum(pp) -1
            else:
                for k1 in range(1, k): # right node has index $k_0$
                    d1.update(count(n-1-i, tuple(pp[:k1]), pp[k1]))
            for r in d1:
                w1 = d1[r]
                t = [0 for _ in range(k)]
                add(t, l)
                add(t, r)
                t.append(1)
                w = w0 * w1
                if n == 1:
                    w =  w-s
                if w > 0:
                    update_dict(d, tuple(t), w)

    memo[(n, p, s)] = d
    return d
###############################################################################
def unrank_pair(r, A, F):
    """
    r : rank inside all pairs (x, y) from AxA such that x!=y and (x, y) does not belong to F
    X: list of elements
    F: list of forbidden pairs of elements 
    returns list F with pair (x, y) inserted at the "right" place and pair (x, y) of rank r
    """
    t = r
    i = 0
    for x in A:
        for y in A:
            if x != y:
                if i < len(F) and (x,y) == F[i]:
                    i += 1
                else:
                    t -= 1
                if t == -1:
                    F = F[:i]+[(x, y)]+F[i:]
                    return F, x, y
    print("error: unrank_pair({}, {}, {})".format(r, A, F))
    
##############################################################################
def decompose(rank, n, pool, target_profile):
    """
    rank: rank within structures with profile target_profile
    n: size of the bdd
    pool_profile:
    target_profile: 
    returns a tuple (r, i, h_left, profile, len_profile) where
    * r is the rank of the left bdd 
    * i is the size of the left bdd
    * h_left is the height of the left bdd
    * profile is the profile of the left bdd
    * ind_max is the index of the node wrt "rooting" the decomposition
    """
    r = rank
    k = len(target_profile)-1
    total_profile = pool.get_profile()
    p, sibling_rank = total_profile[:k], total_profile[k]
    # check bdd with left child bdd having size n-1
    i = n - 1
    while i >= 0:
        d0 = {}
        if i == 0:
            d0[tuple([0]*k)] = sum(p)
        else:
            for k0 in range(1, k): # left node has index $k_0$
                d0.update(count(i, tuple(p[:k0]), p[k0]))
        for lo_profile in d0:
            w0 = d0[lo_profile]
            d1 = {}
            pp = list(p)
            add(pp, lo_profile)
            if (n-1-i == 0): # right subtree is empty
                d1[tuple([0]*k)] = sum(pp) -1
            else:
                for k1 in range(1, k): # right node has index $k_0$
                    d1.update(count(n-1-i, tuple(pp[:k1]), pp[k1]))
            for hi_profile in d1:
                w1 = d1[hi_profile]
                tree_profile = [0 for _ in range(k)]
                add(tree_profile, lo_profile)
                add(tree_profile, hi_profile)
                tree_profile.append(1)
                tree_profile=tuple(tree_profile)
                w = w0 * w1
                if n == 1:
                    w =  w-sibling_rank
                if tree_profile == target_profile:
                    if w > 0:
                        r -= w
                    if r < 0:
                        r += w
                        return (r, i, lo_profile, hi_profile)
        i -= 1
    
##########################################################################
def get_profile(pool):
    profile = []
    for l in pool:
        profile.append(len(l))
    return tuple(profile)
########################################################################
def generate(rank, n, pool, target_profile):
    """
    rank: rank within structures with profile target_profile
    n: size of the bdd (in facts n+2 for the two constants 0 and 1)
    pool: 
    target_profile : profil à atteindre
    """
    k = len(target_profile)-1
    if n ==0:
        i_node = pool.get_nth(rank)
    elif n == 1:
        lo, hi = pool.unrank_pair(rank, k)
        i_node = pool.add_node(lo, hi, k)
    else:
        r, i, lo_profile, hi_profile = decompose(rank, n, pool, target_profile)
        # bdd with left bdd of size i>0 and right bdd of size n-1-i>0
        p = pool.get_profile()
        lo_index = len(lo_profile)-1
        D=count(i, p[:lo_index], p[lo_index])
        max_left = D[lo_profile]
        r_left = r % max_left
        r_right = r // max_left
        # generate left bdd
        i_lo = generate(r_left, i,  pool, lo_profile)
        if n-i-1 == 0 and pool.get_rank(i_lo) <= r_right:
            r_right += 1
        i_hi = generate(r_right, n-1-i, pool, hi_profile)
        i_node = pool.add_node(i_lo, i_hi, k)
    return  i_node
#################################################################################
def rand_bdd(n, k, rank = None):
    sibling_rank = 0
    d = count(n,tuple([2]+[0]*(k-1)), 0)
    total = sum(d.values())
    r = rank
    if r == None:
        r = random.randint(0,total-1)
    elif r >= total:
        r = rank % total
    for tree_profile in d:
        r -= d[tree_profile]
        if r < 0:
            r += d[tree_profile]
            pool = Pool(tree_profile)
            generate(r, n, pool, tree_profile)
            return pool
    return None
###########################################################################
def gen_bdd(size, k):
    sibling_rank = 0
    d = count(n,tuple([2]+[0]*(k-1)), 0)
    for tree_profile in d:
        for rank in range(d[tree_profile]):
            pool = Pool(tree_profile) 
            generate(rank, size, pool, tree_profile)
            yield pool

###########################################################################
def printAll(n, nb_vars):
    r =0
    for b in gen_bdd(n, nb_vars):
        print("*** List({}):\n{}".format(r,b.to_list()))
        r += 1
###########################################################################



## n: size of the BDD (at least 3 for True and False)
n = int(sys.argv[1])
## k: number of variables, the root of the BDD is labeled by x_k; x_k is an essential variable
k = int(sys.argv[2])

## enumerate for each valid profile the number of BDDs having this profile
d = count(n-2, tuple([2]+[0]*(k-1)), 0)

## print the result of the enumeration
## the pool-profiles that are printed to not contain the profile of the pool
## thus for example, the 2 constant True and False do not appear in the first component that is thus reduced to 0
for t in sorted(d):
    print("{} : {}".format(t,d[t]))


## print all BDDs in list-format
printAll(n-2, k)
## number of elements stored for the memoization
print("len(memo)={}".format(len(memo)))
## number of pool-BBDs of size n and  
print("#pool-BDD = {} (#profiles={})".format(sum(d.values()), len(d)))


# sample uniformly at random a BDD of size n and k variables
b = rand_bdd(n-2, k)
# save the dot file corresponding to the sampled BDD
dot_save(b.to_list(), n, k, 'rand')