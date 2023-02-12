#-*- coding: utf-8 -*-
import math
import sys as sys
import time
import random as random
import bisect as bisect
import bdd_func_minimal as bdd
import functools
from math import ceil, sqrt
from poly_func import add
from itertools import zip_longest

from timeit import default_timer as timer
dot = True
###################################
def to_TT(L):
    """
    return a string with the table truth
    CAUTION: To debug !!!
    """
    n = len(L)
    bead = dict()
    bead[n-1]="0"
    bead[n-2]="1"
    for i_node in range(n-3, -1, -1):
        if i_node not in bead:
            level, delta  = L[i_node]
            lo, hi = delta
            lo, hi = abs(lo), abs(hi)
            if level > 0:
                lo_bead = bead[lo]
                hi_bead = bead[hi]
                l0, _ = L[lo]
                l1, _ = L[hi]
                d0 = level - l0 - 1
                d1 = level - l1 - 1
                bead[i_node] = lo_bead * (2**d0) +hi_bead*(2**d1)
    return bead[0]
###################################
def getProfile(L):
    x = []
    for l, _ in L:
        x.append(l)
    p = [0 for _ in range(max(x)+1)]
    for i in x:
        p[i] += 1
    return tuple(p[::-1])
###################################
def is_terminal(n):
    return n is None or n < 0 # or n >= len(L)-2 ?
######################
def weight(L):
    """ returns the weight of the profile and
    the list the one of each node that is terminal, in the inorder traversal"""
    p = getProfile(L)
    w = 1
    D = []
    S = []
    current = 0
    pool = [0 for _ in range(len(p))]
    pool[0] = 2
    while len(S) > 0 or not is_terminal(current): # postorder traversal
        if not is_terminal(current):  # if current is not a leaf
            S.append(current)
            level, delta = L[current]
            lo, hi = delta
            current = lo
        else:
            i_node = S.pop()
            level, delta = L[i_node]
            lo, hi = delta
            M = sum(pool[:level])
            if is_terminal(lo):
                if is_terminal(hi):
                    res = (M*(M-1) - pool[level])
                    w *= res
                    D.append(res)
                else:
                    w *= M
                    D.append(M)
            elif is_terminal(hi):
                w *= M - 1
                D.append(M - 1)
            pool[level] += 1 # current node is not new and belongs to the pool
            current = hi
    D.reverse()
    return w, D
###########################
def to_dot(L, k, max_level=None):
    if max_level is None:
        max_level = k
    omit = False
    # choose your colors !!!
    back_edge_color = "red" # "gray" or "red": color of non tree edges
    background_color = "gray" # color of removed part of the robdd (above max_level excluded)
    multi_color= "red"
    def to_dot_aux(root):
        s = ""
        if root >= 0:
            level, delta = L[root]
            lo, hi = delta
            left= ""
            right = ""
            id = "t{}".format(root)
            N[level].append(id)
            constraint= "false" if (lo>0 and hi>0) else "false"
            color = "black" if level <= max_level else background_color
            s += "{} [label=\"{}\", constraint={}, color={}, fontcolor={}]\n".format(id, root, constraint, color, color)
            if abs(lo) < len(L)-2 or not omit:
                id0="t{}".format(abs(lo))
            if abs(lo) < len(L)-2:
                s+= to_dot_aux(lo)
            constraint="true"
            color = "black" if (lo >= 0) else back_edge_color
            level_child, _ = L[abs(lo)]
            if level> max_level and level_child<= max_level: color = multi_color
            if level> max_level and level_child> max_level: color = background_color

            s += "{}:sw -> {} [style=dotted, color={}, constraint={}]\n".format(id,  id0, color, constraint)
            #mid = "m{}".format(root)
            #s +="{} [label=\" \", style=invis]\n{} -> {} [style=invis]\n".format(mid,id, mid)
            if abs(hi) < len(L)-2 or not omit:
                id1 = "t{}".format(abs(hi))
            if abs(hi) < len(L)-2:
                s+= to_dot_aux(hi)
            constraint="true"
            color = "black" if (hi >= 0) else back_edge_color
            level_child, _ = L[abs(hi)]
            if level> max_level and level_child<= max_level: color = multi_color
            if level> max_level and level_child> max_level: color = background_color
            s += "{}:se -> {} [style=filled, color={}, constraint={}]\n".format(id,  id1, color, constraint)
        return s
    N = [[] for _ in range(k+1)]
    s = "digraph G {\nrankdir=TB;\nnewrank=new;\n\nordering=\"out\";\n\nnodesep=0.2;\nranksep=.2;\nmargin=0;\nnode [shape=circle,style=\"\",color=\"black\",margin=0.001,height=0.05,width=0.05,\nlabel=\"\",\npenwidth=1.5,\nfontsize=12]\nedge [arrowsize=.5];\n"
    for level in range(k+1):
        s+= "L{} [label=\"Layer {}\", style=invis,shape=none]\n".format(level, level)
    for level in range(k):
        s+= "L{} -> L{} [style=invis; minlen=2]\n".format(level+1, level)
    s += "subgraph{\n"
    s += to_dot_aux(0)
    if not(omit):
        s += "t{} [label=\"⊥\",margin=0.03, shape = box, constraint=true]\nt{} [label=\"⊤\", margin=0.03, shape = box, constraint=true]\n".format(len(L)-1, len(L)-2)
        N[0].append("t{}".format(len(L)-1))
        N[0].append("t{}".format(len(L)-2))
    s += "}\n"
    for level in range(k+1):
        if len(N[level])>0:
            s +="{{rank=same; L{}".format(level)
            for id in N[level]:
                s +="-> {}".format(id)
            s = s + " [style= invis, constraint={}];\n}}\n".format("true")
    s+="}\n"
    return s
############################################
def abs(x):
    if x >= 0:
        return x
    return -x
############################################
def rank_to_pair(r, M):
    """
    return the r-th pair in {(i, j) : 0 <= i, j < M and i neq j}
    lex order
    """
    R = r+(r//M)+1
    i = R//M
    j = R % M
    return (i, j)
############################################
def pair_to_rank(i, j, M):
    """
    return the rank of (i, j) in {(i, j) : 0 <= i, j < M and i neq j}
    lex order
    (unused)
    """
    r = i*M+j
    R = r-1- ((r-1)//(M+1))
    return R
############################################
def multi_profile(a):
    if len(a) == 0:
        return ()
    m = [0 for _ in range(a[0]+1)]
    for x in a:
        m[x] += 1
    return tuple(m)
##########################################################
def remove_empty_levels(p, a):
    q = []
    offset = [0 for _ in range(len(p)+1)]
    for i in range(len(p)):
        if p[i] == 0:
            offset[i] = 1
        else:
            q.append(p[i])
    for i in range(len(offset)-1):
        offset[i] += offset[i-1]
    A = tuple(x -offset[x] for x in a)
    return tuple(q), A
#####################################
def spine_rank_to_rank_list(L, r):
    _, D = weight(L)
    Dr = []
    for d in D:
        Dr.append(r % d)
        r //= d
    return Dr
#####################################
def adjust(L, D, p):
    """
    Adjust relative pointers to the pool
    *once the spine is complete*
    by a infix traversal of the tree (so the pool is known at each step)
    profile in order (p_0=2, ..., p_k=1)
    """
    ### HELPER FUNCTIONS: current pool and final profile are known
    ##########
    def rank2node(r): #, pool, profile):
        """
        return the index in L of node of rank r in the pool
        """
        level = 0
        while (True):
            Delta = len(pool[level])
            if (r -Delta ) < 0:
                break
            level += 1
            r = r - Delta
        i_node = -pool[level][r]
        return i_node
    #########
    def node2rank(i_node):
        """
        returns the rank of node of index i_node in L relatively to the pool
        """
        level, _ = L[i_node]
        r = sum([len(pool[i]) for i in range(level)])
        return r+bisect.bisect_left(pool[level], -i_node)
    #########
    def node2absrank(i_node):
        """
        returns the absolute rank of node of index i_node in L relatively to the set of all nodes
        """
        level, _ = L[i_node]
        r = sum(p[:level])
        return r+bisect.bisect_left(pool[level], -i_node)
    #########
    def unrank_leaf(r, M, absolute_pairs):
        shift = 0
        while True:
            i, j = rank_to_pair(r+shift, M)
            lo, hi = rank2node(i), rank2node(j)
            I, J = node2absrank(lo), node2absrank(hi)
            R = bisect.bisect(absolute_pairs, (I, J))
            if shift == R:
                break
            shift =  R
        return (-lo, -hi)
    #########
    # add constants
    k, _ = L[0] # index of the root
    """
    we store -i_node in pool and (-lo, -hi) in pairs
    (used for algorithm bisect of python which need data to be sorted in INCREASING order, a bit lame)
    """
    # pairs= [[] for _ in range(k+1)]
    abs_pairs= [[] for _ in range(k+1)]
    pool = [[] for _ in range(k+1)]
    i_True = len(L)
    L.append((0,(None, None))) # True
    i_False = len(L)
    L.append((0,(None, None))) # False
    pool[0] = [-i_False, -i_True]
    S = []
    i_node = 0 # root since first node in postorder
    while len(S) > 0 or not is_terminal(i_node):
        if not is_terminal(i_node):
            S.append(i_node)
            _, delta = L[i_node]
            lo, _ = delta
            i_node = lo
        else:
            i_node = S.pop()
            level, delta =  L[i_node]
            lo, hi = delta
            pool[level].append(-i_node) # current node is not new and belongs to the pool
            if lo == None:
                r = D.pop()
                #assert(x== i_node)
                if hi == None:
                    M =sum([len(pool[i]) for i in range(level)])
                    lo, hi = unrank_leaf(r, M, abs_pairs[level])
                else:
                    lo = -rank2node(r)
            elif hi == None:
                r = D.pop()
                # assert(x== i_node)
                i = node2rank(lo) # index of lo within the pool
                shift = 1 if r >= i else 0
                hi = -rank2node(r+shift) # hi cannot point to lo kid
            delta = lo, hi
            L[i_node] = level, delta
            i_node = hi
            I, J = node2absrank(abs(lo)), node2absrank(abs(hi))
            assert((I, J) not in abs_pairs[level])
            bisect.insort(abs_pairs[level], (I, J))  # we store pairs of aboslute indices (order is breadth first left to right)
##########################################################
def unrank_bdd_from_profile(r, p):
    L, D = unrank_bdd_from_profile_incomplete(r, p)
    adjust(L, D, p)# all nodes have been distributed, we can unrank edges correctly
    return L
##########################################################
def unrank_bdd_from_profile_incomplete(r, p):
    L = []
    k = len(p)-1
    S = []
    D = []
    E = (k,)
    while len(E) > 0:
        l = E[-1]
        E = E[:-1]
        while p[l] == 0:
            l -= 1
        if l > 0:
            q = tuple(p[:l])+(p[l]-1,)+p[l+1:]
            Delta = count_anchored_extensions(q, E + (l-1, l-1)) - count_anchored_extensions(q, E + (l-1,)) - (p[l] - 1)* count_anchored_extensions(q, E)
            if r - Delta < 0: #create a node at level next_level
                i_node = len(L) # index of node to be created
                L.append((l, (None, None)))
                if len(S) > 0:
                    i_parent, kid = S.pop()
                    level, delta = L[i_parent]
                    L[i_parent] = level, (delta[0], i_node) if kid == 1 else (i_node, delta[1])
                    if kid == 1:
                        S.append((i_parent, 0))
                S.append((i_node, 1))
                p = q
                E = E + (l-1, l-1)
                continue
            r -= Delta
            E = E + (l-1,)
        else:   # no node can be created: setup pointer
            i_node, kid = S.pop()
            level, _ = L[i_node]
            w = sum(p[:level])
            if kid == 1:
                nb_choices = w*(w-1) - p[level]
                Delta = nb_choices * count_anchored_extensions(p, E[:-1]) # how many considering the current node is a leaf?
                if r - Delta < 0:
                    D.append(r % nb_choices)
                    r = r // nb_choices
                    E = E[:-1]
                    continue
                r -= Delta # not a leaf
                nb_choices = w-1 # hi cannot point to lo
                S.append((i_node, 0)) # push back node to compute lo kid later
            else:
                nb_choices = w
            D.append(r % nb_choices)
            r //= nb_choices
    return L, D
###########################
def complete_profile(p, n, k):
    P = (0,1)
    for x in p:
        # if x > 0: # does not happen often, not worth the test
        P = iter(x, P)
    bottom = n - sum(p)
    #print("p ={} P={} n={} k={} bottom={}".format(p, P, n, k, bottom))
    r = 0
    for d, c in enumerate(P):
        if c != 0:
            count = bdd.expand_monomial(n, d, k-len(p))
            #print("c={} d={}\t{}".format(c, d, count))
            if len(count) > bottom:
                r += c*count[bottom]
    return r
##########################
def unrank_profile(rank,  n, nb_vars):
    r = rank
    p =()
    for k in range(nb_vars):
        top = sum(p)
        bottom = n - top
        for m in range(1 +min(top+1, bottom-ceil(sqrt(bottom)))):
            Delta = complete_profile(p+(m,), n, nb_vars)
            if r - Delta  < 0 :
                p = p +(m,)
                break
            r -= Delta
    return r, p
########################
########################################
def calc(p, m):
    """
    p is an incomplete profile
    a is the profile of the multiset entries
    return: a polynomial in Z[X]
    caution: profile p is incomplete (without the last value 2 for constants
    """
    P = (1,)
    for x, y in zip_longest(p, m, fillvalue=0): # padding if needed
        P = iter(x, (0,)*y + P)
    return P
########################################
#@functools.lru_cache(maxsize=None)
def count_multi_BDDs(p, m):
    """
    p is a incomplete profile (without the constants level)
    """
    s = calc(p, m)
    return eval(s)
###########################
#@functools.lru_cache(maxsize=None)
def count_anchored_extensions(p, a):
    """
    in module bdd: profiles are (p_k, ..., p_1)
    Here we use (p_0=2, p_1, ..., p_k) and use an anchor list instead of anchor profile
    this is a "wrapper" around count_multi_BDDs for unranking
    """
    p, a  = remove_empty_levels(p, a)
    m = multi_profile(a)
    res =  count_multi_BDDs(tuple(p[:0:-1])+(0,),m[::-1])
    return res
###############################
#@functools.lru_cache(maxsize=None)
def eval(P):
    i = 1
    s = 0
    for x in P:
        s += i*x
        i *= 2
    return s
###############################
#@functools.lru_cache(maxsize=None)
def iter(x, P):
    """
    x integer (number of nodes)
    P polynomial in Z[x]
    """
    Q = tuple()
    for k in range(x, len(P)):
        coeff = P[k]
        #        if coeff != 0: # does not happen often, not worth the test
        Q = add(Q,  tuple(coeff*x for x in bdd._R(x, k)))
    return Q
################################
def unrank_bdd(rank, size, nb_vars, essential=False):
    """
    return the BDD with rank among BDDs of size with nb_vars
    """
    r = rank
    start = timer()
    r, p = unrank_profile(rank, size, nb_vars)
    end = timer()
    print("#time unranking profile: {} s".format(end - start), flush=True)
    p = (2,)+p[::-1] # translate to have (p_0=2, p_1, ..., p_k)
    start = timer()
    L =  unrank_bdd_from_profile(r, p)
    end = timer()
    print("#time unranking BDD from profile: {} s".format(end - start), flush=True)
    return L
#####################################
if __name__ == '__main__':
    max_level = None
    nb_vars = int(sys.argv[1])
    n = int(sys.argv[2])
    if len(sys.argv) == 3:
        r = None
    elif len(sys.argv) == 4:
        r = int(sys.argv[3])
    elif len(sys.argv) == 5:
        r = int(sys.argv[3])
        max_level = int(sys.argv[4])
    else:
        print("usage: python3 unrank_bdd nb_vars size [rank [max_level]]\n# rank optional, max_level used only for dot output")
    if n > bdd.max_size(nb_vars):
        print(" size too big !!! setting to {}".format( bdd.max_size(nb_vars)))
        n = bdd.max_size(nb_vars)
    print("#maximum size ({} vars): {}".format(nb_vars, bdd.max_size(nb_vars)))
    bdd.init(n, nb_vars)
    start = timer()
    total = bdd.expand_monomial(n, 1, nb_vars)
    total = sum(total)
    print("# precomputation step: {} s".format(timer()-start))
    print("#BDDs = {}, size {} with {} vars".format(total, n, nb_vars))
    if total == 0:
        print("No BDD!")
        exit(-1)
    if r is None:
        r = random.randint(0, total-1) # 0<= r <= total-1
    print("#********* rank = {} **************".format(r))
    print("# r={} / {}".format(r, total))
    assert(r<total)
    L = unrank_bdd(r, n, nb_vars)
    if dot:
        print("{}".format(to_dot(L, nb_vars, max_level)))
    print("# ROBDD: {}".format(L))
    print("# profile={}".format(getProfile(L)))
    #print("# truth table=\"{}\"".format(to_TT(L)))
    print("# detailed weight={}".format(weight(L)))
#    print("# cache count_anchored_extensions: {}".format(count_anchored_extensions.cache_info()))
    
    print("# cache expand_monomial: {}".format(bdd.expand_monomial.cache_info()))
    print("# cache R: {}".format(bdd._R.cache_info()))
