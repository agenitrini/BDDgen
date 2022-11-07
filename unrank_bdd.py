#-*- coding: utf-8 -*-
import math
import sys as sys
import time
import random as random
import bisect as bisect
import bdd_func as bdd
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
def truth_table(L):
    """
    return the truth table of the BDD L; node : n-2 -> True and n-1 -> False
    """
    T = dict()
    i = 0
    for a in L[:-2]:
        T[i] = (a[0], abs(a[1][0]), abs(a[1][1]))
        i += 1
    i -= 1
    T[i+1] = (0,)
    T[i+2] = (0,)

    TT = dict()
    TT[i+1] = [True]
    TT[i+2] = [False]
    while i>=0:
        TT[i] = TT[T[i][1]]*(2**(T[i][0] - T[T[i][1]][0] - 1)) + TT[T[i][2]]*(2**(T[i][0] - T[T[i][2]][0] - 1))
        if len(TT[i]) != 2**(T[i][0]):
            print('pb ', i, L)
        i -= 1
    return TT[0]
###########################
def getProfile(L):
    x = []
    for l, _ in L:
        x.append(l)
    p = [0 for _ in range(max(x)+1)]
    for i in x:
        p[i] += 1
    return tuple(p[::-1])
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
###########################
def count_multi_BDDs(p, a):
    """
    in module bdd profiles are (p_k, ..., p_1)
    Here we use (p_0=2, p_1, ..., p_k)
    """
    p, a  = remove_empty_levels(p, a)
    #print(p, a)
    m = multi_profile(a)
    res =  bdd.count_multi_BDDs(tuple(p[:0:-1])+(0,),m[::-1])
    return res
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

def adjust(L, D, p):
    """
    Adjust relative pointers to the pool
    *once the spine is complete*
    by a infix traversal of the tree (so the pool is known at each step)
    profile in order (p_0=2, ..., p_k=1)
    """
    ### HELPER FUNCTIONS: current pool and final profile are known
    def is_terminal(n):
        return n is None or n < 0 or n >= len(L)-2
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
        #print("#i_node={} pool={} rank={}".format(i_node, pool[level], bisect.bisect_left(pool[level], -i_node)))

        return r+bisect.bisect_left(pool[level], -i_node)
    #########
    def unrank_leaf(r, M,  P):
        rel = []
        for lo, hi in P:
            lo, hi = -lo, -hi
            i, j = node2rank(lo), node2rank(hi)
            R = pair_to_rank(i, j, M)
            rel.append(R)
        rel.sort()
        shift = 0
    
        while True:
          R = bisect.bisect(rel, r+shift) # R = #values(rel) <= r
          if shift == R:
              break
          shift =  R 
        i, j = rank_to_pair(r+shift, M)
        lo, hi  = rank2node(i), rank2node(j)
        return (-lo, -hi)
    #########
    # add constants
    k, _ = L[0] # index of the root
    """
    we store -i_node in pool and (-lo, -hi) in pairs
    (used for algorithm bisect of python which need data to be sorted in INCREASING order, a bit lame)
    """
    pairs= [[] for _ in range(k+1)]
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
                    lo, hi = unrank_leaf(r, M, pairs[level])
                    assert((lo, hi) not in pairs[level])
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
            pairs[level].append((-abs(delta[0]), -abs(delta[1])))  # we store pairs (-lo, -hi) in P
##########################################################
def unrank_bdd_from_profile(r, p):
    L, D = unrank_bdd_from_profile_incomplete(r, p)
    adjust(L, D, p)# all nodes have been distributed, we can unrank edges correctly
    return L

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
            Delta = count_multi_BDDs(q, E + (l-1, l-1)) - count_multi_BDDs(q, E + (l-1,)) - (p[l] - 1)* count_multi_BDDs(q, E)
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
                Delta = nb_choices * count_multi_BDDs(p, E[:-1]) # how many considering the current node is a leaf?
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
def unrank_bdd(rank, size, nb_vars, essential=False):
    """
    return the BDD with rank among BDDs of size with nb_vars
    """
    r = rank
    r, p = bdd.unrank_profile(rank, size, nb_vars)
    p = (2,)+p[::-1] # translate to have (p_0=2, p_1, ..., p_k)
    L =  unrank_bdd_from_profile(r, p)
    return L
#####################################
if __name__ == '__main__':
    max_level = None
    if len(sys.argv) == 3:
        n = int(sys.argv[1])
        nb_vars = int(sys.argv[2])
        r = None
    elif len(sys.argv) == 4:
        r = int(sys.argv[3])
        n = int(sys.argv[1])
        nb_vars = int(sys.argv[2])
    elif len(sys.argv) == 5:
        r = int(sys.argv[3])
        n = int(sys.argv[1])
        nb_vars = int(sys.argv[2])
        max_level = int(sys.argv[4])
        
    total = bdd.count(n, nb_vars)
    print("#maximum size ({} vars): {}".format(nb_vars, bdd.max_size(nb_vars)))
    print("#BDDs = {}, size {} with {} vars".format(total, n, nb_vars))
    if total == 0:
        print("No BDD!")
        exit(-1)
    if r is None:
        r = random.randint(0, total-1) # 0<= r <= total-1
    print("#********* rank = {} **************".format(r))
    print("# r={} / {}".format(r, total))
    L = unrank_bdd(r, n, nb_vars)
    print("# ROBDD: {}".format(L))
    print("# profile={}".format(getProfile(L)))
    #    print("# pretty tree representation=\n"+pt.drawTree2(True)(True)(ns.get_tree(ns.root)))
    #    print("# weight={}".format(ns.weight()))
    print("# cache count_multi_BDDs: {}".format(bdd.count_multi_BDDs.cache_info()))
    print("# cache R: {}".format(bdd.R.cache_info()))
    print("{}".format(to_dot(L, nb_vars, max_level)))
