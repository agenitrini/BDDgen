import sys as sys
import math as math
######################################################
# Auxiliary functions

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

########################################################################
def cmp_pair(a, b):
    a1,a2 = a
    b1,b2 = b
    r = a1-b1
    if r ==0:
        r = a2-b2
    return r
########################################################################
def pair_to_rank(i, j, M):
    r = i*M+j
    R = r-1- ((r-1)//(M+1))
    return R

########################################################################
def rank_to_pair(r, M):
    R = r+(r//M)+1
    i = R//M
    j = R % M
    return (i, j)

########################################################################
def enumerate_profiles(n, k):
    L = k+1 # len(deb)
    if n > 1 and L == 2:
        return set()
    E = [[2]]
    for i in range(1,L-2):
        F =[]
        for f in E:
            m0 = sum(f)*(sum(f)-1)
            m1 = 2**(L-1-i)
            m2 = (n - sum(f) + 3) // 2
            M = min(m0, m1, m2)+1
            for l in range(M):
                p = f+[l]
                F.append(p)
        E = F
    i = L - 2
    if L > 2:
        G = []
        for f in E:
            l = n-sum(f)+1
            if l <= 2 and l >=0:
                p = f + [l, 1]
                G.append(p)
        E = G
    return E
########################################################################
def count_total_profiles(k):
    i = 3
    nb = 1
    total = 0
    while (nb> 0):
        print("{} {}".format(i, nb))
        i += 1
        P = count_profiles(i, k)
        nb = sum(P.values())
        total += nb
    return total
########################################################################
def count_profiles(size, k):
    """
    size is the size of the BDD
    k  is the index of the root (number of variables)
    """
    n = size - 2 # size of the spine
    L = k+1 # 
    if n > 1 and k == 1:
        return {}
    E = {2:1}
    for i in range(1,k-1):
        F ={}
        for S in E:
            m0 = S*(S-1)
            m1 = 2**(k-i)
            m2 = ((n - S + 1) // 2) + 1
            M = min(m0, m1, m2)+1
            for l in range(M):
                C = S+l
                if C in F:
                    F[C]=F[C]+E[S]
                else:
                    F[C] = E[S]
        E = F
    if L > 2:
        F = {}
        for S in E:
            l = n-S+1
            if l <= 2 and l >=0:
                C = S+l+1
                if C in F:
                    F[C]=F[C]+E[S]
                else:
                    F[C] = E[S]
        E = F
    return E

########################################################################
def rank(t, base):
    n = len(t)
    r = 0
    for i in range(n-1, -1, -1):
        r = r*base[i]+t[i]
    return r

########################################################################
def unrank(r, base):
    i = 0
    t = []
    while (r > 0 and i < len(base)):
        t.append(r % base[i])
        r = r // base[i]
        i += 1
    return t

########################################################################
def max_profile(k):
    L = []
    S= 2
    i = 0
    x = S*(S-1)
    while  S*(S-1) < 2**(k-i-1):  
        L.append(S*(S-1))
        S += S*(S-1)
        i +=1
    for j in range(i, k):
        L.append(2**(k-j-1))
    return L

########################################################################
def max_size(k):
    return sum(max_profile(k))

########################################################################
#####################     
def dot_from_tree(tree):
    global dot
    ### rajouter label noeud
    if len(tree)>2:
        if tree[0] > 1:
            dot += '  ' + str(tree[0]) + ' [label= "X_' + str(tree[-1]) + '"]\n'
    if len(tree) == 4 and tree[1] != []:
        if len(tree[1]) == 4:
            dot += '  ' + str(tree[0]) + ' -> ' +  str(tree[1][0]) + '[style=dotted] \n'
        else:
            dot += '  ' + str(tree[0]) + ' -> ' +  str(tree[1][0]) + '[style=dotted, color=red] \n'
            
        if len(tree[2]) == 4:
            dot += '  ' + str(tree[0]) + ' -> ' +  str(tree[2][0]) + '\n'
        else:
            dot += '  ' + str(tree[0]) + ' -> ' +  str(tree[2][0]) + '[color=red] \n'
        dot_from_tree(tree[1])
        dot_from_tree(tree[2])
  
def level(t,L):
    if len(t) == 4:
        L[t[-1]].append(t[0])
        L = level(t[1],L)
        L = level(t[2],L)
    return L
    
def dot_save(t, n, k, j=''):
    global dot

    dot = 'digraph {\n   graph [ordering=out]; node [shape=circle]; edge []; \n'
    dot +=   'subgraph{'
    dot_from_tree(t)
    dot += '  ' + str(0) + ' [label= "False", shape=rectangle]\n'
    dot += '  ' + str(1) + ' [label= "True", shape=rectangle]\n'
    
    L = [[0,1]]
    for i in range(k):
        L.append([]) 
    L = level(t,L)
    for l in L:
        dot+= '{rank = same;'
        for e in l:
             dot += str(e) + ';'
        dot += '}'
    dot += '}'
    fic = open("tests_dot/BDD_"+j+"_"+str(n)+"_"+str(k)+".dot", "w")
    
    fic.write(dot + '}')
    fic.close()
#
    






if __name__ == '__main__':
    k= int(sys.argv[1])
    M = max_profile(k)
    print("max profile: {}\nMax #nodes ={}".format(M, 2+sum([x for x in M])))
    base = [x+1 for x in M]
    size = sum([math.log(x, 2) for x in base])
    print("size in bits: {}".format(size))
        
    r_max = rank(M, base)
    print("{}".format(r_max))
    print(unrank(r_max,base))
    for i in range(1, k+1):
        L = [0]*(i-1)+[1]
        print("{} -> {}: {}".format(i, L, rank(L, base)))
    r = 191988225
    print(unrank(r,base))

    # test max_size
    # for k in range(30):
    #     L = max_size(k)
    #     print("k={} max size={}".format(k, L))

    # test max_profile
    # for k in range(30):
    #     L = max_profile(k)
    #     print("k={} max profile={}".format(k, L))
    pass
