#import numpy as numpy
from array import array
from Utils import *

#################################################
# Class Spine
class Spine:
    def __init__(self, target_profile):
        self.nb_vars = len(target_profile)
        self.base  = [x+1 for x in max_profile(self.nb_vars)]
        self.target_profile = tuple([2] + list(target_profile))
        self.target_size = sum(self.target_profile)
        self.height = len(self.target_profile)
        self.node = array("i", [0]*(4*self.target_size)) # a node = (local_rank, level, lo, hi)
        self.table = array("i", [0]*self.target_size) # list of nodes  by level and creation order
        self.rank = array("i", [0]*self.target_size) # sorted list according to pairs of positions of children by level
        self.rank_offset =  array("i", [0]*self.height) # used to compute rank
        self.offset = array("i", [0]*(self.height+1)) # starting index of levels in table
        start = 0
        i = 0
        for x in self.target_profile:
            self.offset[i] = start
            i += 1
            start += x
        self.offset[i] = start # sentinel
        self.size = 0
        self.profile = array("i", [0]*len(self.target_profile)) # initial profile
        self.add_node(0, 0, 0) # lo=0, hi=0, level=0: False
        self.add_node(1, 1, 0) # lo=1, hi=1, level=0: True
        
    def get_spine_profile(self):
        return tuple(self.profile[1:])

    def decode(self, r):
        return [2]+ unrank(r, self.base)

    def encode(self):
        return rank(self.profile[1:], self.base)

    def iter_level(self, level):
        for i in range(self.offset[level], self.offset[level]+self.profile[level]):
            i_node= self.table[i]
            yield i_node
    def print(self):
        u = "size={}\nprofile={}\nnode={}\noffset={}\ntable={}\nrank={}\n".format(self.size,self.profile, self.node, self.offset,self.table, self.rank)
        for l in range(self.height):
            u+= "level {}:\n".format(l)
            for i_node in self.iter_level(l):
                local_rank, level, lo, hi = self.get_node(i_node)  
                u+="   #{}=({}, {})\n".format(i_node, lo, hi)
        return u
    def __repr__(self):
        u = "size={}\nprofile={}\nnode={}\noffset={}\ntable={}\nrank={}\n".format(self.size,self.profile, self.node, self.offset,self.table, self.rank)
        for l in range(self.height):
            u+= "level {}:\n".format(l)
            for i_node in self.iter_level(l):
                local_rank, level, lo, hi = self.get_node(i_node)  
                u+="   #{}=({}, {})\n".format(i_node, lo, hi)
        return u

    """ 
    accessors to nodes 
    """
    def get_node(self, i):
        """
        i: index of creation
        returns (local_rank, level, lo, hi)
        """
        return tuple([self.node[4*i+j] for j in range(4)])
    
    def get_nth(self, rank, level=None):
        """ order is rank bottom-up, left-to-right 
        return (rank_index, node_index), level is optional 
        """
        r = rank
        if level is None:
            level = 0
            while r >= self.profile[level]:
                r -= self.profile[level]
                level += 1
        i_node = self.table[self.offset[level]+r]
        return i_node

    
    def set_node(self, i, local_rank, level, lo, hi):
        """
        set hi and lo for node of index i
        """
        self.node[4*i] = local_rank
        self.node[4*i+1] = level
        self.node[4*i+2] = lo
        self.node[4*i+3] = hi
        
    def find_rank(self, lo, hi, level):
        i, j  = self.get_rank(lo), self.get_rank(hi)
        pos = 0
        while pos < self.profile[level]:
            i_node = self.rank[self.offset[level]+pos]
            l_rank, l_level, l_lo, l_hi = self.get_node(i_node)
            k, l = self.get_rank(l_lo), self.get_rank(l_hi) 
            # if cmp_pair((k, l), (i, j)) == 0:
            #     print("!!! find_rank: (lo, hi)=({}, {}) level={}".format(lo, hi, level)) # equality should not be possible  
            if cmp_pair((k, l), (i, j)) < 0: # equality should not be possible  
                pos += 1
            else:
                break
        return pos
            
    def insert_rank(self, i_node, level):
        """
        inserts a node indexed by i_node
        in self.node a in ordered list of nodes
        (order given by struct x=(lo, hi) for x in L)
        return (rank of r, new list)
        """
        l_rank, l_level, lo, hi = self.get_node(i_node)
        i_rank = self.find_rank(lo, hi, level)
        for i in range(self.profile[level], i_rank, -1): #beware profile before insertion
            self.rank[self.offset[level]+i] = self.rank[self.offset[level]+i-1]
        self.rank[self.offset[level]+i_rank] = i_node

    def insert_table(self, i_node, level):
        """
        inserts a node indexed by i_node at the end the list at level
        """
        self.table[self.offset[level]+self.profile[level]] = i_node

    def add_node(self, lo, hi, level):
        """
        add a node at level with children (lo,hi)
        returns (rank_index, node_index)
        """
        i_node = self.size
        local_rank = self.profile[level]
        self.set_node(i_node, local_rank, level, lo, hi)
        self.insert_table(i_node, level)
        self.insert_rank(i_node, level)
        self.profile[level] += 1
        for l in range(level+1, self.height):
            self.rank_offset[l] += 1
        self.size += 1
        return i_node

    def to_list(self):
        return self.to_list_aux(self.get_nth(0, self.height-1), [0, 1])[1]

    def to_list_aux(self, root, visited):
        # labeling of variables according to level
        # recursive traversing
        if root in visited:
            return visited, [root]
        local_rank, level, lo, hi = self.get_node(root)
        new_visited, left = self.to_list_aux(lo, visited)
        new_visited, right = self.to_list_aux(hi, new_visited)
        new_visited.append(root)
        return new_visited, [root, left, right, level] 

    def get_rank(self, i_node):
        l_rank, l_level, lo, hi = self.get_node(i_node)
        return l_rank+self.rank_offset[l_level]
    def unrank_singleton(self, rank):
        return self.get_nth(rank)
    
    def unrank_pair(self, rank, level):
        M = sum(self.profile[:level])
        shift = 0
        previous_shift = 0
        r = rank
        while True:
            i, j = rank_to_pair(rank+shift, M)
 #           print("unrank rank={} shift={} M={} pair={}".format(rank+shift, shift, M, (i, j)))
            while shift <  self.profile[level]:
                i_node = self.rank[self.offset[level]+shift]
                local_rank, level, lo, hi = self.get_node(i_node)
                (k, l) = (self.get_rank(lo), self.get_rank(hi))
#                print("cmp with {}".format((k,l)))
                if cmp_pair((i, j), (k, l)) >= 0:
                    shift += 1
                else:
                    break
            if shift > previous_shift:
                previous_shift = shift
            else:
                break
#        print("final unrank_pair initial rank={} M={} pair={}".format(rank,M, (i, j)))
        return (self.get_nth(i), self.get_nth(j))
