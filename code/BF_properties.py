"""
====================================================================
Program to obtain perform operations on Boolean functions and get their
various properties
Written by: Ajay Subbaroyan
Reference: 
====================================================================
"""

import numpy as np
from string import ascii_lowercase
from itertools import chain, permutations
from operator import itemgetter

##Abbreviations
##BF: Boolean Function, EF: Effective Functions, CF: Canalyzing functions
##UF: Unate Function, NCF: Nested Canalyzing Functions

class bf:
    '''
    #functionality
    This class returns the various properties of Boolean functions.
    '''
    def __init__(self, k, f):
        '''
        #arguments
        k: number of inputs to a BF
        f: BF as a binary string : the string from left to right
        is the output of the truth table from top to bottom
        ''' 
        self.f = f
        self.k = k
        self.chars = ascii_lowercase[:self.k]
        assert np.log2(len(f)) == self.k, "Number of inputs do not match size of the truth table"

    def indices(self):
        '''
        returns two dictionaries of the indices of the zeros
        and ones respectively for every input
        
        Note:'1' is the input column closest to the output column
        of the truthtable and k the farthest

        #instance
        >>> bf(3, '11001101').indices()
        >>> {0: {1: [0, 2, 4, 6], 2: [0, 1, 4, 5], 3: [0, 1, 2, 3]}, 1: {1: [1, 3, 5, 7], 2: [2, 3, 6, 7], 3: [4, 5, 6, 7]}}
        
        '''
        z = {}
        o = {}
        L = [i for i in range(2**self.k)]
        temp = []
        for num in range(self.k):
            for i in range(0,(2**self.k),2**(num+1)):
                temp += L[i:i+2**num]
            z[num+1] = temp
            o[num+1] = list(set(L)-set(temp))
            temp = []
        return {0:z, 1:o}

    def inps_neibs(self):
        '''
        returns a list all neighbours 1 Hamming distance away,
        for each vertex of the hypercube

        #instance
        >>> bf(3, '10000101').inps_neibs()
        >>> {0: [4, 2, 1], 1: [5, 3, 0], 2: [6, 0, 3], 3: [7, 1, 2], 4: [0, 6, 5], 5: [1, 7, 4], 6: [2, 4, 7], 7: [3, 5, 6]}
        '''
        n = 2**self.k
        Dct = {}
        for line in range(n):
            b = bin(line)[2:].zfill(self.k)
            x = []
            for i in range(len(b)):
                if i == 0 :
                    x += [int(bin(not int(b[i]))[-1] + b[1:], 2)]
                elif i == len(b)-1:
                    x += [int (b[:-1] + bin(not int(b[i]))[-1], 2)]
                else:
                    x += [int (b[:i] + bin(not int(b[i]))[-1] + b[i+1:], 2)]
            Dct[int(b,2)] = x
        return Dct

    def avg_sensitivity(self):
        '''
        returns the average sensitivity of the BF (between 1 and k)

        #instance
        >>> bf(3, '10000101').avg_sensitivity()
        >>> 1.75
        '''
        I = bf.inps_neibs(self)
        tot = 0
        for pos in I:
            for neib in I[pos]:
                if self.f[neib] != self.f[pos]:
                    tot += 1
        return tot/(2**self.k)

    def bias(self, norm=False):
        '''
        returns the Hamming weight or bias of the BF (if norm = False),
        else it returns the unnormalized bias

        #instance
        >>> bf(3, '10000101').bias()
        >>> 3
        '''
        if norm == False:
            return self.f.count('1')
        else:
            return self.f.count('1')/(2**self.k)

    def is_cana_in_input (self,i):
        '''
        returns:
        '0' if the canalyzing input for input 'i' is '0'
        '1' if the canalyzing input for input 'i' is '1'
        '01' if the canalyzing input for input 'i' is both '0' and '1'
        False if input is not canalyzing
        
        #arguments
        i: Integer from 1 to k

        #instance
        >>> bf(3, '11100000').is_cana_in_input(3)
        >>> '1'
        '''
        self.z = bf.indices(self)[0]   #indices of zeros 
        self.o = bf.indices(self)[1]   #indices of ones

        z_ele = list(itemgetter(*self.z[i])(self.f))
        o_ele = list(itemgetter(*self.o[i])(self.f))
        
        is_c_z = all(bit == z_ele[0] for bit in z_ele) 
        is_c_o = all(bit == o_ele[0] for bit in o_ele)
        
        if is_c_z and is_c_o:
            return '01'
        elif is_c_z:
            return '0'
        elif is_c_o:
            return '1'
        else:
            return False

    def cana_depth (self):
        '''
        returns the canalyzing depth of the BF (an integer from 0 to k)

        #instance
        >>> bf(3, '11001010').cana_depth()
        >>> 0
        '''
        if self.k == 1:
            if self.f =='01' or self.f == '10':
                return 1
            else:
                return 0

        if self.f.count('0') == 2**self.k or self.f.count('1') == 2**self.k:
            return 0
        
        else:
            for i in range(1,self.k+1):
                if bf.is_cana_in_input(self, i) == '0': # canalyzing input is '0'
                    self.f = ''.join(itemgetter(*self.o[i])(self.f))
                    self.k -= 1
                    return bf.cana_depth(self) + 1
                
                elif bf.is_cana_in_input(self, i) == '1': #canalyzing input is '1'
                    self.f = ''.join(itemgetter(*self.z[i])(self.f))
                    self.k -= 1
                    return bf.cana_depth(self) + 1
                
                elif bf.is_cana_in_input(self, i) == '01': #canalyzing inputs are '0' and '1'
                                                           #This implies that the remaining
                                                           #inputs are non-canalyzing
                    return 1
        return 0        

    def right_shift (self):
        '''
        returns a single 'right' cyclic permutation of the  truth table
        i.e  3 | 2 | 1 | output ---> 1 | 3 | 2 | output

        #instance
        >>> bf(3, '11001010').right_shift()
        >>> '11100100'
        '''
        m = int(len(self.f)/2)
        upper_half, lower_half = self.f[:m], self.f[m:]
        perm = ''
        for i in range(m):
            perm += upper_half[i] + lower_half[i]
        return perm

    def swap_rows (self,cols):
        '''
        returns a BF with certain negated literals

        #arguments
        cols: the set of inputs to be negated as a list

        #instance
        >>> bf(3, '11100000').swap_rows([1,3])
        >>> '00001101'
        '''
        I = bf.indices(self)
        str_lst = list(self.f)    
        for col in cols:
            I0, I1 = I[0][col], I[1][col]  
            for i in range(len(I0)):
                str_lst[I0[i]], str_lst[I1[i]] =  str_lst[I1[i]], str_lst[I0[i]] 
        return ''.join(str_lst)

    def some_perms(self, n):
        '''
        returns the permutations of the BF for a given 'n'

        #arguments
        n: number of inputs to be permuted

        #instance
        >>> bf(3, '11100000').some_perms(2)
        >>> ['11100000']
        '''
        if n == 1:
            return ''
        
        else:
            parts = [self.f[i:i+2**n] for i in range(0,2**self.k,2**n)]
            perms_k = [self.f]
            for times in range(1,n):
                perms_k += [''.join([bf(int(np.log2(len(part))),part).right_shift() for part in parts])]
                self.f = perms_k[-1]
                parts = [self.f[i:i+2**n] for i in range(0,2**self.k,2**n)]
            Q = list(set(perms_k))
            return Q
    
    def all_perms (self):
        '''
        returns all the permutations of the BF (given that all inputs have
        same sign i.e positive or nagetive)

        #instance
        >>> bf(3, '11100000').all_perms()
        >>> ['10101000', '11001000', '11100000']
        '''
        n = self.k
        L_tot = bf(self.k,self.f).some_perms(n)
        L_new = L_tot
        for n in range(self.k-1,1,-1):
            L_new = list(set(chain(*[bf(self.k,f).some_perms(n) for f in L_new])))
            L_tot = list(set(L_tot + L_new))
        return L_tot
