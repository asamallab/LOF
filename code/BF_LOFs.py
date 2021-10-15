"""
====================================================================
Program to obtain perform operations on Boolean functions and get their
various properties
Written by: Ajay Subbaroyan
Reference: 
====================================================================
"""
from BF_generators import *
from pyeda.inter import *
from itertools import combinations
from string import ascii_lowercase

def is_LOF (signs, f, kind):
    '''
    signs: signs of the inputs as a string i.e 'aaia' where
    'a' represents activator and 'i' represents inhibitor
    f: BF as a binary string ; the string from left to right
    is the output of the truth table from top to bottom
    kind: 'AND_NOT' or 'OR_NOT' or 'AND-pairs' or 'OR-pairs'

    returns True or False

    #instance
    >>> is_LOF('aia','00001011', 'OR-NOT')
    >>> False
    
    '''
    k, bias, bias_z = len(signs), f.count('1'), f.count('0')
    act, inh = [], []
    for index, sign in enumerate(signs):
        if sign == 'a':
            act += [k-index]
        else:
            inh += [k-index]
    I = bf(k,f).indices()
    zeros, ones = I[0], I[1]

    if len(act) >= 1 and len(inh) >=1 :
        
        #Check AND_NOT
        if kind == 'AND-NOT':
            if bias == 2**(len(act)) - 1:
                inh_zeros = [zeros[ele] for ele in inh]
                ones_output_index = list(set.intersection(*map(set, inh_zeros)))
                ones_output = itemgetter(*ones_output_index)(f)
                if (ones_output[0] == '0' and ones_output.count('1') == bias):
                    return True
                else:
                    return False
            else:
                return False

        #Check OR_NOT
        elif kind == 'OR-NOT':
            if bias_z == 2**(len(inh)) - 1:        
                act_zeros = [zeros[ele] for ele in act]
                zeros_output_index = list(set.intersection(*map(set, act_zeros)))
                zeros_output = itemgetter(*zeros_output_index)(f)
                if (zeros_output[0] == '1' and zeros_output.count('0') == bias_z):
                    return True
                else:
                    return False
            else:
                return False
        
        #Check AND-pairs function
        elif kind == 'AND-pairs':
            Z = [zeros[a] for a in act]
            O = [ones[i] for i in inh]
            rows = list(set(list(set(Z[0]).intersection(*Z[1:])) + list(set(O[0]).intersection(*O[1:]))))
            if len(rows) == bias_z:
                for row in rows:
                    if f[row] == '1':
                        return False
                return True
            else:
                return False

        #Check OR-pairs function
        elif kind == 'OR-pairs':
            Z = [zeros[i] for i in inh]
            O = [ones[a] for a in act]
            rows = list(set(list(set(Z[0]).intersection(*Z[1:])) + list(set(O[0]).intersection(*O[1:]))))
            if len(rows) == bias:
                for row in rows:
                    if f[row] == '0':
                        return False
                return True
            else:
                return False

    else:
        return ("LOF is not possible")

def generate_link_func (sign, kind):
    '''
    signs: signs of the inputs as a string i.e 'aaia' where
    'a' represents activator and 'i' represents inhibitor
    kind: 'AND_NOT' or 'OR_NOT' or 'AND-pairs' or 'OR-pairs'

    returns the LOF of the given type corresponding to the
    given sign combination

    #instance
    >>> generate_link_func ('aaia', 'AND-pairs')
    >>> '0100110011001100'
  
    '''
    var = ascii_lowercase
    exp = '('
    m, n = sign.count('a'), sign.count('i')
    
    if m>=1 and n>=1:
        if kind == 'AND-NOT':
            for i in range(m+n):
                if i < m-1:
                    exp += var[i]+'|'
                elif i == m-1:
                    exp += var[i]+')&~('
                elif i < m+n-1:
                    exp += var[i]+'|'
                elif i == m+n-1:
                    exp += var[i]+')'
            
        elif kind == 'OR-NOT':
            for i in range(m+n):
                if i < m-1:
                    exp += var[i]+'|'
                elif i == m-1:
                    exp += var[i]+')|~('
                elif i < m+n-1:
                    exp += var[i]+'|'
                elif i == m+n-1:
                    exp += var[i]+')'

        elif kind == 'AND-pairs':
            for i in range(m+n):
                if i < m-1:
                    exp += var[i]+'|'
                elif i == m-1:
                    exp += var[i]+')&('
                elif i < m+n-1:
                    exp += '~'+var[i]+'|'
                elif i == m+n-1:
                    exp += '~'+var[i]+')'

        elif kind == 'OR-pairs':
            for i in range(m+n):
                if i < m-1:
                    exp += var[i]+'&'
                elif i == m-1:
                    exp += var[i]+')|('
                elif i < m+n-1:
                    exp += '~'+var[i]+'&'
                elif i == m+n-1:
                    exp += '~'+var[i]+')'
    else:
        return ("Enter a valid 'sign' for generating LOFs")

    tt = str(expr2truthtable(expr(exp)))
    index = [ind-1 for ind, char in enumerate (tt) if char == '\n']
    f = list(itemgetter(*index)(tt))[1:]
    bf_str = ''.join(f)
    perms = bf(m+n, bf_str).all_perms()

    for term in perms:
        if check_if(m+n, term).is_UF() == sign:
            return term
