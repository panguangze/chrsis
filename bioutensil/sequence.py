# -*- coding: utf-8 -*-
"""
    simulate.sequence
    ~~~~~~~~~~~~~~~~~

    DNA sequence manipulation. 

    @Copyright: (c) 2017 by Lingxi Chen (chanlingxi@gmail.com). 
    @License: LICENSE_NAME, see LICENSE for more details.
"""
complements = {
    'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
    'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n',
    }
def reverse_complement(seq):
    return ''.join([complements[x] for x in seq[::-1]])
