# -*- coding: utf-8 -*-
"""
    simulate.fasta
    ~~~~~~~~~~~~~~

    

    @Copyright: (c) 2017-08 by Lingxi Chen (chanlingxi@gmail.com). 
    @License: LICENSE_NAME, see LICENSE for more details.
"""
from Bio import SeqIO


def get_fasta_total_length(fasta_fn):
    length = 0 
    for line in SeqIO.parse(open(fasta_fn, 'rU'), 'fasta'):
        length += len(line.seq)
    return length

   
    
