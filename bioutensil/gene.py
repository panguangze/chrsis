# -*- coding: utf-8 -*-
"""
    tenxtools.gene
    ~~~~~~~~~~~~~~

    function for gene

    @Copyright: (c) 2018-04 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import os

from bioutensil import basic


def read_gene_list(fn, sep='\t'):
    if not fn:
        return []
    gene_list = []
    for row in basic.safe_open(fn, 'r'):
        if row.startswith('#'):
            continue
        gene = row.strip().split(sep)[0]
        if gene:
            gene_list.append(gene)

    return gene_list

def read_gene_list_from_dir(in_dir):
    gene_set = set()
    for gs in read_gene_dict_from_dir(in_dir).values():
        gene_set = set(gs) | gene_set
    return list(gene_set)

def read_gene_dict_from_dir(in_dir):
    if not in_dir:
        return {}
    if os.path.isfile(in_dir):
        _, fn = os.path.split(in_dir)
        return {fn: read_gene_list(fn)}
    gene_dict = {}
    for fn in os.listdir(in_dir):
        in_fn = os.path.join(in_dir, fn)
        gene_dict[fn] = read_gene_list(in_fn)

    return gene_dict
