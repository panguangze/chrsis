# -*- coding: utf-8 -*-
"""
    utils.constants
    ~~~~~~~~~~~~~~~

    Useful bio constants specification.

    @Copyright: (c) 2017 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""


chrs = ['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']


hg19_fai_bp = {
    'chr1': 249250621,
    'chr2': 243199373,
    'chr3': 198022430,
    'chr4': 191154276,
    'chr5': 180915260,
    'chr6': 171115067,
    'chr7': 159138663,
    'chr8': 146364022,
    'chr9': 141213431,
    'chr10': 135534747,
    'chr11': 135006516,
    'chr12': 133851895,
    'chr13': 115169878,
    'chr14': 107349540,
    'chr15': 102531392,
    'chr16': 90354753,
    'chr17': 81195210,
    'chr18': 78077248,
    'chr19': 59128983,
    'chr20': 63025520,
    'chr21': 48129895,
    'chr22': 51304566,
    'chrX': 155270560,
    'chrY': 59373566,
    'chrM': 16571,
}

hg19_arm = {
    'chr1': 124535434,
    'chr2': 95326171,
    'chr3': 93504854,
    'chr4': 52660117,
    'chr5': 49405641,
    'chr6': 61830166,
    'chr7': 61054331,
    'chr8': 46838887,
    'chr9': 50367679,
    'chr10': 42254935,
    'chr11': 54644205,
    'chr12': 37856694,
    'chr13': 19000000,
    'chr14': 19000000,
    'chr15': 20000000,
    'chr16': 38335801,
    'chr17': 25263006,
    'chr18': 18460898,
    'chr19': 27681782,
    'chr20': 29369569,
    'chr21': 14288129,
    'chr22': 16000000,
    'chrX': 61632012,
    'chrY': 13104553,
}


def get_arm(chrom, start, end=None):
    if chrom not in hg19_arm:
        return ''
    middle = hg19_arm[chrom]
    start_arm = 'p' if start <= middle else 'q'
    if end:
        end_arm = 'p' if end <= middle else 'q'
    else:
        end_arm = ''

    arm = start_arm + end_arm
    if arm in 'pp' or 'qq':
        arm = arm[0]
    return arm
