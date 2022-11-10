#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    tenxtools.significant_test
    ~~~~~~~~~~~~~~~~~~~~~~~~~~

    @Copyright: (c) 2018-08 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
from scipy import stats
from skidmarks import wald_wolfowitz
from statsmodels.stats.multitest import fdrcorrection

import numpy as np
np.random.seed(12345678)

def fdr(pvals):
    return fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False)

def exponent_distribution(x):
    '''
    null hypothesis: samples from exponent distribution
    '''
    if len(x) <= 2:
        return None, None
    param = stats.expon.fit(x)
    D, p = stats.kstest(x, 'expon', args=param)
    return p < 0.05, p


def multinomial_distribution(x):
    '''
    null hypothesis: samples from multinominal distribution
    '''
    # the expected frequencies are uniform and given by the mean of the observed frequencies.
    chisq, p = stats.chisquare(x)
    return p < 0.05, p


def wald_wolfowitz_test(x):
    '''
    null hypothesis: samples from alternating sequence
    010101010101

    '''
    if len(set(x)) == 1:
        return None, None, None
    if len(x) <= 2:
        return None, None, None
    z, mean, sd, p ,n_runs = wald_wolfowitz(x).values()
    return p < 0.05, p, int(n_runs)/float(len(x))


def rank_sum_test(x1, x2):

    '''
    statistic : float
    The test statistic under the large-sample approximation that the rank sum statistic is normally distributed
    pvalue : float
    The two-sided p-value of the test
    '''
    s, p = stats.ranksums(x1, x2)

    return p < 0.05, p


def t_test(x1, x2, equal_var=True, ind=True):

    if ind:
        s, p = stats.ttest_ind(x1, x2, equal_var=equal_var)
    else:
        s, p = stats.ttest_rel(x1, x2)

    return p < 0.05, p


def anova_test(x1, x2):
    s, p = stats.f_oneway(x1, x2)
    return p < 0.05, p

