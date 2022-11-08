# -*- coding: utf-8 -*-
"""
    simulate.math
    ~~~~~~~~~~~~~

    self-defined math function

    @Copyright: (c) 2017 by Lingxi Chen (chanlingxi@gmail.com). 
    @License: LICENSE_NAME, see LICENSE for more details.
"""
from numpy import random

def randint(low, high):
    if low == high: return low
    else:           return random.randint(low, high)
