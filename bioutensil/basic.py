# -*- coding: utf-8 -*-
"""
    tenxtools.basic
    ~~~~~~~~~~~~~~~

    basic operation on linux

    @Copyright: (c) 2018-04 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""

import gzip


def safe_open(fn, mod):
    if fn.endswith('gz'):
        return gzip.open(fn, mod)
    else:
        return open(fn, mod)
