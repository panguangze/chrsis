#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    tenxtools.bash
    ~~~~~~~~~~~~~~

    bash tools to access directory

    @Copyright: (c) 2018-07 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""

import os


def get_files(in_dir, fn_suffix=''):

    in_fns = []
    for root, _, fns in os.walk(in_dir):
        if not fns:
            continue
        for fn in fns:
            if not fn.endswith(fn_suffix):
                continue
            in_fns.append(os.path.join(root, fn))

    return in_fns
