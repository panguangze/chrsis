# -*- coding: utf-8 -*-
"""
    utils.myprint
    ~~~~~~~~~~~~~

    my pretty print functions

    @Copyright: (c) 2017 by Lingxi Chen (chanlingxi@gmail.com). 
    @License: LICENSE_NAME, see LICENSE for more details.
"""

#######################################################################
#                            Matrix Print                             #
#######################################################################
import numpy as np

def mprint(matrix):
    for row in matrix:
        for item in row:
            item = [ [x for x in line] for line in item.split('\n') ]
            combine = np.hstack((combine, np.array(item)))
            print combine 
            print '----------------'



if __name__ == "__main__":
    m =[[
"""-----
| a b |
|  +  |
-----""",
"""DEL
---"""
            ],
    ]
    mprint(m)

