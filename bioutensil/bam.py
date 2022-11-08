# -*- coding: utf-8 -*-
"""
    utils.bam
    ~~~~~~~~~~~~~~~~~

    Function to process BAM file     

    @Copyright: (c) 2017 by Lingxi Chen (chanlingxi@gmail.com). 
    @License: LICENSE_NAME, see LICENSE for more details.
"""


# check if hex_x is in dec_flag
def in_flag(hex_x, dec_flag):
    # convert hex_x into decimal 
    dec_x = int(hex_x, 16)
    # operation
    dec_res = dec_x & dec_flag

    return dec_res == dec_x

'''
print in_flag('0x1',3)
print in_flag('0x10',81), bin(81)
print in_flag('0x100',81), bin(81)
print in_flag('0x10',82), bin(82)
'''


