# -*- cod:  -*-
"""
    utils.sv
    ~~~~~~~~

    Structure variation definition.

    @Copyright: (c) 2017 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""


class Deletion(object):

    """Docstring for Deletion. """

    def __init__(self, l_bp, r_bp, seq='', supports=0):
        self.l_bp = l_bp  # left breakpoint
        self.r_bp = r_bp  # right breakpoint

        self.seq = seq      # seq deleted
        self.supports = supports    # support read count

    def __str__(self):
        pass


class Insertion(object):

    """Docstring for Insertion. """

    def __init__(self, seq, l_bp, r_bp, supports=0):
        self.seq = seq  # seq inserted
        self.l_bp = l_bp    # left breakpoint
        self.r_bp = r_bp    # right breakpoint
        self.supports = supports    # support read count

    def __str__(self):
        pass


class TandemDuplication(object):

    """Docstring for TandemDuplication. """

    def __init__(self, seq, l_bp, r_bp, supports=0):
        self.seq = seq  # seq duplicated
        self.l_bp = l_bp    # left breakpoint
        self.r_bp = r_bp    # right breakpoint
        # self.dup_num = dup_num  # seq duplication times
        self.supports = supports    # support read count

    def __str__(self):
        pass


class Inversion(object):

    """Docstring for Inversion. """

    def __init__(self):
        """TODO: to be defined1. """
        pass


class ITX(object):

    """Docstring for ITX. """

    def __init__(self):
        """TODO: to be defined1. """


class CTX(object):

    """Docstring for CTX. """

    def __init__(self):
        """TODO: to be defined1. """
