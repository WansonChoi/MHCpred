#-*- coding: utf-8 -*-

import os, sys, re



class MHCpredError(Exception):

    """
    Base Error for MHCpred.
    """

    pass



class MHCpredInputPreparationError(MHCpredError):

    """
    Error related to preparing input.
    """

    def __init__(self, _message):
        print(_message)
