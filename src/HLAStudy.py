#-*- coding: utf-8 -*-

import os, sys, re
from os.path import exists, dirname, basename, join
import pandas as pd

import src.MHCpredError as MHCpredError
import src.bashPLINK as bashPLINK
from src.MakeSS import MakeSS
from src.Study import Study

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))



class HLA_Study(Study):
    """

    - refineBP.py
    - GTtrick.py
    - AA | HLA | SNPS | ...

    """

    def __init__(self, _input, _out, _label="", _a1_allele=None,
                 _pheno=None, _pheno_name=None, _covar=None, _covar_name=None, _pcs=None):
        # super().__init__(_input, _out, _a1_allele, _pheno, _pheno_name, _covar, _covar_name, _pcs)
        Study.__init__(_input, _out, _a1_allele, _pheno, _pheno_name, _covar, _covar_name, _pcs)

        ### Main variables
        self.label = _label if bool(_label) else basename(_input)  # Override

    # def refine