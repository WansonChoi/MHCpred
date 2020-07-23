#-*- coding: utf-8 -*-

import os, sys, re
from os.path import exists, dirname, basename, join
import pandas as pd

import src.MHCpredError as MHCpredError
import src.bashPLINK as bashPLINK
from src.MakeSS import MakeSS
from src.HLAStudy import Study, HLA_Study

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

# Patterns to use
p_HLA_marker = re.compile(r'^(AA_|HLA_|SNPS_|INS_)')



class SummaryStats(object):

    """

    """

    def __init__(self, _out, _label="", _study=None, _ss=None, _pval=(1, 1e-1, 1e-3, 1e-5, 1e-7, 5e-8)):

        ### Main variables (Initialization)
        self.study = None
        self.ss = None # a list of logistic regression results.
        self.f_haveHLAbmarkers = False


        if bool(_ss) and exists(_ss) and (not bool(_study)):

            self.ss = [_ss]

            self.f_haveHLAbmarkers = pd.read_csv(self.ss[0], header=0, dtype=str).iloc[:, 1] \
                                        .str.match(p_HLA_marker).any()

        else:

            self.study = self.checkStudy(_study)

            self.ss = [(MakeSS(
                self.PLINK_LogisticRegression(self.study.out2+'.p{:.0E}'.format(pval), _a1_allele=self.study.a1_allele,
                                              _allow_no_sex=(not self.study.haveSexInfo), _pfilter=pval),
                self.study.bim,
                self.study.n_case,
                self.study.n_control,
                self.study.out2+'.p{:.0E}.ss'.format(pval)
            )) for pval in _pval]

            self.f_haveHLAbmarkers = isinstance(self.study, HLA_Study)



    def checkStudy(self, _study):

        if not (isinstance(_study, Study) or isinstance(_study, HLA_Study)):
            raise MHCpredError.MHCpredInputPreparationError(
                std_ERROR_MAIN_PROCESS_NAME +
                "Wrong object for TRAIN study."
            )
        else:
            return _study



    def PLINK_LogisticRegression(self, _out, _a1_allele=None, _extract=None, _exclude=None, _keep=None, _remove=None,
                                 _allow_no_sex=True, _ci=0.95, _pfilter=1):

        assoc_result = \
            bashPLINK.LogisticRegression(
                _out, _bed=self.study.bed, _bim=self.study.bim, _fam=self.study.fam,
                _a1_allele=_a1_allele, _extract=_extract, _exclude=_exclude, _keep=_keep, _remove=_remove,
                _allow_no_sex=_allow_no_sex, _ci=_ci, _pfilter=_pfilter)

        print("logsitic regression result : {}".format(assoc_result))

        return assoc_result



    def haveHLAbmarkers(self):
        return self.f_haveHLAbmarkers