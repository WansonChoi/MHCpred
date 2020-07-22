#-*- coding: utf-8 -*-

import os, sys, re
from os.path import exists, dirname, basename, join
import argparse, textwrap
import pandas as pd

import src.Study as s
import src.MHCpredError as MHCpredError
from src.ManualCoord import ManualCoord

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))



def MHCPred(_validation, _out, _hg, _train=None, _train_ss=None):

    """

    """

    ##### < [] Preprocessing > #####

    ## Loading Genotype data

    # Validation set
    __VAL__ = s.Study(_validation, _out, _label="raw validation data set").WhetherToHLAStudy(_validation)


    # Training set
    if bool(_train_ss) and exists(_train_ss):

        __TRAIN__ = s.TRAIN_Study(_out, _ss=_train_ss, _label="TRAIN given as Summary statistics at first by argument.")

    else:
        __TRAIN__ = s.Study(_train, _out, _label="raw train data set").WhetherToHLAStudy(_train)

        ## ManualCoord
        MC_toExclude, MC_toExtract = ManualCoord( __VAL__.bim, __TRAIN__.bim, __VAL__.bim, join(dirname(_out), 'ManualCoord'))
        print(MC_toExclude)
        print(MC_toExtract)

        __TRAIN__ = __TRAIN__.PLINK_make_bed(join(dirname(_out), basename(_train)+'.MC'), _extract=MC_toExtract, _exclude=MC_toExclude)

        __VAL__ = __VAL__.PLINK_make_bed(join(dirname(_out), basename(_validation)+'.MC'), _a1_allele=__TRAIN__.a1_allele, # set a1 allele to that of TRAIN set.
                                         _extract=MC_toExtract, _exclude=MC_toExclude)


        __TRAIN__ = s.TRAIN_Study(_out, _study=__TRAIN__, _label="raw train data set")


    ##### < [] LDpred2 > #####




    ##### < [] Scoring > #####


    return 0



if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        MHCpred

        (Created by Wanson Choi, 2020/04/03)

        


    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    ### Common arguments to share over the modules.

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')
    parser.add_argument("--out", help="\nOutput file name prefix\n\n", required=True)
    parser.add_argument("--hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"], required=True)

    # parser.add_argument("--use-bmarker", help="\nWhen bmarker is included.\n\n", action="store_true")

    TRAIN = parser.add_mutually_exclusive_group(required=True)
    TRAIN.add_argument("--train", help="\nTraining set file prefix.\n\n")
    TRAIN.add_argument("--train-ss", help="\nSummary statistics file of Training set.\n\n")
    parser.add_argument("--val", help="\nValidation set file prefix.\n\n", required=True)





    ##### < for Testing > #####

    args = parser.parse_args('--out tests/test3/test3 '
                             '--hg 19 '
                             '--train tests/example/RA/TRAIN.RA+58C+NBS.06.hg19.27892021-34892022 '
                             '--val tests/example/RA/VALIDATION.RA+58C+NBS.06.hg19.27892021-34892022'.split(' '))

    # args = parser.parse_args('--out tests/test2/test2 '
    #                          '--hg 19 '
    #                          '--ref-bfile tests/example/RA2/wtccc_filtered_58C+NBS_06.hg19.27892021-34892022.CONTROL '
    #                          '--train-ss tests/example/RA2/wtccc_filtered_NBS+RA_06.hg19.27892021-34892022.TRAIN.assoc.logistic.NoNA.ss '
    #                          '--train-N 2818 '
    #                          '--test-bfile tests/example/RA2/wtccc_filtered_58C+RA_06.hg19.27892021-34892022.TEST '
    #                          '--test-pheno tests/example/RA2/wtccc_filtered_58C+RA_06.hg19.27892021-34892022.TEST.phe'.split(' '))


    ##### < for Publish > #####
    # args = parser.parse_args()
    print(args)

    MHCPred(args.val, args.out, args.hg, _train=args.train, _train_ss=args.train_ss)