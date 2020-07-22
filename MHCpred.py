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


# Patterns to use
p_HLA_marker = re.compile(r'^(AA_|HLA_|SNPS_|INS_)')



def MHCPred(_train, _validation, _out, _hg):

    """

    """

    ##### < [] Preprocessing > #####

    ## Loading Genotype data

    # contains HLA bmarker?
    f_bmarker_TRAIN = useBmarker(_train+'.bim')
    f_bmarker_VAL = useBmarker(_validation+'.bim')

    if (not f_bmarker_TRAIN) and (not f_bmarker_VAL):

        # Plain genotype data set without HLA bmarkers.

        __TRAIN__ = s.Study(_train, _out)
        __VAL__ = s.Study(_validation, _out)

    elif f_bmarker_TRAIN and f_bmarker_VAL:

        # Genotype data set with HLA bmarkers.

        __TRAIN__ = s.HLA_Study(_train, _out)
        __VAL__ = s.HLA_Study(_validation, _out)

    else:
        raise MHCpredError.MHCpredInputPreparationError(
            std_ERROR_MAIN_PROCESS_NAME +
            "For TRAIN and VALIDATION data sets, Both of them or Neither of them have HLA binary markers "
        )


    ## ManualCoord




    ##### < [] LDpred2 > #####




    ##### < [] Scoring > #####


    return 0



def useBmarker(_bim):

    return pd.read_csv(_bim, sep='\s+', header=None, dtype=str).iloc[:, 1] \
                .str.match(p_HLA_marker) \
                .any()



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

    parser.add_argument("--train", help="\nTraining set file prefix.\n\n", required=True)
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

    MHCPred(args.train, args.val, args.out, args.hg)