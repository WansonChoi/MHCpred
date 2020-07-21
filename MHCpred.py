#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


def MHCPred(_out, _hg):

    """

    """

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

    parser.add_argument("--use-bmarker", help="\nWhen bmarker is included.\n\n", action="store_true")





    ##### < for Testing > #####

    # args = parser.parse_args('--out tests/test1/test1 '
    #                          '--hg 19 '
    #                          '--ref-bfile tests/example/RA/wtccc_filtered_58C+NBS_06.hg18.28mb-35mb.MC.CONTROL '
    #                          '--train-ss tests/example/RA/wtccc_filtered_NBS+RA_06.hg18.28mb-35mb.MC.TRAIN.NoNA.ss '
    #                          '--train-N 2818 '
    #                          '--test-bfile tests/example/RA/wtccc_filtered_58C+RA_06.hg18.28mb-35mb.MC.TEST '
    #                          '--test-pheno tests/example/RA/wtccc_filtered_58C+RA_06.hg18.28mb-35mb.MC.TEST.phe'.split(' '))

    # args = parser.parse_args('--out tests/test2/test2 '
    #                          '--hg 19 '
    #                          '--ref-bfile tests/example/RA2/wtccc_filtered_58C+NBS_06.hg19.27892021-34892022.CONTROL '
    #                          '--train-ss tests/example/RA2/wtccc_filtered_NBS+RA_06.hg19.27892021-34892022.TRAIN.assoc.logistic.NoNA.ss '
    #                          '--train-N 2818 '
    #                          '--test-bfile tests/example/RA2/wtccc_filtered_58C+RA_06.hg19.27892021-34892022.TEST '
    #                          '--test-pheno tests/example/RA2/wtccc_filtered_58C+RA_06.hg19.27892021-34892022.TEST.phe'.split(' '))


    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)

    MHCPred(args.out, args.hg)