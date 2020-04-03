#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


def MHCPred(_REF, _TRAIN, _TEST, _out, _hg):

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



    ### Arguments for each data set.

    # Reference data set
    parser.add_argument("--ref-bfile", help="\nPrefix of reference genotype data set. (*.{bed,bim,fam})\n\n")

    # Train data set
    """
    *(1) Summary statistic file(*.ss; Primarily Logistic regression result) => for now, 
    (2) Genotype data (*.{bed,bim,fam}) with other information (ex. *.phe, *.a1_allele, Phenotype name)
    (3) All required files are given separately.
    
    (2), (3) => Introduce this later when 'ManualCoord' is integrated.
    """
    # parser.add_argument("--input-train", help="\nPrefix of reference genotype data set. (*.{bed,bim,fam,a1_allele})\n\n")
    # parser.add_argument("--train-bfile", help="\nPrefix of train genotype data set. (*.{bed,bim,fam})\n\n")
    parser.add_argument("--train-ss", help="\nPrefix of train summary statistic file.\n\n")

    # Test data set
    parser.add_argument("--test-bfile", help="\nPrefix of test genotype data set. (*.{bed,bim,fam})\n\n")
    parser.add_argument("--test-pheno", help="\nPrefix of test phenotype file. (*.phe)\n\n")
    # parser.add_argument("--pheno-name-test", help="\nPrefix of train genotype data set. (*.{bed,bim,fam})\n\n")





    ##### < for Testing > #####

    # args = parser.parse_args('--omnibus --out /Users/wansun/Git_Projects/HATK/tests/_4_HLA_Analysis/Omnibus/merged/20190618_merged_covar --input /Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/data/HLA_Analysis/sample/Merged/merged --covar-name GWAS --pheno-name All_UC'.split(' '))


    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)