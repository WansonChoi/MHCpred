#-*- coding: utf-8 -*-

import os, sys, re
import subprocess
import argparse, textwrap

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


def MHCPred(_REF, _TRAIN, _TEST, _out, _hg,
            _N, _TEST_phe=None):



    """

    """


    ##### < 1. Coord step > #####

    print(std_MAIN_PROCESS_NAME + "##### < 1. Coord step > #####")

    command = "python -m ldpred --debug coord " \
              "--gf {} " \
              "--ssf {} " \
              "--out {} " \
              "--N {} --ssf-format CUSTOM --eff OR --rs SNP --pos BP --pval P".format(_REF, _TRAIN, _out+'.coord', _N)
    # print(command)

    try:
        subprocess.run(command.split(' '), check=True, stdout=open(_out+'.coord.debug.log', 'w'), stderr=open(_out+'.coord.err.log', 'w'))
        # os.system(command)
    except subprocess.CalledProcessError:
        print(std_ERROR_MAIN_PROCESS_NAME + "Coord step failed. Please check '{}' file.".format(_out+'.coord.err.log'))
        return -1
    else:

        os.system("rm {}".format(_out+'.coord.err.log'))


        # (1) *.coord file.
        __coord__ = _out+'.coord'


        # (2) LD radius value
        ldr = getLDR(_out+'.coord.debug.log')

        if ldr < 0:
            print(std_ERROR_MAIN_PROCESS_NAME + "Failed to get LDR value.")
        else:
            print("LD radius : {}".format(ldr))



    ##### < 2. Gibbs step > #####

    print(std_MAIN_PROCESS_NAME + "##### < 2. Gibbs step > #####")

    command = "python -m ldpred --debug gibbs " \
              "--cf {} " \
              "--ldr {} " \
              "--ldf {} " \
              "--out {} " \
              "--N {}".format(__coord__, ldr, _out+'.LD', _out+'.GIBBS', _N)
    print(command)

    try:
        subprocess.run(command.split(' '), check=True, stdout=open(_out+'.GIBBS.log', 'w'), stderr=open(_out+'.GIBBS.err.log', 'w'))

    except subprocess.CalledProcessError:
        print(std_ERROR_MAIN_PROCESS_NAME + "Gibbs step failed. Please check '{}' file.".format(_out+'.GIBBS.err.log', 'w'))
        return -1
    else:

        os.system("rm {}".format(_out+'.GIBBS.err.log'))


        # (1) *.GIBBS
        __rf__ = _out+'.GIBBS'



    ##### < 3. Score step > #####

    print(std_MAIN_PROCESS_NAME + "##### < 3. Score step > #####")

    command = "python -m ldpred --debug score " \
              "--gf {} " \
              "--rf {} " \
              "--out {} " \
              "--pf {} " \
              "--pf-format LSTANDARD " \
              "--summary-file {}".format(_TEST, __rf__, _out+'.SCORE', _TEST_phe, _out+'.SCORE.summaryfile.txt')
    print(command)

    try:
        subprocess.run(command.split(' '), check=True, stdout=open(_out+'.SCORE.log', 'w'), stderr=open(_out+'.SCORE.err.log', 'w'))

    except subprocess.CalledProcessError:
        print(std_ERROR_MAIN_PROCESS_NAME + "Score step failed. Please check '{}' file.".format(_out+'.SCORE.err.log'))
        return -1
    else:

        os.system("rm {}".format(_out+'.SCORE.err.log'))



    print(std_MAIN_PROCESS_NAME + "MHCpred done.")

    return 0



def getLDR(_coord_log):

    total = -1
    ambiguous = -1
    unknown_unsupported = -1
    other_discrepency = -1
    MAF = -1
    AF_discrepency = -1

    p_total = re.compile(r'^SNPs retained after filtering:\s+(\d+)$')
    p_ambiguous = re.compile(r'^SNPs w ambiguous nucleotides filtered:\s+(\d+)$')
    p_unknown_unsupported = re.compile(r'^SNPs w unknown/unsupported nucleotides filtered:\s+(\d+)$')
    p_other_discrepency = re.compile(r'^SNPs w other nucleotide discrepancies filtered:\s+(\d+)$')
    p_MAF = re.compile(r'^SNPs w MAF<0.010 filtered:\s+(\d+)$')
    p_AF_discrepency = re.compile(r'^SNPs w allele freq discrepancy > 0.100 filtered:\s+(\d+)$')


    with open(_coord_log, 'r') as f_coord_log:

        for line in f_coord_log:

            m_total = p_total.match(line)
            m_ambiguous = p_ambiguous.match(line)
            m_unknown_unsupported = p_unknown_unsupported.match(line)
            m_other_discrepency = p_other_discrepency.match(line)
            m_MAF = p_MAF.match(line)
            m_AF_discrepency = p_AF_discrepency.match(line)



            if bool(m_total):
                total = int(m_total.group(1))

            if bool(m_ambiguous):
                ambiguous = int(m_ambiguous.group(1))

            if bool(m_unknown_unsupported):
                unknown_unsupported = int(m_unknown_unsupported.group(1))

            if bool(m_other_discrepency):
                other_discrepency = int(m_other_discrepency.group(1))

            if bool(m_MAF):
                MAF = int(m_MAF.group(1))

            if bool(m_AF_discrepency):
                AF_discrepency = int(m_AF_discrepency.group(1))


    flag_fine = True

    if total == -1:
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to catch the value of \"SNPs retained after filtering\".")
        flag_fine = False

    if ambiguous == -1:
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to catch the value of \"SNPs w ambiguous nucleotides filtered\".")
        flag_fine = False

    if unknown_unsupported == -1:
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to catch the value of \"SNPs w unknown/unsupported nucleotides filtered\".")
        flag_fine = False

    if other_discrepency == -1:
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to catch the value of \"SNPs w other nucleotide discrepancies filtered\".")
        flag_fine = False

    if MAF == -1:
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to catch the value of \"SNPs w MAF<0.010 filtered\".")
        flag_fine = False

    if AF_discrepency == -1:
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to catch the value of \"SNPs w allele freq discrepancy > 0.100 filtered\".")
        flag_fine = False



    if flag_fine:
        return (total - ambiguous - unknown_unsupported - other_discrepency - MAF - AF_discrepency)
    else:
        return -1





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
    parser.add_argument("--train-N", help="\nThe number of samples(patients) used in train summary statistic file.\n\n")

    # Test data set
    parser.add_argument("--test-input", help="\nPre-consented prefix of test data set. (*.{bed,bim,fam,phe})\n\n")
    parser.add_argument("--test-bfile", help="\nPrefix of test genotype data set. (*.{bed,bim,fam})\n\n")
    parser.add_argument("--test-pheno", help="\nPrefix of test phenotype file. (*.phe)\n\n")
    # parser.add_argument("--pheno-name-test", help="\nPrefix of train genotype data set. (*.{bed,bim,fam})\n\n")





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

    MHCPred(args.ref_bfile, args.train_ss, args.test_bfile, args.out, args.hg,
            args.train_N, _TEST_phe=args.test_pheno)