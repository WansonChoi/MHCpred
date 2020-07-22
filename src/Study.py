#-*- coding: utf-8 -*-

import os, sys, re
from os.path import exists, dirname, basename, join
import pandas as pd

import src.MHCpredError as MHCpredError

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))



class Study(object):

    def __init__(self, _input, _out, _a1_allele=None, _pheno=None, _pheno_name=None, _covar=None, _covar_name=None,
                 _pcs=None):

        """

        Data

        """

        self.out2 = join(dirname(_out), basename(_input))

        ### Main variable (Initialization)
        self.bed, self.bim, self.fam = self.checkInput(_input)

        self.df_bim = pd.read_csv(self.bim, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'al1', 'al2'])
        self.df_fam = pd.read_csv(self.fam, sep='\s+', header=None, dtype=str, names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phe'])

        self.a1_allele = self.checkA1Allele(_input, _a1_allele)
        self.pheno, self.pheno_name, self.havePheInfoinFAM = self.checkPheInfo(_input, _pheno, _pheno_name)
        self.covar = _covar
        self.covar_name = _covar_name
        self.pcs = _pcs

        self.haveSexInfo = self.checkSexInfo()


        # Printing summary
        self.str_summary = self.getSummary()
        print('\n'.join(self.str_summary))



    def checkInput(self, _input):

        if bool(_input):

            if not exists(_input+'.bed'):
                raise MHCpredError.MHCpredInputPreparationError(std_ERROR_MAIN_PROCESS_NAME + "Input PLINK bed file('{}') can't be found.".format(self.bed))

            if not exists(_input+'.bim'):
                raise MHCpredError.MHCpredInputPreparationError(std_ERROR_MAIN_PROCESS_NAME + "Input PLINK bim file('{}') can't be found.".format(self.bim))

            if not exists(_input+'.fam'):
                raise MHCpredError.MHCpredInputPreparationError(std_ERROR_MAIN_PROCESS_NAME + "Input PLINK fam file('{}') can't be found.".format(self.fam))

            return _input+'.bed', _input+'.bim', _input+'.fam'

        else:

            raise MHCpredError.MHCpredInputPreparationError(std_ERROR_MAIN_PROCESS_NAME + "No value is given for Input PLINK file prefix('{}').".format(_input))



    def checkA1Allele(self, _input, _a1_allele):

        if bool(_a1_allele):

            if exists(_a1_allele):
                return _a1_allele
            else:
                raise MHCpredError.MHCpredInputPreparationError(std_ERROR_MAIN_PROCESS_NAME + "Given a1-allele file ('{}') can't be found.".format(_a1_allele))

        else:

            if exists(_input + '.a1_allele'):
                # Given with `_input` as prefix.
                return _input + '.a1_allele'

            else:
                # Generate a new one.
                self.df_bim[['Label', 'al1']].to_csv(self.out2+'.a1_allele', sep='\t', header=False, index=False)

                return self.out2+'.a1_allele'



    def checkPheInfo(self, _input, _pheno, _pheno_name):

        if bool(_pheno):

            # if exists(_pheno):
            #
            #     if bool(_pheno_name):
            #         # pheno_name given
            #
            #     else:
            #         # pheno_name NOT given
            # else:
            #     raise MHCpredError.MHCpredInputPreparationError(std_ERROR_MAIN_PROCESS_NAME + "Given phenotype file ('{}') can't be found.".format(_pheno))

            pass

        else:

            # Not given (ex. all '-9' or '0')
            f_NotGiven = (self.df_fam['Phe'] == '-9').all() or (self.df_fam['Phe'] == '0').all()

            if f_NotGiven:
                raise MHCpredError.MHCpredInputPreparationError(std_ERROR_MAIN_PROCESS_NAME + "No Phenotype information is in the fam file('{}').".format(self.fam))


            # All cases or controls
            f_homogenous = (self.df_fam['Phe'] == '1').all() or (self.df_fam['Phe'] == '2').all()

            if  f_homogenous:
                raise MHCpredError.MHCpredInputPreparationError(std_ERROR_MAIN_PROCESS_NAME + "All samples are either Cases or Controls.".format(_pheno))


            return None, None, True



    def checkCovarInfo(self):
        return 0



    def checkSexInfo(self):
        return ((self.df_fam['Sex'] == '1') | (self.df_fam['Sex'] == '2')).all()



    def PLINK_LogisticRegression(self):
        return 0



    def PLINK_make_bed(self, _out, _extract=None, _exclude=None, _keep=None, _remove=None, _allow_no_sex=True):

        return 0



    def getSummary(self):

        str_summary = ["bed : {}\nbim : {}\nfam : {}".format(self.bed, self.bim, self.fam)]

        if bool(self.a1_allele):
            str_summary.append("a1_allele : {}".format(self.a1_allele))

        if not self.havePheInfoinFAM:
            str_summary.append("phenotype : {}".format(self.pheno))
            str_summary.append("phenotype-name: {}".format(self.pheno_name))
        else:
            str_summary.append("Phenotype info in fam file.".format(self.pheno_name))

        if bool(self.covar):
            str_summary.append("covariate : {}".format(self.covar))
        if bool(self.covar_name):
            str_summary.append("covariate-name : {}".format(self.covar_name))

        if bool(self.haveSexInfo):
            str_summary.append("Sex info in fam : {}".format(self.haveSexInfo))
        else:
            str_summary.append("No Sex info.")


        return str_summary





class HLA_Study(Study):

    """

    - refineBP.py
    - GTtrick.py
    - AA | HLA | SNPS | ...

    """

    def __init__(self, _input, _out, _a1_allele=None, _pheno=None, _pheno_name=None, _covar=None, _covar_name=None,
                 _pcs=None):

        super().__init__(_input, _out, _a1_allele, _pheno, _pheno_name, _covar, _covar_name, _pcs)






# class TRAIN_study(Study):
#
#     """
#
#     """
#
#     def __init__(self, _input):
#         super().__init__(_input)
#
#         self.checkA1Allele()
#
#
