#-*- coding: utf-8 -*-

import os, sys, re
import pandas as pd

def MakeSS(_assoc_logistic_beta, _bim_ManualCoord, _n_case, _n_control, _out):


    # (1) Loading '*.assoc.logistic' file (with BETA not OR).
    df_assoc_logistic = pd.read_csv(_assoc_logistic_beta, sep='\s+', header=0, dtype=str)
    # print("df_assoc_logistic:\n{}\n".format(df_assoc_logistic))

    # (2) Loading '*.bim' file (preprocessed by ManualCoord.)
    df_bim = pd.read_csv(_bim_ManualCoord, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'al1', 'al2'])
    # print("df_bim:\n{}\n".format(df_bim))

    df_merge0 = df_assoc_logistic.merge(df_bim[['Label', 'al1', 'al2']], left_on='SNP', right_on='Label')
    # print("df_merge0:\n{}\n".format(df_merge0))


    """
    Needed column
    
    1. Chr
    2. Label
    3. BP
    4. al1
    5. al2 (bim_ManualCoord)
    6. beta
    7. beta_se
    8. n_case
    9. n_control
    10. P-value
    
    => chromosome,marker.ID,physical.pos,allele1,allele2,beta,beta_se,n_case,n_control,p
    """

    df_ss = pd.concat([
        df_merge0[['CHR', 'SNP', 'BP', 'al1', 'al2', 'BETA', 'SE', 'P']],
        pd.Series([str(_n_case)]*df_merge0.shape[0], name='n_case'),
        pd.Series([str(_n_control)]*df_merge0.shape[0], name='n_control'),
    ], axis=1) \
        .reindex(['CHR', 'SNP', 'BP', 'al1', 'al2', 'BETA', 'SE', 'n_case', 'n_control', 'P'], axis=1)

    df_ss.columns = 'chromosome,marker.ID,physical.pos,allele1,allele2,beta,beta_se,n_case,n_control,p'.split(',')
    # print("df_ss:\n{}\n".format(df_ss))

    df_ss.to_csv(_out, header=True, index=False)
    return _out



if __name__ == '__main__':

    """
    MakeSS.py
    
    """

    # [_assoc_logistic_beta, _bim_ManualCoord, _n_case, _n_control, _out] = [
    #     'RA/TRAIN.NBS+RA.06.hg18.28mb-35mb.assoc.logistic',
    #     'RA/TRAIN.NBS+RA.06.hg18.28mb-35mb.bim',
    #     1360,
    #     1458,
    #     'RA/TRAIN.NBS+RA.06.hg18.28mb-35mb.ss'
    # ]

    [_assoc_logistic_beta, _bim_ManualCoord, _n_case, _n_control, _out] = sys.argv[1:]

    MakeSS(_assoc_logistic_beta, _bim_ManualCoord, _n_case, _n_control, _out)