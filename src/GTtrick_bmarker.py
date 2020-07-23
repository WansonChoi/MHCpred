import os, sys, re
import pandas as pd

p_bmarker = re.compile(r'^(AA_|HLA_|SNPS_|INS_)')


def GTtrick_bmarker(_bim, _out=None):
    
    __bim__ = pd.read_csv(_bim, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'al1', 'al2'])
    # print("__bim__:\n{}\n".format(__bim__))
    
    
    f_bmarker = __bim__['Label'].str.match(p_bmarker)
    # print("Before:\n{}\n".format(__bim__[f_bmarker]))
    
    __bim__.loc[f_bmarker, 'al1'] = 'G'
    __bim__.loc[f_bmarker, 'al2'] = 'T'
    # print("AFter:\n{}\n".format(__bim__[f_bmarker]))
    
    
    if bool(_out):
        __bim__.to_csv(_out, sep='\t', header=False, index=False)
        return _out
    else:
        return __bim__
    
    
    
if __name__ == '__main__':
    
    
    [_bim, _out] = sys.argv[1:]
    GTtrick_bmarker(_bim, _out)