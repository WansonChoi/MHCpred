#-*- coding: utf-8 -*-

import os, sys, re
from os.path import exists
import subprocess
from shutil import which

import src.MHCpredError as MHCpredError

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))



p_PLINK = None

if which('plink'):
    p_PLINK = which('plink')
elif exists('dependency/plink'):
    p_PLINK = 'dependency/plink'
else:
    raise MHCpredError.MHCpredInputPreparationError(
        std_ERROR_MAIN_PROCESS_NAME +
        'No plink is found. Please install PLINK(1.9b) in your system or \'dependency/\' folder.')



def MakeBed(_out, _bfile=None, _bed=None, _bim=None, _fam=None, _a1_allele=None,
            _extract=None, _exclude=None, _keep=None,
            _remove=None, _allow_no_sex=True):


    command = [p_PLINK, '--make-bed']

    if _bfile:
        command.append('--bfile {}'.format(_bfile))
    elif _bed and _bim and _fam:
        command.append('--bed {} --bim {} --fam {}'.format(_bed, _bim, _fam))
    else:
        raise MHCpredError.MHCpredBashExecutionError(
            std_ERROR_MAIN_PROCESS_NAME +
            "Inappropriate input for PLINK '--make-bed' command.\n"
            "('--bifle' : {}, '--bed' : {}, '--bim' : {}, '--fam' : {})".format(
                _bfile, _bed, _bim, _fam
            )
        )

    if _out:
        command.append('--out {}'.format(_out))
    else:
        raise MHCpredError.MHCpredBashExecutionError(
            std_ERROR_MAIN_PROCESS_NAME +
            "Inappropriate output prefix for PLINK '--make-bed' command.\n"
            "('--out' : {})".format(_out)
        )

    if _a1_allele:
        command.append('--a1-allele {}'.format(_a1_allele))

    if _extract:
        command.append('--extract {}'.format(_extract))

    if _exclude:
        command.append('--exclude {}'.format(_exclude))

    if _keep:
        command.append('--keep {}'.format(_keep))

    if _remove:
        command.append('--remove {}'.format(_remove))

    if _allow_no_sex:
        command.append('--allow-no-sex')


    # print(command)
    command = ' '.join(command)

    try:
        subprocess.run(re.split('\s+', command), check=True, stdout=subprocess.DEVNULL)

    except subprocess.CalledProcessError:
        raise MHCpredError.MHCpredBashExecutionError(
            std_ERROR_MAIN_PROCESS_NAME +
            "PLINK '--make-bed' failed.\n"
            "(command : '{}')".format(command)
        )
    else:
        return _out



def LogisticRegression(_out, _bfile=None, _bed=None, _bim=None, _fam=None, _a1_allele=None,
                       _extract=None, _exclude=None, _keep=None, _remove=None, _allow_no_sex=True,
                       _ci=0.95, _pfilter=1):


    command = [p_PLINK, '--logistic hide-covar beta', '--pfilter {}'.format(_pfilter), '--ci {}'.format(_ci)]

    if _bfile:
        command.append('--bfile {}'.format(_bfile))
    elif _bed and _bim and _fam:
        command.append('--bed {} --bim {} --fam {}'.format(_bed, _bim, _fam))
    else:
        raise MHCpredError.MHCpredBashExecutionError(
            std_ERROR_MAIN_PROCESS_NAME +
            "Inappropriate input for PLINK '--make-bed' command.\n"
            "('--bifle' : {}, '--bed' : {}, '--bim' : {}, '--fam' : {})".format(
                _bfile, _bed, _bim, _fam
            )
        )

    if _out:
        command.append('--out {}'.format(_out))
    else:
        raise MHCpredError.MHCpredBashExecutionError(
            std_ERROR_MAIN_PROCESS_NAME +
            "Inappropriate output prefix for PLINK '--make-bed' command.\n"
            "('--out' : {})".format(_out)
        )

    if _a1_allele:
        command.append('--a1-allele {}'.format(_a1_allele))

    if _extract:
        command.append('--extract {}'.format(_extract))

    if _exclude:
        command.append('--exclude {}'.format(_exclude))

    if _keep:
        command.append('--keep {}'.format(_keep))

    if _remove:
        command.append('--remove {}'.format(_remove))

    if _allow_no_sex:
        command.append('--allow-no-sex')



    # print(command)
    command = ' '.join(command)

    try:
        subprocess.run(re.split('\s+', command), check=True, stdout=subprocess.DEVNULL)

    except subprocess.CalledProcessError:
        raise MHCpredError.MHCpredBashExecutionError(
            std_ERROR_MAIN_PROCESS_NAME +
            "PLINK '--make-bed' failed.\n"
            "(command : '{}')".format(command)
        )
    else:
        return _out+'.assoc.logistic'