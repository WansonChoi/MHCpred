import os, sys, re
import pandas as pd

# p_bmarker = re.compile(r'^(AA_|HLA_|SNPS_|INS_)')
p_AA = re.compile(r'^AA_')
p_HLA = re.compile(r'^HLA_')
p_SNPS = re.compile(r'^SNPS_')
p_INS = re.compile(r'^INS_')

VALID_Nucleotides = ['A', 'C', 'G', 'T']
AMBIGUOUS_ALLELE_pair = [{'A', 'T'}, {'G', 'C'}]
reverse_nt = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}


def ManualCoord(_REFERENCE, _TRAIN, _TEST, _out):
    
    REFERENCE_bim = pd.read_csv(_REFERENCE, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'al1', 'al2'], usecols=[1, 4, 5])
#     print("REFERENCE_bim:\n{}\n".format(REFERENCE_bim))
    
    TRAIN_bim = pd.read_csv(_TRAIN, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'al1', 'al2'], usecols=[1, 4, 5])
#     print("TRAIN_bim:\n{}\n".format(TRAIN_bim))
    
    TEST_bim = pd.read_csv(_TEST, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'al1', 'al2'], usecols=[1, 4, 5])
#     print("TEST_bim:\n{}\n".format(TEST_bim))
    
    # Joining 3 tables.
    df_nt_table = REFERENCE_bim.merge(TRAIN_bim, on='Label', suffixes=['_REF', '_TRAIN']) \
                                .merge(TEST_bim, on='Label')
    print("df_nt_table:\n{}\n".format(df_nt_table))
    
    
    
    
    
    
    ### Main iteration ###
    
    l_normal = []
    
    l_a1_allele_TRAIN = []
    l_a1_allele_TEST = []
    l_toFlip_TRAIN = []
    l_toFlip_TEST = []
    
    l_OtherDiscrep = [] # Trash bin.
    l_ambiguous = []
    
    count = 0
    
    
    
    for line in df_nt_table.itertuples():
        
#         print(line)
        index = line[0]
        rs_id = line[1]
        REF_al1, REF_al2 = line[2], line[3]
        TRAIN_al1, TRAIN_al2 = line[4], line[5]
        TEST_al1, TEST_al2 = line[6], line[7]
        
        
        
        # Checking bmarker
        f_AA = bool(p_AA.match(rs_id))
        f_HLA = bool(p_HLA.match(rs_id))
        f_SNPS = bool(p_SNPS.match(rs_id)) and (not (REF_al1 in VALID_Nucleotides and REF_al2 in VALID_Nucleotides))
        f_INS = bool(p_INS.match(rs_id))
        f_bmarker = f_AA or f_HLA or f_SNPS or f_INS
        
        
        
        if not f_bmarker:
            
            ### [1] Normal SNP(ex. rs1234) or [2] HLA SNP(ex. SNP_A_12345678)

            if (not (((REF_al1 in VALID_Nucleotides) and (REF_al2 in VALID_Nucleotides)) and 
                    ((TRAIN_al1 in VALID_Nucleotides) and (TRAIN_al2 in VALID_Nucleotides)) and 
                    ((TEST_al1 in VALID_Nucleotides) and (TEST_al2 in VALID_Nucleotides)))):
                
                # (1) Valid allele?

                l_OtherDiscrep.append(index)
                
                

            elif ((set([REF_al1, REF_al2]) in AMBIGUOUS_ALLELE_pair) or (set([TRAIN_al1, TRAIN_al2]) in AMBIGUOUS_ALLELE_pair) or (set([TEST_al1, TEST_al2]) in AMBIGUOUS_ALLELE_pair)):
                
                # (2) Ambiguous allele? 

                l_ambiguous.append(index)

                

            else:
                # (3-1) Ref / Alt relationship
                # (3-2) Flip relationship


                g_nt = [REF_al1, REF_al2] # REF allele pair
                os_g_nt = [reverse_nt[REF_al1], reverse_nt[REF_al2]]

                ss_nt = [TRAIN_al1, TRAIN_al2]
                t_nt = [TEST_al1, TEST_al2]


                # Naive version of LDpred conditions.

    #             if ((g_nt == ss_nt) or (os_g_nt == ss_nt)) and ((g_nt == t_nt) or (os_g_nt == t_nt)):
    #                 l_normal.append(index)
    #             else:
    #                 # FLIP ?

    #                 # (1) Flip between reference and Train
    #                 if (((g_nt[1] == ss_nt[0]) and (g_nt[0] == ss_nt[1])) or ((os_g_nt[1] == ss_nt[0]) and (os_g_nt[0] == ss_nt[1]))):
    #                     l_normal.append(index)
    #                     l_toFlip_TRAIN.append(index)
    #                 else:
    #                     l_OtherDiscrep.append(index)

    #                 # (2) Flip between reference and Test
    #                 if (((g_nt[1] == t_nt[0]) and (g_nt[0] == t_nt[1])) or ((os_g_nt[1] == t_nt[0]) and (os_g_nt[0] == t_nt[1]))):
    #                     l_normal.append(index)
    #                     l_toFlip_TEST.append(index)
    #                 else:
    #                     l_OtherDiscrep.append(index)



                if ((g_nt == ss_nt) and (g_nt == t_nt)):
                    # No problem
                    l_normal.append(index)

                elif ((not (g_nt == ss_nt)) and (g_nt == t_nt)):

                    # Check TRAIN data

                    if (g_nt[1] == ss_nt[0]) and (g_nt[0] == ss_nt[1]):
                        # Change between REF and ALT relationship
                        l_a1_allele_TRAIN.append(index)
                        l_normal.append(index)
                    elif (os_g_nt[1] == ss_nt[0]) and (os_g_nt[0] == ss_nt[1]):
                        # Reversed Flip relationship
                        l_toFlip_TRAIN.append(index)
                        l_normal.append(index)
    #                     l_a1_allele_TRAIN.append(index)
                    else:
                        l_OtherDiscrep.append(index) # Trash can

                elif ((g_nt == ss_nt) and (not (g_nt == t_nt))):

                    # Check Test data

                    if (g_nt[1] == t_nt[0]) and (g_nt[0] == t_nt[1]):
                        # Change between REF and ALT relationship
                        l_a1_allele_TEST.append(index)
                        l_normal.append(index)
                    elif (os_g_nt[1] == t_nt[0]) and (os_g_nt[0] == t_nt[1]):
                        # Reversed Flip relationship
                        l_toFlip_TEST.append(index)
                        l_normal.append(index)
                    else:
                        l_OtherDiscrep.append(index) # Trash can


                else:

                    if ((g_nt[1] == ss_nt[0]) and (g_nt[0] == ss_nt[1])) and ((g_nt[1] == t_nt[0]) and (g_nt[0] == t_nt[1])):
                        # Change between REF and ALT relationship
                        l_a1_allele_TRAIN.append(index)
                        l_a1_allele_TEST.append(index)
                        l_normal.append(index)

                    elif ((os_g_nt[1] == ss_nt[0]) and (os_g_nt[0] == ss_nt[1])) and ((os_g_nt[1] == t_nt[0]) and (os_g_nt[0] == t_nt[1])):
                        l_toFlip_TRAIN.append(index)
                        l_toFlip_TEST.append(index)
                        l_normal.append(index)

                    else:
                        l_OtherDiscrep.append(index) # Trash can

                        
                        
        else:
            
            ### bmarker (AA, HLA, INS)

            """
            - Ref / Alt relationship
            """

            g_nt = [REF_al1, REF_al2] # REF allele pair
            ss_nt = [TRAIN_al1, TRAIN_al2]
            t_nt = [TEST_al1, TEST_al2]



            if ((g_nt == ss_nt) and (g_nt == t_nt)):
                # No problem
                l_normal.append(index)

            elif ((not (g_nt == ss_nt)) and (g_nt == t_nt)):

                # Check TRAIN data

                if (g_nt[1] == ss_nt[0]) and (g_nt[0] == ss_nt[1]):
                    # Change between REF and ALT relationship
                    l_a1_allele_TRAIN.append(index)
                    l_normal.append(index)
                else:
                    l_OtherDiscrep.append(index) # Trash can

            elif ((g_nt == ss_nt) and (not (g_nt == t_nt))):

                # Check Test data

                if (g_nt[1] == t_nt[0]) and (g_nt[0] == t_nt[1]):
                    # Change between REF and ALT relationship
                    l_a1_allele_TEST.append(index)
                    l_normal.append(index)
                else:
                    l_OtherDiscrep.append(index) # Trash can


            else:

                # Check both TRAIN and TEST
                if ((g_nt[1] == ss_nt[0]) and (g_nt[0] == ss_nt[1])) and ((g_nt[1] == t_nt[0]) and (g_nt[0] == t_nt[1])):
                    # Change between REF and ALT relationship
                    l_a1_allele_TRAIN.append(index)
                    l_a1_allele_TEST.append(index)
                    l_normal.append(index)
                else:
                    l_OtherDiscrep.append(index) # Trash can

            

                
        
        count += 1
#         if count > 30 : break
            
            
#     print("Other Discrepency:\n{}".format(l_OtherDiscrep))
#     print("Ambiguous:\n{}".format(l_ambiguous))
#     print("To Flip (TRAIN):\n{}".format(l_toFlip_TRAIN))
#     print("To Flip (TEST):\n{}".format(l_toFlip_TEST))
#     print("a1-allele (TRAIN):\n{}".format(l_a1_allele_TRAIN))
#     print("a1-allele (TEST):\n{}".format(l_a1_allele_TEST))
#     print("Normal:\n{}".format(l_normal))
    

#     print(":\n{}\n".format())

    if len(l_OtherDiscrep) > 0:
        sr_OtherDiscrepency = df_nt_table.iloc[l_OtherDiscrep, 0]
        print("sr_OtherDiscrepency:\n{}\n".format(sr_OtherDiscrepency.head()))
        sr_OtherDiscrepency.to_csv(_out+'.OtherDiscrepency.ToExclude.txt', header=False, index=False)
    else:
        print("No other discrepency")
        
    if len(l_ambiguous) > 0:
        sr_Ambiguous = df_nt_table.iloc[l_ambiguous, 0]
        print("sr_Ambiguous:\n{}\n".format(sr_Ambiguous.head()))
        sr_Ambiguous.to_csv(_out+'.Ambiguous.ToExclude.txt', header=False, index=False)
    else:
        print("No ambiguous")
        
    if len(l_toFlip_TRAIN) > 0:
        sr_toFlip_TRAIN = df_nt_table.iloc[l_toFlip_TRAIN, 0]
        print("sr_toFlip_TRAIN:\n{}\n".format(sr_toFlip_TRAIN.head()))
        sr_toFlip_TRAIN.to_csv(_out+'.TRAIN.ToFLIP.txt', header=False, index=False)
    else:
        print("No toFlip_TRAIN")

    if len(l_toFlip_TEST) > 0:
        sr_toFlip_TEST = df_nt_table.iloc[l_toFlip_TEST, 0]
        print("sr_toFlip_TEST:\n{}\n".format(sr_toFlip_TEST.head()))
        sr_toFlip_TEST.to_csv(_out+'.TEST.ToFLIP.txt', header=False, index=False)
    else:
        print("No toFlip_TEST")
        
    if len(l_a1_allele_TRAIN) > 0:
        sr_a1_allele_TRAIN = df_nt_table.iloc[l_a1_allele_TRAIN, [0, 1]]
        print("sr_a1_allele_TRAIN:\n{}\n".format(sr_a1_allele_TRAIN.head()))
        sr_a1_allele_TRAIN.to_csv(_out+'.TRAIN.a1_allele.txt', sep='\t', header=False, index=False)
    else:
        print("No a1_allele_TRAIN")
        
    if len(l_a1_allele_TEST) > 0:
        sr_a1_allele_TEST = df_nt_table.iloc[l_a1_allele_TEST, [0, 1]]
        print("sr_a1_allele_TEST:\n{}\n".format(sr_a1_allele_TEST.head()))
        sr_a1_allele_TEST.to_csv(_out+'.TEST.a1_allele.txt', sep='\t', header=False, index=False)
    else:
        print("No a1_allele_TEST")
        
    if len(l_normal) > 0:
        sr_normal = df_nt_table.iloc[l_normal, 0]
        print("sr_normal:\n{}\n".format(sr_normal.head()))
        sr_normal.to_csv(_out+'.Normal.ToExtract.txt', header=False, index=False)
    else:
        print("No sr_normal")
        
    
    
    
    
        
        

    
if __name__ == '__main__':
    
    
    [_REFERENCE_bim, _TRAIN_bim, _TEST_bim, _out] = sys.argv[1:]
    
    ManualCoord(_REFERENCE_bim, _TRAIN_bim, _TEST_bim, _out)