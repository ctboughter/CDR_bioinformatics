
# coding: utf-8

# # A Nice Little Script to Assign CDR Loops
# Hopefully this scipt will take a heavy or light antibody chain as input, and output assigend CDR loops
# A lot of other groups seem to do something similar, but this might be a more simplistic way of doing it.

# In[133]:


import numpy as np
import pandas
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import regex


# In[264]:


# Load in all of our files
IGHV=pandas.read_csv('IMGT_database/human_IGHV.txt',sep='\s+',header=[1])
IGLV=pandas.read_csv('IMGT_database/human_IGLV.txt',sep='\s+',header=[1])
IGKV=pandas.read_csv('IMGT_database/human_IGKV.txt',sep='\s+',header=[1])
IGHJ=pandas.read_csv('IMGT_database/human_IGHJ.txt',sep='\s+',header=[1])
IGLJ=pandas.read_csv('IMGT_database/human_IGLJ.txt',sep='\s+',header=[1])
IGKJ=pandas.read_csv('IMGT_database/human_IGKJ.txt',sep='\s+',header=[1])

ori_IGHJ=pandas.read_csv('IMGT_database/human_IGHJ_IMGT.txt',sep='\s+',header=[1])
ori_IGLJ=pandas.read_csv('IMGT_database/human_IGLJ_IMGT.txt',sep='\s+',header=[1])
ori_IGKJ=pandas.read_csv('IMGT_database/human_IGKJ_IMGT.txt',sep='\s+',header=[1])


# In[165]:


# Separate out Heavy Chain Gene V Segements
HV_A=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2, 0]]
HV_B=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2, 1]]
HV_BC=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2, 2]]
HV_C=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2, 3]]
HV_Cp=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2, 4]] # Cp = C'
HV_CpCpp=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2,5]] # Cpp = C''
HV_Cpp=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2,6]]
HV_D=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2,7]]
HV_E=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2,8]]
HV_F=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2,9]]
HV_FG=[Seq(a, IUPAC.protein) for a in IGHV.iloc[2::2,10]]


# In[168]:


# Separate out Lambda Light Chain V Gene Segements
LV_A=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2, 0]]
LV_B=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2, 1]]
LV_BC=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2, 2]]
LV_C=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2, 3]]
LV_Cp=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2, 4]] # Cp = C'
LV_CpCpp=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2,5]] # Cpp = C''
LV_Cpp=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2,6]]
LV_D=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2,7]]
LV_E=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2,8]]
LV_F=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2,9]]
LV_FG=[Seq(a, IUPAC.protein) for a in IGLV.iloc[2::2,10]]


# In[169]:


# Separate out Kappa Light Chain V Gene Segements
KV_A=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2, 0]]
KV_B=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2, 1]]
KV_BC=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2, 2]]
KV_C=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2, 3]]
KV_Cp=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2, 4]] # Cp = C'
KV_CpCpp=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2,5]] # Cpp = C''
KV_Cpp=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2,6]]
KV_D=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2,7]]
KV_E=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2,8]]
KV_F=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2,9]]
KV_FG=[Seq(a, IUPAC.protein) for a in IGKV.iloc[2::2,10]]


# In[265]:


# Separate out Heavy Chain J Gene Segements
HJ=[Seq(a, IUPAC.protein) for a in IGHJ.iloc[0::2, 0]]
ori_HJ=[Seq(a, IUPAC.protein) for a in ori_IGHJ.iloc[0::2, 0]]


# In[266]:


# Separate out Lambda Light Chain J Gene Segements
LJ=[Seq(a, IUPAC.protein) for a in IGLJ.iloc[0::2, 0]]
ori_LJ=[Seq(a, IUPAC.protein) for a in ori_IGLJ.iloc[0::2, 0]]


# In[268]:


# Separate out Kappa Light Chain J Gene Segements
KJ=[Seq(a, IUPAC.protein) for a in IGKJ.iloc[0::2, 0]]
ori_KJ=[Seq(a, IUPAC.protein) for a in ori_IGKJ.iloc[0::2, 0]]


# In[141]:


pre_heavy=pandas.read_csv('heavy_sequences.txt',sep='\s+',header=None)
pre_light=pandas.read_csv('light_sequences.txt',sep='\s+',header=None)

input_heavy=[Seq(a, IUPAC.protein) for a in pre_heavy.iloc[:, 0]]
input_light=[Seq(a, IUPAC.protein) for a in pre_light.iloc[:, 0]]

# # Work in Progress Note:
# So for some reason, including the truncated version of the genes results in a bad fit... Specifically, the J segment was essentially drowned out to noise so the aligner couldn't properly fit it. This makes sense in retrospect... More residues means we are less susceptible to noise. Does this mean we should also include the IMGT CDR loop sections in these alignments? Probably... but let's leave them out for now
# 
# # Note as I'm thinking about... uh... the above note:
# Might want to do a two-step alignment. Align to get just the CDR1 and CDR2, THEN pick out the region that should ~roughly~ have the CDR3 region, THEN align that to the J segment... That's for a future implementation

# In[363]:


# Here we look at heavy chain comparisons, try to find the BEST fit of some given gene usage.

for k in range(len(input_heavy)):
    heavy_t=input_heavy[k]
    score=0
    best=0
    for i in range(len(HV_A)):
        for j in range(len(HJ)):
            aa=pairwise2.align.globalxx(HV_A[i]+HV_B[i]+HV_C[i]+HV_Cp[i]+HV_Cpp[i]+HV_D[i]+HV_E[i]+HV_F[i]+ori_HJ[j],heavy_t)
            if aa[0][2] > score:
                # Bring back these options to look at the second best
                #best_old=best
                #score_old=score
                score=aa[0][2]
                best=[i,j]
    #print(score)
    #print('V gene # and J gene # (IMGT List Number)')
    #print(best)
    final=pairwise2.align.globalxx(HV_A[best[0]]+HV_B[best[0]]+HV_C[best[0]]+HV_Cp[best[0]]+HV_Cpp[best[0]]+HV_D[best[0]]+HV_E[best[0]]+HV_F[best[0]]+HJ[best[1]],heavy_t)
    #print(format_alignment(*final[0]))
# So here we leverage the fact that the alignments are the exact same length.
# Print the in betweens... CDR1 = Between B and C regions... CDR2= Between C' and C''... 
# CDR3= Between F domain and entire J gene usage.

# Need to have this regex ditty in there because of some slight disagreements for
# 
    print('Heavy Sequence '+str(k))
    find_B_tol=regex.findall(str('('+HV_B[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
    error=0
    while find_B_tol==[]:
        error=error+1
        find_B_tol=regex.findall(str('('+HV_B[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
        if error >= 4:
            print('High Mismatch with Conserved Regions, Check for Problems (H_B)')
    find_C_tol=regex.findall(str('('+HV_C[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
    error=0
    while find_C_tol==[]:
        error=error+1
        find_C_tol=regex.findall(str('('+HV_C[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
        if error >= 4:
            print('High Mismatch with Conserved Regions, Check for Problems (H_C)')
    CDR1_start=final[0][0].find(find_B_tol[0])+len(find_B_tol[0])
    CDR1_end=final[0][0].find(find_C_tol[0])
    CDR1=final[0][1][CDR1_start:CDR1_end].replace('-','')
    print(CDR1)

    find_Cp_tol=regex.findall(str('('+HV_Cp[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
    error=0
    while find_Cp_tol==[]:
        error=error+1
        find_Cp_tol=regex.findall(str('('+HV_Cp[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
        if error >= 4:
            print('High Mismatch with Conserved Regions, Check for Problems (H_Cp)')
    find_Cpp_tol=regex.findall(str('('+HV_Cpp[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
    error=0
    while find_Cpp_tol==[]:
        error=error+1
        find_Cpp_tol=regex.findall(str('('+HV_Cpp[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
        if error >= 4:
            print('High Mismatch with Conserved Regions, Check for Problems (H_Cpp')
    CDR2_start=final[0][0].find(find_Cp_tol[0])+len(find_Cp_tol[0])
    CDR2_end=final[0][0].find(find_Cpp_tol[0])
    CDR2=final[0][1][CDR2_start:CDR2_end].replace('-','')
    print(CDR2)

    find_F_tol=regex.findall(str('('+HV_F[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
    error=0
    while find_F_tol==[]:
        error=error+1
        find_F_tol=regex.findall(str('('+HV_F[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
        if error >= 4:
            print('High Mismatch with Conserved Regions, Check for Problems (H_F)')
    find_J_tol=regex.findall(str('('+HJ[best[1]])+'){i<=0,s<=0,e<=0}',final[0][0])
    error=0
    while find_J_tol==[]:
        error=error+1
        find_J_tol=regex.findall(str('('+HJ[best[1]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
        if error >= 4:
            print('High Mismatch with Conserved Regions, Check for Problems (HJ)')
    CDR3_start=final[0][0].find(find_F_tol[0])+len(find_F_tol[0])
    CDR3_end=final[0][0].find(find_J_tol[0])
    CDR3=final[0][1][CDR3_start:CDR3_end].replace('-','')
    print(CDR3)


# In[346]:


# Alright now do the same stuff for the light chain
# lambda time
for k in range(len(input_light)):
    light_t=input_light[k]
    score=0
    best=0
    for i in range(len(LV_A)):
        for j in range(len(LJ)):
            aa=pairwise2.align.globalxx(LV_A[i]+LV_B[i]+LV_C[i]+LV_Cp[i]+LV_Cpp[i]+LV_D[i]+LV_E[i]+LV_F[i]+ori_LJ[j],light_t)
            if aa[0][2] > score:
                # Bring back these options to look at the second best
                #best_old=best
                #score_old=score
                score=aa[0][2]
                best=[i,j]
                lc='lambda'
    lambdascore=score
    # Light is a little bit trickier, because you have to check the kappa 
    for i in range(len(KV_A)):
        for j in range(len(KJ)):
            aa=pairwise2.align.globalxx(KV_A[i]+KV_B[i]+KV_C[i]+KV_Cp[i]+KV_Cpp[i]+KV_D[i]+KV_E[i]+KV_F[i]+ori_KJ[j],light_t)
            if aa[0][2] >= score:
                # Bring back these options to look at the second best
                #best_old=best
                #score_old=score
                score=aa[0][2]
                best=[i,j]
                lc='kappa'
    if score==lambdascore:
        print('We have a tie between lambda and kappa genes! This should not really happen')

    if lc=='lambda':
        final=pairwise2.align.globalxx(LV_A[best[0]]+LV_B[best[0]]+LV_C[best[0]]+LV_Cp[best[0]]+LV_Cpp[best[0]]+LV_D[best[0]]+LV_E[best[0]]+LV_F[best[0]]+LJ[best[1]],light_t)
    if lc=='kappa':
        final=pairwise2.align.globalxx(KV_A[best[0]]+KV_B[best[0]]+KV_C[best[0]]+KV_Cp[best[0]]+KV_Cpp[best[0]]+KV_D[best[0]]+KV_E[best[0]]+KV_F[best[0]]+KJ[best[1]],light_t)


    print('Light Sequence '+str(k))

    if lc=='lambda':
        find_B_tol=regex.findall(str('('+LV_B[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_B_tol==[]:
            error=error+1
            find_B_tol=regex.findall(str('('+LV_B[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (L_B)')
        find_C_tol=regex.findall(str('('+LV_C[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_C_tol==[]:
            error=error+1
            find_C_tol=regex.findall(str('('+LV_C[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (L_C)')
        CDR1_start=final[0][0].find(find_B_tol[0])+len(find_B_tol[0])
        CDR1_end=final[0][0].find(find_C_tol[0])
        CDR1=final[0][1][CDR1_start:CDR1_end].replace('-','')
        print(CDR1)

        find_Cp_tol=regex.findall(str('('+LV_Cp[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_Cp_tol==[]:
            error=error+1
            find_Cp_tol=regex.findall(str('('+LV_Cp[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (L_Cp)')
        find_Cpp_tol=regex.findall(str('('+LV_Cpp[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_Cpp_tol==[]:
            error=error+1
            find_Cpp_tol=regex.findall(str('('+LV_Cpp[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (L_Cpp)')
        CDR2_start=final[0][0].find(find_Cp_tol[0])+len(find_Cp_tol[0])
        CDR2_end=final[0][0].find(find_Cpp_tol[0])
        CDR2=final[0][1][CDR2_start:CDR2_end].replace('-','')
        print(CDR2)

        find_F_tol=regex.findall(str('('+LV_F[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_F_tol==[]:
            error=error+1
            find_F_tol=regex.findall(str('('+LV_F[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (L_F)')
        find_J_tol=regex.findall(str('('+LJ[best[1]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_J_tol==[]:
            error=error+1
            find_J_tol=regex.findall(str('('+LJ[best[1]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (LJ)')
        CDR3_start=final[0][0].find(find_F_tol[0])+len(find_F_tol[0])
        CDR3_end=final[0][0].find(find_J_tol[0])
        CDR3=final[0][1][CDR3_start:CDR3_end].replace('-','')
        print(CDR3)


# In[278]:


# So here we leverage the fact that the alignments are the exact same length.
# Print the in betweens... CDR1 = Between B and C regions... CDR2= Between C' and C''... 
# CDR3= Between F domain and entire J gene usage.

# Need to have this regex ditty in there because of some slight disagreements for
# 
    if lc=='kappa':
        find_B_tol=regex.findall(str('('+KV_B[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_B_tol==[]:
            error=error+1
            find_B_tol=regex.findall(str('('+KV_B[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (K_B)')
        find_C_tol=regex.findall(str('('+KV_C[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_C_tol==[]:
            error=error+1
            find_C_tol=regex.findall(str('('+KV_C[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (K_C)')
        CDR1_start=final[0][0].find(find_B_tol[0])+len(find_B_tol[0])
        CDR1_end=final[0][0].find(find_C_tol[0])
        CDR1=final[0][1][CDR1_start:CDR1_end].replace('-','')
        print(CDR1)

        find_Cp_tol=regex.findall(str('('+KV_Cp[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_Cp_tol==[]:
            error=error+1
            find_Cp_tol=regex.findall(str('('+KV_Cp[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (K_Cp)')
        find_Cpp_tol=regex.findall(str('('+KV_Cpp[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_Cpp_tol==[]:
            error=error+1
            find_Cpp_tol=regex.findall(str('('+KV_Cpp[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (K_Cpp)')
        CDR2_start=final[0][0].find(find_Cp_tol[0])+len(find_Cp_tol[0])
        CDR2_end=final[0][0].find(find_Cpp_tol[0])
        CDR2=final[0][1][CDR2_start:CDR2_end].replace('-','')
        print(CDR2)

        find_F_tol=regex.findall(str('('+KV_F[best[0]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_F_tol==[]:
            error=error+1
            find_F_tol=regex.findall(str('('+KV_F[best[0]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (K_F)')
        find_J_tol=regex.findall(str('('+KJ[best[1]])+'){i<=0,s<=0,e<=0}',final[0][0])
        error=0
        while find_J_tol==[]:
            error=error+1
            find_J_tol=regex.findall(str('('+KJ[best[1]])+'){i<='+str(error)+',s<='+str(error)+',e<='+str(error)+'}',final[0][0])
            if error >= 4:
                print('High Mismatch with Conserved Regions, Check for Problems (KJ)')
        CDR3_start=final[0][0].find(find_F_tol[0])+len(find_F_tol[0])
        CDR3_end=final[0][0].find(find_J_tol[0])
        CDR3=final[0][1][CDR3_start:CDR3_end].replace('-','')
        print(CDR3)

