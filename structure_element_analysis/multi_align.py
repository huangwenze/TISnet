#-*- coding:utf-8 -*-
import sys
import numpy as np
import os
import pdb
from itertools import product
from collections import deque
import statistics as stat
#from scipy.spatial.distance import pdist
#from scipy import stats


def read_table(tablefile1):
    """
    path of shape file
    AT1G01010.1|132|232 0.42574257425742573 GAGGAUCAAGUUGGGUUUGGGUUCCGUCCGAACGACGAGGAGCUCGUUGGUCACUAUCUCCGUAACAAAAUCGAAGGAAACACUAGCCGCGACGUUGAAGU   .....((((((((((((((..((((...(((..((((((...))))))((.........)).........)))..)))).))..)))).)))).))))...   -43.86  0   0.999939    1   5   97  27  4   UCAAGUUGGGUUUGUUCCCGAGACGAG AGUUCAGCCCGAACAAGGGCUUUGCUC GAG
    """
    sdict1 = {}
    sdict2 = {}
    sdict3 = {}
    f = open(tablefile1, 'r')
    shape_list = []
    i = 0
    trx_id = ""
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent1 = line.split('\t')
        trx_id = sent1[0]
        sdict1[sent1[1]] = 0
        sdict2[sent1[1]] = sent1[2]
        sdict3[sent1[3]] = 0
    f.close()
    #sdict = {}
    n1 = len(sdict1)
    print(str(n1))
    return(sdict1, sdict2, sdict3)


def kmer(sdict, index, length1 = 3):
    kdict = {}
    for ke in sdict: 
        schr = sdict[ke][index]
        for i in range(len(schr) - length1):
            ker = schr[i:(i+length1)]
            if ker in kdict:
                kdict[ker] = kdict[ker] + 1
            else:
                kdict[ker] = 1
    return(kdict)


def seqscore(schr, kmer1, len1):
    kdict = {}
    score1 = 0
    for i in range(len(schr) - len1):
        ker = schr[i:(i+len1)]
        if ker in kmer1:
            score1 = score1 + kmer1[ker]
    return(score1) 


def kscore(sdict, kmer1, len1, kmer2, len2, kmer3, len3):
    kdict = {}
    for ke in sdict: 
        schr1, schr2, schr3 = sdict[ke][0], sdict[ke][1], sdict[ke][2]
        score1 = seqscore(schr1, kmer1, len1)
        score2 = seqscore(schr2, kmer2, len2)
        score3 = seqscore(schr3, kmer3, len3)
        sid = schr1 + "|" + schr2 + "|" + schr3
        kdict[sid] = score1 + score2 + score3
    return(kdict) 


def kmer2(sdict, length = 3):
    kdict = {}
    pat = {"AU":0, "UA":1, "UG":2, "GU":3, "CG":4, "GC":5}
    for ke in sdict: 
        dchr1 = sdict[ke][0][::-1]
        dchr2 = sdict[ke][1][::-1]
        for i in range(len(dchr1)):
            ker = dchr1[i] + dchr2[i]
            if i in kdict:
                kdict[i][pat[ker]] = kdict[i][pat[ker]] + 1
            else:
                kdict[i] = [0, 0, 0, 0, 0, 0]
                kdict[i][pat[ker]] = kdict[i][pat[ker]] + 1
    return(kdict)


def needleman_wunsch(x, y):
    """Run the Needleman-Wunsch algorithm on two sequences.

    x, y -- sequences.

    Code based on pseudocode in Section 3 of:

    Naveed, Tahir; Siddiqui, Imitaz Saeed; Ahmed, Shaftab.
    "Parallel Needleman-Wunsch Algorithm for Grid." n.d.
    https://upload.wikimedia.org/wikipedia/en/c/c4/ParallelNeedlemanAlgorithm.pdf
    """
    N, M = len(x), len(y)
    s = lambda a, b: int(a == b)

    DIAG = -1, -1
    LEFT = -1, 0
    UP = 0, -1

    # Create tables F and Ptr
    F = {}
    Ptr = {}

    F[-1, -1] = 0
    for i in range(N):
        F[i, -1] = -i
    for j in range(M):
        F[-1, j] = -j

    option_Ptr = DIAG, LEFT, UP
    for i, j in product(range(N), range(M)):
        option_F = (
            F[i - 1, j - 1] + s(x[i], y[j]),
            F[i - 1, j] - 1,
            F[i, j - 1] - 1,
        )
        F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))

    # Work backwards from (N - 1, M - 1) to (0, 0)
    # to find the best alignment.
    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        if direction == DIAG:
            element = i, j
        elif direction == LEFT:
            element = i, None
        elif direction == UP:
            element = None, j
        alignment.appendleft(element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
    while j >= 0:
        alignment.appendleft((None, j))
        j -= 1

    return list(alignment)


def cal_identity(x, y):
    alist = needleman_wunsch(x, y)
    n2, n3 = 0, len(alist)
    for i, j in alist:
        if (i is not None) and (j is not None) and (x[i] == y[j]) and (x[i] != "-"):
            n2 = n2 + 1 
    return(n2/n3)


def compare_dict(ke1, dict1):
    list1 = []
    dict2 = dict1
    flag = 0
    for ke2 in dict1:
        rate = cal_identity(ke1, ke2)
        dict2[ke2].append(rate)
        list1.append(rate)
    dict2[ke1] = list1
    return(dict2)


def cal_similar(dict1):
    ll = list(dict1.keys())
    ke1 = ll[0]
    dict2 = {}
    dict2[ke1] = []
    ll = ll[1:]
    for i in range(len(ll)):
        dict2 = compare_dict(ll[i], dict2)
    return(dict2)


def change_seq(chr1, chr2, alist):
    chr01, chr02 = "", ""
    for i,j in alist:
        if (i is not None):
            chr01 = chr01 + chr1[i]
        else:
            chr01 = chr01 + "-"
        if (j is not None):
            chr02 = chr02 + chr2[j]
        else:
            chr02 = chr02 + "-"      
    return(chr01, chr02)


def malign(dict1):
    dict2 = {}
    for ke in dict1:
        #pdb.set_trace()
        dict2[ke] = stat.mean(dict1[ke])
    sl = sorted(dict2.items(), key=lambda d: d[1], reverse = True)
    ke1 = sl[0][0]
    clist = [ke1]
    clist2 = []
    ddist = {}
    index1 = 1
    ddist[index1] = [ke1]
    for ke2, num in sl[1:]:
        clist = ddist[index1]
        ke1 = clist[0]
        index1 = index1 + 1
        ddist[index1] = [kk for kk in clist]
        alist = needleman_wunsch(ke1, ke2)
        ke3 = ke2
        #pdb.set_trace()
        for k in range(len(ddist[index1])):
            chr1, chr2 = change_seq(ddist[index1][k], ke2, alist)
            ddist[index1][k] = chr1
        ddist[index1].append(chr2)
        #clist = clist2 
        #ddist[index1].append(ke3)
    return(ddist[index1])


def talign(llist1, dict1):
    llist2 = []
    for ke1 in llist1:
        ke2 = ke1.replace("-", "")
        if ke2 in dict1:
            ke3 = dict1[ke2]
            ke4 = ""
            j = 0
            for i in range(len(ke1)):
                if ke1[i] == "-":
                    ke4 = ke4 + "-"
                else:
                    ke4 = ke4 + ke3[j]
                    j = j + 1
            llist2.append(ke4)
        else:
            print(ke1 + "\t" + ke2)
    #llist2.append(ke4)
    return(llist2)


def align_mat(list1):
    mat1 = [[0, 0, 0, 0] for i in range(len(list1[0]))]
    for ke in list1:
        for i in range(len(ke)):
            if ke[i] == "A":
                mat1[i][0] = mat1[i][0] + 1
            elif ke[i] == "C":
                mat1[i][1] = mat1[i][1] + 1
            elif ke[i] == "G":
                mat1[i][2] = mat1[i][2] + 1
            elif ke[i] == "U":
                mat1[i][3] = mat1[i][3] + 1
    return(mat1)


if __name__ == '__main__':
    dat_file1 = sys.argv[1] 
    out_file1 = sys.argv[2]
    out_file2 = sys.argv[3]
    #out_file3 = sys.argv[4]
    tdict1, tdict2, tdict3= read_table(dat_file1)
    tdict01 = cal_similar(tdict1)
    #tdict02 = cal_similar(tdict2)
    tdict03 = cal_similar(tdict3)
    clist1 = malign(tdict01)
    #clist2 = malign(tdict02)
    clist2 = talign(clist1, tdict2)
    clist3 = malign(tdict03)
    mat1 = align_mat(clist1)
    mat2 = align_mat(clist2)
    mat3 = align_mat(clist3)
    fw1 = open(out_file1, 'w')
    fw2 = open(out_file2, 'w')
    outstr1 = ">p1\n" + "\n".join(clist1) + "\n" + ">p2\n" + "\n".join(clist2) + "\n" + ">loop\n" + "\n".join(clist3) + "\n"
    fw1.write(outstr1)
    for i in range(len(mat1)):
        outstr2 = "1\t" + str(i)
        for j in range(4):
            outstr2 = outstr2 + "\t" + str(mat1[i][j])
        outstr2 = outstr2 + "\n"
        fw2.write(outstr2)
    for i in range(len(mat2)):
        outstr2 = "2\t" + str(i)
        for j in range(4):
            outstr2 = outstr2 + "\t" + str(mat2[i][j])
        outstr2 = outstr2 + "\n"
        fw2.write(outstr2)
    for i in range(len(mat3)):
        outstr2 = "3\t" + str(i)
        for j in range(4):
            outstr2 = outstr2 + "\t" + str(mat3[i][j])
        outstr2 = outstr2 + "\n"
        fw2.write(outstr2)
    fw1.close()
    fw2.close()


