#-*- coding:utf-8 -*-
import sys
import numpy as np
import os
from itertools import product
from collections import deque
#from scipy.spatial.distance import pdist
#from scipy import stats


def read_table(tablefile1, tablefile2):
    """
    path of shape file
    AT1G01010.1|132|232 0.42574257425742573 GAGGAUCAAGUUGGGUUUGGGUUCCGUCCGAACGACGAGGAGCUCGUUGGUCACUAUCUCCGUAACAAAAUCGAAGGAAACACUAGCCGCGACGUUGAAGU   .....((((((((((((((..((((...(((..((((((...))))))((.........)).........)))..)))).))..)))).)))).))))...   -43.86  0   0.999939    1   5   97  27  4   UCAAGUUGGGUUUGUUCCCGAGACGAG AGUUCAGCCCGAACAAGGGCUUUGCUC GAG
    """
    sdict = {}
    f = open(tablefile1, 'r')
    shape_list = []
    i = 0
    trx_id = ""
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent1 = line.split('\t')
        sent2 = sent1[0].split('|')
        trx_id = sent2[0] + "|mAUG"
        if (float(sent1[6]) > 0.90) and (float(sent1[10]) >= 15):
            sdict[trx_id] = sent1[12:] #p1, p2, loopseq
    f.close()
    #sdict = {}
    f = open(tablefile2, 'r')
    shape_list = []
    i = 0
    trx_id = ""
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent1 = line.split('\t')
        sent2 = sent1[0].split('|')
        trx_id = sent2[0] + "|uAUG"
        if (float(sent1[6]) > 0.90) and (float(sent1[10]) >= 15):
            sdict[trx_id] = sent1[12:] #p1, p2, loopseq
    f.close()
    n1 = len(sdict)
    print(str(n1))
    return sdict


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


def cal_identity2(x, y):
    alist = needleman_wunsch(x, y)
    n1, n2, n3 = 0, 0, 0
    for i, j in alist:
        if (i is not None) and (j is not None) and (x[i] == y[j]) and (x[i] != "-"):
            n2 = n2 + 1
    for i in range(len(x)):
        if x[i] != "-":
            n1 = n1 + 1
    for i in range(len(y)):
        if y[i] != "-":
            n3 = n3 + 1    
    return(n1, n2, n3)


def cal_identity(x, y):
    alist = needleman_wunsch(x, y)
    n2, n3 = 0, len(alist)
    for i, j in alist:
        if (i is not None) and (j is not None) and (x[i] == y[j]) and (x[i] != "-"):
            n2 = n2 + 1 
    return(n2, n3)


def compare_str(ke1, ke2):
    sent1 = ke1.split("|")
    sent2 = ke2.split("|")
    n02, n03 = cal_identity(sent1[0], sent2[0])
    n12, n13 = cal_identity(sent1[1], sent2[1])
    n22, n23 = cal_identity(sent1[2], sent2[2])
    return((1.0 * n02 / n03 + 1.0 * n12 / n13 + 2.0 * n22 / n23)/4)


def compare_list(ke1, list1, rate1):
    flag = 0
    for ke2 in list1:
        if compare_str(ke1, ke2) > rate1:
            flag = 1
            break
    return(flag)



def seq_cluster(sdict1):
    #sl = sorted(sdict1.items(), key=lambda d: d[1], reverse = True)
    ddict1 = {}
    ddict2 = {}
    print(str(len(sdict1)))
    i = 0
    sl = list(sdict1.keys())
    for i in range(len(sl)):
        ddict1[sl[i]] = ["0" for j in range(len(sl))]
        ddict2[sl[i]] = sdict1[sl[i]][0] + "|" + sdict1[sl[i]][1] + "|" + sdict1[sl[i]][2]
    for i in range(len(sl)):
        if i % 100 == 0:
            print(str(i))
        for j in range(i + 1, len(sl)):
            ke1 = sdict1[sl[i]][0] + "|" + sdict1[sl[i]][1] + "|" + sdict1[sl[i]][2] 
            ke2 = sdict1[sl[j]][0] + "|" + sdict1[sl[j]][1] + "|" + sdict1[sl[j]][2] 
            ddict1[sl[i]][j] = str(round(compare_str(ke1, ke2),4))
    return(ddict1, ddict2)



if __name__ == '__main__':
    dat_file1 = sys.argv[1]
    dat_file2 = sys.argv[2] 
    out_file1 = sys.argv[3]
    out_file2 = sys.argv[4]
    #out_file3 = sys.argv[4]
    tdict = read_table(dat_file1, dat_file2)
    ddict1, ddict2 = seq_cluster(tdict)
    fw1 = open(out_file1, 'w')
    fw2 = open(out_file2, 'w')
    dl = list(ddict1.keys())
    for pke in dl:
        outstr1 = pke + "\t" + ddict2[pke] + "\n"
        outstr2 = pke + "\t" + "\t".join(ddict1[pke]) + "\n"
        fw1.write(outstr1)
        fw2.write(outstr2)
    fw1.close()
    fw2.close()


