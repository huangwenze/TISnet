#-*- coding:utf-8 -*-
import sys
import numpy as np
import os
#from scipy.spatial.distance import pdist
#from scipy import stats



def read_table(tablefile):
    """
    path of shape file
    A       AT1G01010.1|132|232     GAGGATCGCCGCGACGTTGAAGT   1,1,1,0,0,1,0,0,0,0,1,1,1     1    1 0.42574257425742573     GAGGGAAGGAAACACUAGCCGCGACGUUGAAGU   .....((((((((((((((..((((...(((..((((((...))))))((.........)).........)))..)))).))..)))).)))).))))... -43.86  0       0.999939        0       20
    """
    sdict = {}
    f = open(tablefile, 'r')
    shape_list = []
    i = 0
    trx_id = ""
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent = line.split('\t')
        trx_id = sent[1]
        sdict[trx_id] = sent[2:]
    f.close()
    return sdict


def get_posi_sample(sdict):
    sdict2 = {}
    for ke in sdict:
        sent = sdict[ke]
        if (int(sent[3]) > 0) and (float(sent[9]) > 0.9):
            sdict2[ke] = sent[4:10]
            #average shape, sequence, structure, energy
    num1 = len(sdict2)
    return sdict2


def str_to_pair(str1):
    pdict = {}
    list1 = []
    for i in range(len(str1)):
        if str1[i] == '(':
            list1.append(i)
        elif str1[i] == ')':
            j = list1.pop()
            pdict[(j, i)] = i-j
        else:
            j = i
    return(pdict)


def define_motif(pdict):
    ll = sorted(pdict.items(), key=lambda pdict: pdict[1], reverse=True)
    pdict2 = {}
    for i in range(len(ll)):
        pr = ll[i][0]
        p1, p2 = pr[0], pr[1]
        flag = 0
        for ke in pdict2:
            sent1 = pdict2[ke][-1]
            if ke[0] < p1 and p2 < ke[1]:
                if sent1[0] < p1 and p2 < sent1[1]:
                    pdict2[ke].append((p1,p2))
                flag = 1
                break
        if flag == 0:
            pdict2[(p1,p2)] = [(p1,p2)]
    return(pdict2)


def motif_pattern(pdict1, seq1):
    pdict2 = {}
    number2 = 0
    mp1, mp2, maxdlen, maxslen, dchr1, dchr2, schr1 = 0, 0, 0, 0, "", "", ""
    for ke in pdict1:
        #dlen1, slen1 = len(sent1), 0
        d1, d2, s1 = "", "", ""
        sent1 = pdict1[ke]
        dlen1, slen1 = len(sent1), sent1[-1][1] - sent1[-1][0] + 1
        for i in range(len(sent1)):
            d1 = d1 + seq1[sent1[i][0]]
            d2 = d2 + seq1[sent1[i][1]]
            t1, t2 = sent1[i][0] + 1, sent1[i][1]
            s1 = seq1[t1:t2]
        pdict2[ke] = [str(dlen1), str(slen1), d1, d2, s1]
        if dlen1 > maxdlen:
            mp1, mp2, maxdlen, maxslen, dchr1, dchr2, schr1 = ke[0], ke[1], dlen1, slen1, d1, d2, s1
        number2 = number2 + 1
    dlist = [str(number2), str(mp1), str(mp2), str(maxdlen), str(maxslen), dchr1, dchr2, schr1]
    return(pdict2, dlist)


def take_motif_pattern(infile, outfile):
    tdict1 = read_table(infile)
    tdict2 = get_posi_sample(tdict1) #average shape, sequence, structure, energy
    fw1 = open(outfile, 'w')
    for ke in tdict2:
        str1 = tdict2[ke][2]
        pdict = str_to_pair(str1)
        pdict2 = define_motif(pdict)
        pdict3 = motif_pattern(pdict2, tdict2[ke][1])
        for pke in pdict3:
            outstr = ke + "\t" + "\t".join(tdict2[ke]) + "\t" + str(pke[0]) + "\t" + str(pke[1]) + "\t" + "\t".join(pdict3[pke]) + "\n"
            fw1.write(outstr)
    fw1.close()


def take_motif_pattern2(infile, outfile):
    tdict1 = read_table(infile)
    tdict2 = get_posi_sample(tdict1) #average shape, sequence, structure, energy
    fw1 = open(outfile, 'w')
    for ke in tdict2:
        str1 = tdict2[ke][2]
        pdict = str_to_pair(str1)
        pdict2 = define_motif(pdict)
        pdict3, dlist3 = motif_pattern(pdict2, tdict2[ke][1])
        #for pke in pdict3:
        outstr = ke + "\t" + "\t".join(tdict2[ke]) + "\t" + "\t".join(dlist3) + "\n"
        fw1.write(outstr)
    fw1.close()


def take_motif_pattern3(infile, outfile):
    tdict1 = read_table(infile)
    tdict2 = get_posi_sample(tdict1) #average shape, sequence, structure, energy
    sdict3 = {}
    fw1 = open(outfile, 'w')
    for ke in tdict2:
        str1 = tdict2[ke][2]
        pdict = str_to_pair(str1)
        pdict2 = define_motif(pdict)
        pdict3, dlist3 = motif_pattern(pdict2, tdict2[ke][1])
        #mp1, mp2, maxdlen, maxslen, dchr1, dchr2, schr1
        if int(dlist3[2]) >= 15:
            sdict3[ke] = [dlist3[5], dlist3[6], dlist3[7]]
        #for pke in pdict3:
        outstr = ke + "\t" + "\t".join(tdict2[ke]) + "\t" + "\t".join(dlist3) + "\n"
        fw1.write(outstr)
    fw1.close()
    return(sdict3)



def take_motif_pattern4(sdict, outfile1, outfile2):
    pat = ["AU", "UA", "UG", "GU", "CG", "GC"]
    kdict1 = kmer(sdict)
    kdict2 = kmer2(sdict)
    fw1 = open(outfile1, 'w')
    fw2 = open(outfile2, 'w')
    for ke in kdict1:
        outstr = ke + "\t" + str(kdict1[ke]) + "\n"
        fw1.write(outstr)
    outstr = "index\t" + "\t".join(pat) + "\n"
    fw2.write(outstr)
    for ke in kdict2:
        outstr = str(ke)
        for j in range(6):
            outstr = outstr + "\t" + str(kdict2[ke][j])
        outstr = outstr + "\n"
        fw2.write(outstr)
    fw1.close()
    fw2.close()


def kmer(sdict, length = 3):
    kdict = {}
    for ke in sdict: 
        schr = sdict[ke][2]
        for i in range(len(schr) - 3):
            ker = schr[i:(i+3)]
            if ker in kdict:
                kdict[ker] = kdict[ker] + 1
            else:
                kdict[ker] = 1
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

def gene_str(data_file, tmpfile1, tmpfile2, outfile):
    """
    dict, shape_dict       -- Two dictionary of float values
    A   AT1G04440.1   GGGGTTTGAGCCCCGGAGGCAGAGCGGCTGCCATGGCCAA NULL,NULL,0.144,0.216    0.284414195    1
    Return -1 if failed
    """
    dat1 = {}
    f = open(data_file, 'r')
    i = 0
    for line in f:
        i = i + 1
        fw1 = open(tmpfile1, 'w')
        fw2 = open(tmpfile2, 'w')
        line = line.strip('\n')
        sent = line.split('\t')
        sent1 = sent[3].split(',')
        outstr = ">" + sent[1] + "\n" + sent[2] + "\n"
        fw1.write(outstr)
        outstr = ""
        for j in range(len(sent1)):
            #if sent1[j] == "NULL":
            #    outstr = str(j+1) + "\t" + str(-1) + "\n"
            #    fw2.write(outstr)
            #else:
            outstr = str(j+1) + "\t" + str(float(sent1[j])*2) + "\n"
            fw2.write(outstr)
        fw1.close()
        fw2.close()
        cmd = "RNAfold --noPS --shape=" + tmpfile2 + " < " + tmpfile1 + " >> " + outfile
        os.system(cmd)
    f.close()

def str_file(infile1, infile2, outfile3):
    f1 = open(infile1, 'r')
    f2 = open(infile2, 'r')
    fw3 = open(outfile3, 'w')
    i = 0
    dict1 = {}
    trx_id = ""
    sequence = ""
    structure = ""
    for line in f2:
        i = i + 1
        line = line.strip('\n')
        if i % 3 == 1:
            trx_id = line[1:]
        elif i % 3 == 2:
            sequence = line
        else:
            sent = line.split(" ")
            structure = sent[0]
            dict1[trx_id] = [sequence, structure]
            #print(trx_id)
    for line in f1:
        line = line.strip('\n')
        sent = line.split("\t")
        if sent[1] in dict1:
            outstr = line + "\t" + dict1[sent[1]][0] + "\t" + dict1[sent[1]][1] + "\n"
            fw3.write(outstr)
        else:
            print(sent[1])
    f1.close()
    f2.close()
    fw3.close()


if __name__ == '__main__':
    dat_file = sys.argv[1] 
    out_file1 = sys.argv[2]
    out_file2 = sys.argv[3]
    out_file3 = sys.argv[4]
    tdict = take_motif_pattern3(dat_file, out_file1)
    take_motif_pattern4(tdict, out_file2, out_file3)


