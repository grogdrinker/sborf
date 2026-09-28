#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parser.py
#
#  Copyright 2019 u0127326 <u0127326@gbw-l-l0039>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
import pickle, os

SCRIPT_PATH = "/".join(os.path.abspath(__file__).split("/")[:-1])


def leggifasta(database):  # legge un file fasta e lo converte in un dizionario
    f = open(database)
    uniprot = f.readlines()
    f.close()
    return leggifasta_string("".join(uniprot))


def leggifasta_string(text):  # parses fasta content from a string (e.g. copy-pasted)
    dizio = {}
    uniprotid = None
    # lines are trimmed first (spaces, tabs, non-breaking spaces), so a header is recognised even when it is
    # indented, and a byte-order mark at the very start of the text is ignored
    for i in text.lstrip('﻿').splitlines():
        i = i.strip()
        if i == '':
            continue
        if i[0] == '>':
            header = i.lstrip('>').strip()
            if '|' in header:
                uniprotid = header.split('|')[1].strip()
            else:
                uniprotid = header.split()[0] if header.split() else ''
            if not uniprotid:
                uniprotid = 'sequence_' + str(len(dizio) + 1)
            dizio[uniprotid] = ''
        else:
            if uniprotid is None:
                continue
            dizio[uniprotid] = dizio[uniprotid] + i
    return dizio


def parse_featurea(fil=SCRIPT_PATH + '/dataset/db1_1072017.csv'):
    header = True
    cont = 0
    diz = {}
    '''
    hashing={'TKG':[1,0,0,0,0,0,0,0,0], 
             'KG': [0,1,0,0,0,0,0,0,0], 
             'G':  [0,0,1,0,0,0,0,0,0], 
             'K':  [0,0,0,1,0,0,0,0,0], 
             'N':  [0,0,0,0,1,0,0,0,0], 
             'S':  [0,0,0,0,0,1,0,0,0],
             'TK': [0,0,0,0,0,0,1,0,0], 
             'TG': [0,0,0,0,0,0,0,1,0], 
             'T':  [0,0,0,0,0,0,0,0,1]
             }
    '''
    hashing = {'G': [1, 0, 0],
               'N': [0, 0, 0],
               'T': [0, 1, 0],
               'TG': [1, 1, 0],
               'K': [0, 0, 1],
               'KG': [1, 0, 1],
               'TK': [0, 1, 1],
               'TKG': [1, 1, 1],
               'S': [-1, -1, -1]
               }
    for lin in open(fil).readlines():
        if header:
            h = lin.strip().replace('"', '').split('\t')
            header = False
        else:
            a = lin.strip().replace('"', '').split('\t')
            fea = []
            tmp = {}

            for i in range(len(h)):
                if h[i] == 'uid':
                    name = a[i]

                elif h[i] == 'C3':
                    chap = a[i]
                else:

                    # print a[i],h[i]

                    tmp[h[i]] = a[i]
                # continue
            if chap == 'S':
                pass
            # continue
            diz[name] = {}
            diz[name]["C3"] = hashing[chap]
            diz[name].update(tmp)
    f = open('aasa', 'w')
    for i in diz.keys():
        f.write(i + '\n')

    return diz


def parse_dna_rash(fil='dataset/rna_full.tab'):
    diz = {}
    for i in open(fil).readlines():
        if i == '\n':
            continue
        else:
            a = i.split()
            if a[1] == 'XXX':
                continue
            else:
                cod = [a[1][c:c + 3] for c in range(0, len(a[1]), 3)]
                assert len(cod[-1]) == 3
                diz[a[0]] = cod
    return diz

def codon_to_seq(s):
    return "".join([codon2aa[i.replace("T","U")] for i in s])

codon2aa = {"AAA": "K", "AAC": "N", "AAG": "K", "AAU": "N",
            "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T",
            "AGA": "R", "AGC": "S", "AGG": "R", "AGU": "S",
            "AUA": "I", "AUC": "I", "AUG": "M", "AUU": "I",

            "CAA": "Q", "CAC": "H", "CAG": "Q", "CAU": "H",
            "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
            "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R",
            "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L",

            "GAA": "E", "GAC": "D", "GAG": "E", "GAU": "D",
            "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A",
            "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G",
            "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",

            "UAA": "", "UAC": "Y", "UAG": "", "UAU": "T",
            "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S",
            "UGA": "", "UGC": "C", "UGG": "W", "UGU": "C",
            "UUA": "L", "UUC": "F", "UUG": "L", "UUU": "F"}

allowed_codons = {}
aa2codon = {}
for i in codon2aa.keys():
    if not codon2aa[i] in aa2codon:
        aa2codon[codon2aa[i]] = []
    aa2codon[codon2aa[i]]+= [i.replace("U","T")]

    if codon2aa[i]!="":
        allowed_codons[i.replace("U","T")] = True
def parse_raw(fil='dataset/K_full.csv'):
    f = open(fil, 'r')
    cont = 0
    paper = []
    diz = {}
    for i in f.readlines():
        a = i.split('\t')
        if cont == 0:
            for k in a:
                if k == "":
                    break
                else:
                    paper += [k]
        else:
            diz[a[0]] = {}
            # print paper

            for k in range(len(paper[1:])):
                diz[a[0]]['agg_' + paper[k + 1].strip()] = int(a[k + 1])
        cont += 1
    return diz


def parse_raw_reshmi(seqfile='dataset/seqs20cdhit.fasta', dnafile='dataset/rna/rna_full.tab'):
    # seqfile='dataset/seqsFull.fasta'

    g = parse_raw(fil='dataset/G_full.csv')

    t = parse_raw(fil='dataset/T_full.csv')
    k = parse_raw(fil='dataset/K_full.csv')
    df = {}

    dna = parse_dna_rash()
    seq = leggifasta(seqfile)
    dafare = seq.keys()
    final = {}
    saltato = 0
    for i in dafare:
        if i in dna and i in g and i in t and i in k:
            final[i] = {}
            final[i]['G'] = g[i]
            final[i]['T'] = t[i]
            final[i]['K'] = k[i]
            final[i]['seq'] = seq[i]
            final[i]['dna'] = dna[i][1:]
        else:
            saltato += 1
    print('\ttotal of ', len(final.keys()), 'matches. Skipped', saltato, 'proteins')
    return final


def parser_and_dna():
    from Bio import pairwise2

    from Bio.SubsMat import MatrixInfo as matlist

    matrix = matlist.blosum62

    gap_open = -10

    gap_extend = -0.5
    diz = final_parser()
    dna = leggifasta('dataset/rna/dna.fasta')
    fin = {}
    cont = 0
    for id in dna.keys():
        cont += 1
        if not id in diz:
            continue
        prot_seq = diz[id]['seq']
        s = dna[id].replace('T', 'U')
        max_score = -999
        for rf in range(3):
            n = 3
            cods = [s[i:i + n] for i in range(rf, len(s), n)]
            p = ''
            cod_p = []
            prots = []
            dnas = []
            # print cods

            for k in cods:
                if len(k) < 3:
                    continue
                if not k in codon2aa:
                    continue
                c = codon2aa[k]

                if c == '_':
                    if p == '':
                        continue
                    else:
                        prots += [p]
                        p = ''
                        dnas += [cod_p]
                        cod_p = []
                else:
                    p += codon2aa[k]
                    cod_p += [k]

            for pn in range(len(prots)):
                try:
                    dna_a, pro_a, score, begin, end = \
                    pairwise2.align.globalds(prots[pn], prot_seq, matrix, gap_open, gap_extend)[0]
                except:
                    continue
                if score > max_score:
                    max_score = score
                    dfin = dna_a[begin:end]
                    pfin = pro_a[begin:end]

                    conta = 0
                    cods_a = []
                    segmentodna = dnas[pn]
                    for l in range(len(dna_a)):
                        if dna_a[l] != '-':
                            cods_a += [segmentodna[conta]]
                            conta += 1
                        else:
                            cods_a += ['-']
                    dn = cods_a[begin:end]
                    # print len(dn),len(dfin),len(pfin)
                    lkd = []
                    lkp = []
                    lkc = []
                    b = 0
                    for i in range(len(dfin)):
                        if dfin[i] == '-' or pfin[i] == '-':
                            lkd += [dfin[b:i]]
                            lkp += [pfin[b:i]]
                            lkc += [dn[b:i]]
                            b = i + 1
                    if b != len(dfin):
                        lkd += [dfin[b:i]]
                        lkp += [pfin[b:i]]
                        lkc += [dn[b:i]]
                    maxlen = 0
                    candp = []
                    candd = []
                    candf = []
                    for i in range(len(lkd)):
                        if len(lkd[i]) > maxlen:
                            maxlen = len(lkd[i])
                            candd = lkc[i]
                            candp = lkp[i]
                            candf = lkc[i]
                    if len(prot_seq) * 0.8 < len(candp) and score > 0:
                        dfin = candp
                        pfin = candp
                        dn = candf
        if not '-' in dfin and not '-' in pfin:
            s = ''
            for i in dn:
                s += codon2aa[i]

            # assert s==pfin
            fin[id] = diz[id]
            fin[id]['seq'] = s
            fin[id]['dna'] = dn



    pickle.dump(fin, open('marshalled/diz_dna.m', 'w'))


def parse_tests_BurgessBrown2008(fil="dataset/tests_BurgessBrown2008.csv"):
    diz = {}
    for l in open(fil, "r").readlines()[1:]:
        a = l.split("\t")
        name = a[0]
        native = a[3].strip()[:-3]
        mut = a[4].strip()[:-3]
        diz[name] = [native, mut]
    return diz


def parsePaxDB(fol="dataset/paxdb/human/", fileFasta="dataset/paxdb/human_abundance.fasta", mapp="dataset/paxdb/full_uniprot_2_paxdb.04.2015.tsv"):
    diz = {}
    seqs = leggifasta(fileFasta)  # per crearlo guarda dataset/paxdb/protein2dna.py
    for fil in os.listdir(fol):
        for l in open(fol + fil, "r").readlines():
            if l[0] == "#":
                continue
            if l.strip() == '':
                continue
            a = l.split()

            if len(a)==3:
                paxID = a[1].split(".")[1]
                abund = float(a[2])
            else:
                paxID = a[0].split(".")[1]
                abund = float(a[1])


            diz[paxID] = {}
            diz[paxID]["abn_new"] = abund

    mapping = {}
    for l in open(mapp).readlines():
        a = l.split()
        uniprot = a[1].split("|")[0]
        paxid = a[2]
        mapping[paxid] = uniprot

    fin = {}
    for id in diz.keys():
        if id in mapping:
            if mapping[id] in seqs:
                fin[mapping[id]] = diz[id]
                a = [seqs[mapping[id]][i:i + 3] for i in range(0, len(seqs[mapping[id]]), 3)]
                if len(a[-1]) != 3:
                    del fin[mapping[id]]
                    continue
                aaSeq = ""
                for i,k in enumerate(a):
                    if codon2aa[k.replace("T","U")] == "_" and i==len(a)-1:
                        a=a[:-1] # removing last stop codon
                        continue
                    aaSeq += codon2aa[k.replace("T","U")]
                if not "_" in aaSeq:
                    fin[mapping[id]]["aa"] = aaSeq
                    fin[mapping[id]]["dna"] = a
                else:
                    del fin[mapping[id]]


    return fin

def parsePaxDBNewformat(fol="dataset/paxdb/human/", fileFasta="dataset/paxdb/human_abundance.fasta"):
    diz = {}
    seqs = leggifasta(fileFasta)  # per crearlo guarda dataset/paxdb/protein2dna.py
    for fil in os.listdir(fol):
        for l in open(fol + fil, "r").readlines():
            if l[0] == "#":
                continue
            if l.strip() == '':
                continue
            a = l.split()

            if len(a)==3:
                paxID = a[1].split(".")[1]
                abund = float(a[2])
            else:
                paxID = a[0].split(".")[1]
                abund = float(a[1])


            diz[paxID] = {}
            diz[paxID]["abn_new"] = abund

    if "." in list(seqs.items())[0][0]:
        seqTMP={}
        for i in seqs.keys():
            seqTMP[i.split(".")[1]] = seqs[i]
        seqs = seqTMP
    fin = {}
    for id in diz.keys():
        aaSeq = ""
        if not id in seqs:
            continue
        assert len(seqs[id])%3==0
        a = [seqs[id][i:i+3] for i in range(0,len(seqs[id]),3)]
        seqs[id] = a
        for i,k in enumerate(seqs[id]):
            if codon2aa[k.replace("T","U")] == "_" and i==len(seqs[id])-1:
                seqs[id]=seqs[id][:-1] # removing last stop codon
                continue
            aaSeq += codon2aa[k.replace("T","U")]
        fin[id] = {}
        if not "_" in aaSeq:
            fin[id]["aa"] = aaSeq
            fin[id]["dna"] = seqs[id]
            fin[id]["abn_new"] = float(diz[id]["abn_new"])
        else:
            del  fin[id]
    return fin


def parse_featurea(fil='dataset/db1_1072017.csv'):
    header = True
    cont = 0
    diz = {}

    hashing = {'G': [1, 0, 0],
               'N': [0, 0, 0],
               'T': [0, 1, 0],
               'TG': [1, 1, 0],
               'K': [0, 0, 1],
               'KG': [1, 0, 1],
               'TK': [0, 1, 1],
               'TKG': [1, 1, 1],
               'S': [-1, -1, -1]
               }
    for lin in open(fil).readlines():
        if header:
            h = lin.strip().replace('"', '').split('\t')
            header = False
        else:
            a = lin.strip().replace('"', '').split('\t')
            fea = []
            tmp = {}

            for i in range(len(h)):
                if h[i] == 'uid':
                    name = a[i]

                elif h[i] == 'C3':

                    chap = a[i]
                else:

                    # print a[i],h[i]

                    tmp[h[i]] = a[i]

            # continue
            if chap=="S":
                continue
            diz[name] = {}

            if chap=="N":
                diz[name]["C3"] = 0
            else:
                diz[name]["C3"] = 1

            """
            if "K" in chap:#=="N":
                diz[name]["C3"] = 1
            else:
                diz[name]["C3"] = 0
            """
            #diz[name]["C3"] = hashing[chap]
            diz[name].update(tmp)
    f = open('aasa', 'w')
    for i in diz.keys():
        f.write(i + '\n')

    return diz

def final_parser(feafile='dataset/db1_1072017.csv', seqfile='dataset/seqs20cdhit.fasta'):
    fea=parse_featurea(feafile)
    dna=parse_dna_rash()
    seq=leggifasta(seqfile)
    dafare=dna.keys()
    final={}
    saltato=0
    for i in dafare:
        if i in fea and i in dna:
            final[i]=fea[i]
            final[i]['seq']=codon_to_seq(dna[i][1:])
            final[i]['dna']=dna[i][1:]
        else:
            saltato+=1
    print('\ttotal of ', len(final.keys()), 'matches. Skipped', saltato, 'proteins')
    return final