#!/opt/conda/bin/python3
# -*- coding: utf-8 -*-

#タマネギ色のONTデータをデュアルインデックスで振り分ける


from Bio import SeqIO
from Bio.Seq import Seq
import os
import gzip
import re
import Levenshtein
import argparse
import sys


#inputを指定
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
args = parser.parse_args()


print(args,input)

if args.input == "None":
    print("inputのfastq.gzファイルを指定してください")
    sys.exit()

in_fq = args.input

file_name = in_fq.split('.')[0]
print(file_name)

#recs = SeqIO.parse(open(in_fq, 'rt', encoding='utf-8'), 'fastq')
recs = SeqIO.parse(gzip.open(in_fq, 'rt'), 'fastq')
rec1_F = []
rec2_F = []
rec3_F = []
rec4_F = []
rec5_F = []
rec6_F = []
rec7_F = []
rec8_F = []
rec9_F = []
rec10_F = []
rec11_F = []
rec12_F = []
rec13_F = []
rec14_F = []
rec15_F = []
rec16_F = []
not_found = []

for rec in recs:
    sequence = rec.seq

    #sequenceからtail部分を検索
    #if  re.search("ACTTGCCTGTCGCTCTATCTTC",str(sequence)):
    #    tail_seq = re.search("ACTTGCCTGTCGCTCTATCTTC",str(sequence)).span()
    #elif re.search("TTTCTGTTGGTGCTGATATTGC",str(sequence)):
    #    tail_seq = re.search("TTTCTGTTGGTGCTGATATTGC",str(sequence)).span()
    #else:
    #    not_found.append(rec)
    #    continue
    #tail部分の検索を8塩基で行う
    if  re.search("CTATCTTC",str(sequence)):
        tail_seq = re.search("CTATCTTC",str(sequence)).span()
    elif re.search("GATATTGC",str(sequence)):
        tail_seq = re.search("GATATTGC",str(sequence)).span()
    else:
        not_found.append(rec)
        continue

    #tail部分の次の塩基から8塩基がindex配列
    #tail部分が末端の150塩基以内に検出された場合に採用
    if tail_seq[1] < 151:
        index_start = tail_seq[1]
        index_end = index_start + 8
    else:
        not_found.append(rec)
        continue

    index = sequence[index_start:index_end]
    index1_F = "ATTACTCG"
    index2_F = "TCCGGAGA"
    index3_F = "CGCTCATT"
    index4_F = "GAGATTCC"
    index5_F = "AACCTACG"
    index6_F = "GAATTCGT"
    index7_F = "CTGATGAG"
    index8_F = "TAATGCGC"
    index9_F = "CGGCTATG"
    index10_F = "TGCTTGCT"
    index11_F = "TCTCGCGC"
    index12_F = "AGCGATAG"
    index13_F = "ACAGCGTA"
    index14_F = "GTACGTAC"
    index15_F = "GCTGACTA"
    index16_F = "CTAGAGAT"

    #indexを16個に振り分ける.2塩基までのミスを許容する
    if Levenshtein.distance(str(index), index1_F) < 3:
        rec1_F.append(rec.id)

    elif Levenshtein.distance(str(index), index2_F) < 3:
        rec2_F.append(rec.id)

    elif Levenshtein.distance(str(index), index3_F) < 3:
        rec3_F.append(rec.id)

    elif Levenshtein.distance(str(index), index4_F) < 3:
        rec4_F.append(rec.id)

    elif Levenshtein.distance(str(index), index5_F) < 3:
        rec5_F.append(rec.id)

    elif Levenshtein.distance(str(index), index6_F) < 3:
        rec6_F.append(rec.id)

    elif Levenshtein.distance(str(index), index7_F) < 3:
        rec7_F.append(rec.id)

    elif Levenshtein.distance(str(index), index8_F) < 3:
        rec8_F.append(rec.id)

    elif Levenshtein.distance(str(index), index9_F) < 3:
        rec9_F.append(rec.id)

    elif Levenshtein.distance(str(index), index10_F) < 3:
        rec10_F.append(rec.id)

    elif Levenshtein.distance(str(index), index11_F) < 3:
        rec11_F.append(rec.id)

    elif Levenshtein.distance(str(index), index12_F) < 3:
        rec12_F.append(rec.id)

    elif Levenshtein.distance(str(index), index13_F) < 3:
        rec13_F.append(rec.id)

    elif Levenshtein.distance(str(index), index14_F) < 3:
        rec14_F.append(rec.id)

    elif Levenshtein.distance(str(index), index15_F) < 3:
        rec15_F.append(rec.id)

    elif Levenshtein.distance(str(index), index16_F) < 3:
        rec16_F.append(rec.id)

    else:
        not_found.append(rec.id)

#反対向きから調べる
recs = SeqIO.parse(gzip.open(in_fq, 'rt'), 'fastq')
rec1_R = []
rec2_R = []
rec3_R = []
rec4_R = []
rec5_R = []
rec6_R = []
rec7_R = []
rec8_R = []
rec9_R = []
rec10_R = []
rec11_R = []
rec12_R = []
rec13_R = []
rec14_R = []
rec15_R = []
rec16_R = []
unknown = []

for rec in recs:
    sequence = rec.seq
    sequence = sequence.reverse_complement()

    #sequenceからtail部分を検索(一部の不一致を許容するように改善の余地あり)
    #if  re.search("ACTTGCCTGTCGCTCTATCTTC",str(sequence)):
    #    tail_seq = re.search("ACTTGCCTGTCGCTCTATCTTC",str(sequence)).span()
    #elif re.search("TTTCTGTTGGTGCTGATATTGC",str(sequence)):
    #    tail_seq = re.search("TTTCTGTTGGTGCTGATATTGC",str(sequence)).span()
    #else:
    #    unknown.append(rec)
    #    continue

    #tail部分の検索を8塩基で
    if  re.search("CTATCTTC",str(sequence)):
        tail_seq = re.search("CTATCTTC",str(sequence)).span()
    elif re.search("GATATTGC",str(sequence)):
        tail_seq = re.search("GATATTGC",str(sequence)).span()
    else:
        unknown.append(rec)
        continue

    #tail部分の次の塩基から8塩基がindex配列
    #tail部分が末端の150塩基以内に検出された場合に採用
    if tail_seq[1] < 151:
        index_start = tail_seq[1]
        index_end = index_start + 8
    else:
        continue

    index = sequence[index_start:index_end]
    index1_R = "ATTACTCG"
    index2_R = "TCCGGAGA"
    index3_R = "CGCTCATT"
    index4_R = "GAGATTCC"
    index5_R = "AACCTACG"
    index6_R = "GAATTCGT"
    index7_R = "CTGATGAG"
    index8_R = "TAATGCGC"
    index9_R = "CGGCTATG"
    index10_R = "TGCTTGCT"
    index11_R = "TCTCGCGC"
    index12_R = "AGCGATAG"
    index13_R = "ACAGCGTA"
    index14_R = "GTACGTAC"
    index15_R = "GCTGACTA"
    index16_R = "CTAGAGAT"

     #indexを12個に振り分ける
    if Levenshtein.distance(str(index), index1_R) < 3:
        rec1_R.append(rec.id)

    elif Levenshtein.distance(str(index), index2_R) < 3:
        rec2_R.append(rec.id)

    elif Levenshtein.distance(str(index), index3_R) < 3:
        rec3_R.append(rec.id)

    elif Levenshtein.distance(str(index), index4_R) < 3:
        rec4_R.append(rec.id)

    elif Levenshtein.distance(str(index), index5_R) < 3:
        rec5_R.append(rec.id)

    elif Levenshtein.distance(str(index), index6_R) < 3:
        rec6_R.append(rec.id)

    elif Levenshtein.distance(str(index), index7_R) < 3:
        rec7_R.append(rec.id)

    elif Levenshtein.distance(str(index), index8_R) < 3:
        rec8_R.append(rec.id)

    elif Levenshtein.distance(str(index), index9_R) < 3:
        rec9_R.append(rec.id)

    elif Levenshtein.distance(str(index), index10_R) < 3:
        rec10_R.append(rec.id)

    elif Levenshtein.distance(str(index), index11_R) < 3:
        rec11_R.append(rec.id)

    elif Levenshtein.distance(str(index), index12_R) < 3:
        rec12_R.append(rec.id)

    elif Levenshtein.distance(str(index), index13_R) < 3:
        rec13_R.append(rec.id)

    elif Levenshtein.distance(str(index), index14_R) < 3:
        rec14_R.append(rec.id)

    elif Levenshtein.distance(str(index), index15_R) < 3:
        rec15_R.append(rec.id)

    elif Levenshtein.distance(str(index), index16_R) < 3:
        rec16_R.append(rec.id)

    else:
        unknown.append(rec.id)

#ForwardとReverseに共通する要素を格納
rec1_id = list(set(rec1_F) & set(rec1_R))
rec2_id = list(set(rec2_F) & set(rec2_R))
rec3_id = list(set(rec3_F) & set(rec3_R))
rec4_id = list(set(rec4_F) & set(rec4_R))
rec5_id = list(set(rec5_F) & set(rec5_R))
rec6_id = list(set(rec6_F) & set(rec6_R))
rec7_id = list(set(rec7_F) & set(rec7_R))
rec8_id = list(set(rec8_F) & set(rec8_R))
rec9_id = list(set(rec9_F) & set(rec9_R))
rec10_id = list(set(rec10_F) & set(rec10_R))
rec11_id = list(set(rec11_F) & set(rec11_R))
rec12_id = list(set(rec12_F) & set(rec12_R))
rec13_id = list(set(rec13_F) & set(rec13_R))
rec14_id = list(set(rec14_F) & set(rec14_R))
rec15_id = list(set(rec15_F) & set(rec15_R))
rec16_id = list(set(rec16_F) & set(rec16_R))

#idに沿ってリードを配分
recs = SeqIO.parse(gzip.open(in_fq, 'rt'), 'fastq')

rec1 = []
rec2 = []
rec3 = []
rec4 = []
rec5 = []
rec6 = []
rec7 = []
rec8 = []
rec9 = []
rec10 = []
rec11 = []
rec12 = []
rec13 = []
rec14 = []
rec15 = []
rec16 = []

for rec in recs:
    if rec.id in rec1_id:
        rec1.append(rec)

    if rec.id in rec2_id:
        rec2.append(rec)
    if rec.id in rec3_id:
        rec3.append(rec)
    if rec.id in rec4_id:
        rec4.append(rec)
    if rec.id in rec5_id:
        rec5.append(rec)
    if rec.id in rec6_id:
        rec6.append(rec)
    if rec.id in rec7_id:
        rec7.append(rec)
    if rec.id in rec8_id:
        rec8.append(rec)
    if rec.id in rec9_id:
        rec9.append(rec)
    if rec.id in rec10_id:
        rec10.append(rec)
    if rec.id in rec11_id:
        rec11.append(rec)
    if rec.id in rec12_id:
        rec12.append(rec)
    if rec.id in rec13_id:
        rec13.append(rec)
    if rec.id in rec14_id:
        rec14.append(rec)
    if rec.id in rec15_id:
        rec15.append(rec)
    if rec.id in rec16_id:
        rec16.append(rec)


#アダプター配列のトリミング
x = 0
for rec in [rec1, rec2, rec3, rec4, rec5, rec6, rec7, rec8, rec9, rec10, rec11, rec12, rec13, rec14, rec15, rec16]:
    trim_result = []
    x = x + 1
    for read in rec:
        sequence_F = read.seq
        sequence_R = sequence_F.reverse_complement()


        if  re.search("CTATCTTC",str(sequence_F)):
            tail_F = re.search("CTATCTTC",str(sequence_F)).span()
        elif re.search("GATATTGC",str(sequence_F)):
            tail_F = re.search("GATATTGC",str(sequence_F)).span()
        else:
            not_found.append(rec)
            continue

        if  re.search("CTATCTTC",str(sequence_R)):
            tail_R = re.search("CTATCTTC",str(sequence_R)).span()
        elif re.search("GATATTGC",str(sequence_R)):
            tail_R = re.search("GATATTGC",str(sequence_R)).span()
        else:
            unknown.append(rec)
            continue

        #tail部分の次の塩基から8塩基がindex配列
        #tail部分が末端の150塩基以内に検出された場合に採用
        if tail_F[1] < 151:
            index_start_F = tail_F[1]
            index_end_F = index_start_F + 8
        else:
            index_end_F = 0


        #tail部分の次の塩基から8塩基がindex配列
        #tail部分が末端の150塩基以内に検出された場合に採用
        if tail_R[1] < 151:
            index_start_R = tail_F[1]
            index_end_R = index_start_R + 8
        else:
            index_end_R = 0

        index_start_R = len(sequence_F) - index_end_R

        read_trimmed = read[index_end_F:index_start_R]
        trim_result.append(read_trimmed)

    SeqIO.write(trim_result, file_name + "_index" + str(x) + ".fastq","fastq")
