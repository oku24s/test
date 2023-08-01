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

t1_i1_group = []
t1_i2_group = []
t1_i3_group = []
t1_i4_group = []
t1_i5_group = []
t1_i6_group = []
t1_i7_group = []
t1_i8_group = []
t1_i9_group = []
t1_i10_group = []
t1_i11_group = []
t1_i12_group = []
t1_i13_group = []
t1_i14_group = []
t1_i15_group = []
t1_i16_group = []
t1_i17_group = []
t1_i18_group = []
t1_i19_group = []
t1_i20_group = []
t1_i21_group = []
t1_i22_group = []
t1_i23_group = []
t1_i24_group = []
t2_i1_group = []
t2_i2_group = []
t2_i3_group = []
t2_i4_group = []
t2_i5_group = []
t2_i6_group = []
t2_i7_group = []
t2_i8_group = []
t2_i9_group = []
t2_i10_group = []
t2_i11_group = []
t2_i12_group = []
t2_i13_group = []
t2_i14_group = []
t2_i15_group = []
t2_i16_group = []
t2_i17_group = []
t2_i18_group = []
t2_i19_group = []
t2_i20_group = []
t2_i21_group = []
t2_i22_group = []
t2_i23_group = []
t2_i24_group = []
not_found = []


tail1 = "ACTTGCCTGTCGCTCTATCTTC"
tail2 = "TTTCTGTTGGTGCTGATATTGC"
index1 = "ATTACTCG"
index2 = "TCCGGAGA"
index3 = "CGCTCATT"
index4 = "GAGATTCC"
index5 = "AACCTACG"
index6 = "GAATTCGT"
index7 = "CTGATGAG"
index8 = "TAATGCGC"
index9 = "CGGCTATG"
index10 = "TGCTTGCT"
index11 = "TCTCGCGC"
index12 = "AGCGATAG"
index13 = "ACAGCGTA"
index14 = "GTACGTAC"
index15 = "GCTGACTA"
index16 = "CTAGAGAT"
index17 = "TATAGCCT"
index18 = "ATAGAGGC"
index19 = "CCTATCCT"
index20 = "GGCTCTGA"
index21 = "ACGCGAAC"
index22 = "TAGTCTTA"
index23 = "CAGGACGT"
index24 = "GTTCTGTC"

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
rec17_F = []
rec18_F = []
rec19_F = []
rec20_F = []
rec21_F = []
rec22_F = []
rec23_F = []
rec24_F = []
unknown = []

tail_list = [tail1, tail2]
index_list = [index1, index2, index3, index4, index5, index6, index7, index8, index9, index10, index11, index12,
        index13, index14, index15, index16, index17, index18, index19, index20, index21, index22, index23, index24]
index_name = ['index1', 'index2', 'index3', 'index4', 'index5', 'index6', 'index7', 'index8', 'index9', 'index10', 'index11', 'index12',
        'index13', 'index14', 'index15', 'index16', 'index17', 'index18', 'index19', 'index20', 'index21', 'index22', 'index23', 'index24']

tail_dic = {tail1:'tail1', tail2:'tail2'}
index_dic = dict(zip(index_list, index_name))

recs = SeqIO.parse(gzip.open(in_fq, 'rt'), 'fastq')

x=0
for rec in recs:
    #indexの初期化
    index_result = ""
    tail_result = ""
    sequence = rec.seq
    #sequenceからtail部分を検索
    #Levenshtein距離が8以下の場合、検出とする。最小のLevenshtein距離になるフレームをtail配列として検出。
    score = 9
    for start in range(0,100):
        end = start + 30
        candidate = str(sequence[start:end])
        for tail in tail_list:
            for index in index_list:
                subject_seq = tail + index
                tail_index_score = Levenshtein.distance(candidate, subject_seq)
                if tail_index_score < score:
                    score = tail_index_score
                    tail_result = tail_dic[tail]
                    index_result = index_dic[index]
    #print(index_result, tail_result, score)
    if index_result == "index1":
        rec1_F.append(rec.id)
    if index_result == "index2":
        rec2_F.append(rec.id)
    if index_result == "index3":
        rec3_F.append(rec.id)
    if index_result == "index4":
        rec4_F.append(rec.id)
    if index_result == "index5":
        rec5_F.append(rec.id)
    if index_result == "index6":
        rec6_F.append(rec.id)
    if index_result == "index7":
        rec7_F.append(rec.id)
    if index_result == "index8":
        rec8_F.append(rec.id)
    if index_result == "index9":
        rec9_F.append(rec.id)
    if index_result == "index10":
        rec10_F.append(rec.id)
    if index_result == "index11":
        rec11_F.append(rec.id)
    if index_result == "index12":
        rec12_F.append(rec.id)
    if index_result == "index13":
        rec13_F.append(rec.id)
    if index_result == "index14":
        rec14_F.append(rec.id)
    if index_result == "index15":
        rec15_F.append(rec.id)
    if index_result == "index16":
        rec16_F.append(rec.id)
    if index_result == "index17":
        rec17_F.append(rec.id)
    if index_result == "index18":
        rec18_F.append(rec.id)
    if index_result == "index19":
        rec19_F.append(rec.id)
    if index_result == "index20":
        rec20_F.append(rec.id)
    if index_result == "index21":
        rec21_F.append(rec.id)
    if index_result == "index22":
        rec22_F.append(rec.id)
    if index_result == "index23":
        rec23_F.append(rec.id)
    if index_result == "index24":
        rec24_F.append(rec.id)


#反対向きから調べる
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
rec17_R = []
rec18_R = []
rec19_R = []
rec20_R = []
rec21_R = []
rec22_R = []
rec23_R = []
rec24_R = []
recs = SeqIO.parse(gzip.open(in_fq, 'rt'), 'fastq')

for rec in recs:
    sequence = rec.seq
    sequence = sequence.reverse_complement()
    #indexの初期化
    index_result = ""
    tail_result = ""
    #sequenceからtail部分を検索
    #Levenshtein距離が5以下の場合、検出とする。最小のLevenshtein距離になるフレームをtail配列として検出。
    score = 9
    for start in range(0,150):
        end = start + 30
        candidate = str(sequence[start:end])
        for tail in tail_list:
            for index in index_list:
                subject_seq = tail + index
                tail_index_score = Levenshtein.distance(candidate, subject_seq)
                if tail_index_score < score:
                    score = tail_index_score
                    tail_result = tail_dic[tail]
                    index_result = index_dic[index]
    #print(index_result, tail_result, score)
    if index_result == "index1":
        rec1_R.append(rec.id)
    if index_result == "index2":
        rec2_R.append(rec.id)
    if index_result == "index3":
        rec3_R.append(rec.id)
    if index_result == "index4":
        rec4_R.append(rec.id)
    if index_result == "index5":
        rec5_R.append(rec.id)
    if index_result == "index6":
        rec6_R.append(rec.id)
    if index_result == "index7":
        rec7_R.append(rec.id)
    if index_result == "index8":
        rec8_R.append(rec.id)
    if index_result == "index9":
        rec9_R.append(rec.id)
    if index_result == "index10":
        rec10_R.append(rec.id)
    if index_result == "index11":
        rec11_R.append(rec.id)
    if index_result == "index12":
        rec12_R.append(rec.id)
    if index_result == "index13":
        rec13_R.append(rec.id)
    if index_result == "index14":
        rec14_R.append(rec.id)
    if index_result == "index15":
        rec15_R.append(rec.id)
    if index_result == "index16":
        rec16_R.append(rec.id)
    if index_result == "index17":
        rec17_R.append(rec.id)
    if index_result == "index18":
        rec18_R.append(rec.id)
    if index_result == "index19":
        rec19_R.append(rec.id)
    if index_result == "index20":
        rec20_R.append(rec.id)
    if index_result == "index21":
        rec21_R.append(rec.id)
    if index_result == "index22":
        rec22_R.append(rec.id)
    if index_result == "index23":
        rec23_R.append(rec.id)
    if index_result == "index24":
        rec24_R.append(rec.id)

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
rec17_id = list(set(rec17_F) & set(rec17_R))
rec18_id = list(set(rec18_F) & set(rec18_R))
rec19_id = list(set(rec19_F) & set(rec19_R))
rec20_id = list(set(rec20_F) & set(rec20_R))
rec21_id = list(set(rec21_F) & set(rec21_R))
rec22_id = list(set(rec22_F) & set(rec22_R))
rec23_id = list(set(rec23_F) & set(rec23_R))
rec24_id = list(set(rec24_F) & set(rec24_R))

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
rec17 = []
rec18 = []
rec19 = []
rec20 = []
rec21 = []
rec22 = []
rec23 = []
rec24 = []

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
    if rec.id in rec17_id:
        rec17.append(rec)
    if rec.id in rec18_id:
        rec18.append(rec)
    if rec.id in rec19_id:
        rec19.append(rec)
    if rec.id in rec20_id:
        rec20.append(rec)
    if rec.id in rec21_id:
        rec21.append(rec)
    if rec.id in rec22_id:
        rec22.append(rec)
    if rec.id in rec23_id:
        rec23.append(rec)
    if rec.id in rec24_id:
        rec24.append(rec)


#アダプター配列のトリミング
x = 0
for rec in [rec1, rec2, rec3, rec4, rec5, rec6, rec7, rec8, rec9, rec10, rec11, rec12, rec13, rec14, rec15, rec16, rec17, rec18, rec19, rec20, rec21, rec22, rec23, rec24]:
    trim_result = []
    x = x + 1
    for read in rec:
        sequence_F = read.seq
        sequence_R = sequence_F.reverse_complement()
        #sequenceからtail部分を検索
        #Levenshtein距離が5以下の場合、検出とする。最小のLevenshtein距離になるフレームをtail配列として検出。
        score = 5
        pos_F = -1
        for start in range(0,150):
            end = start + 22
            candidate = str(sequence_F[start:end])
            if Levenshtein.distance(candidate, tail1) < score:
                score = Levenshtein.distance(candidate, tail1)
                pos_F = end

            elif Levenshtein.distance(candidate, tail2) < score:
                score = Levenshtein.distance(candidate, tail2)
                pos_F = end
        index_end_F = pos_F + 8
        #sequenceからtail部分を検索
        #Levenshtein距離が5以下の場合、検出とする。最小のLevenshtein距離になるフレームをtail配列として検出。
        score = 5
        pos_R = -1
        for start in range(0,150):
            end = start + 22
            candidate = str(sequence[start:end])
            if Levenshtein.distance(candidate, tail1) < score:
                score = Levenshtein.distance(candidate, tail1)
                pos_R = end

            elif Levenshtein.distance(candidate, tail2) < score:
                score = Levenshtein.distance(candidate, tail2)
                pos_R = end
        index_end_R = pos_R + 8


        index_start_R = len(sequence_F) - index_end_R

        read_trimmed = read[index_end_F:index_start_R]
        trim_result.append(read_trimmed)

    SeqIO.write(trim_result, file_name + "_index" + str(x) + ".fastq","fastq")
