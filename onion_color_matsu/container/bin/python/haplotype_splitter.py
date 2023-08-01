#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import vcf
import pysam
from Bio.Seq import Seq
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", required = True)
parser.add_argument("-b", "--bam", required = True)
parser.add_argument("-t", "--target", required = True)

args = parser.parse_args()
vcf_name = args.vcf
bam_name = args.bam
target_name = args.target

vcffile = vcf.Reader(open(vcf_name, "r"))

#ターゲットの変異をどうするか？まずはSNPのみとする
if target_name == "DFR":
    target_snp = [606, 1280]
    target_indel = [180, 466, 488, 1205, 1992]
    target_list = target_snp + target_indel
elif target_name == "ANS":
    target_snp = [1812, 1934]
    target_indel = [1537, 1635]
    target_list = target_snp + target_indel
else:
    sys.exit("target gene is required !")

#vcffileでターゲットの変異がヘテロだった場合にリストに追加
snp_list = []
indel_list = []
for record in vcffile:
    sample = record.samples
    if record.POS in target_snp and sample[0]["GT"] == "0/1":
        snp_list.append(record.REF)
        snp_list.append(str(record.ALT[0]))
        snp_pos = record.POS
    elif record.POS in target_indel and sample[0]["GT"] == "0/1":
        indel_pos = record.POS


#snp_listとindel_listが両方とも空の場合、以下の動作は実施せず、プログラムを強制終了
if len(snp_list) == 0 and len(indel_list) == 0:
    #print("homo")
    sys.stdout.write("homo")
    sys.exit("homo sample")

#snp_listとindel_listが2個より多い場合（ターゲットが複数カ所発生した場合）、フェージングさせる？今後アップデートする。
if len(snp_list) + len(indel_list) > 4:
    #print("multi")
    sys.stdout.write("multi")
    sys.exit("more than 2 mutations were detected!")

bamfile = pysam.AlignmentFile(bam_name, "rb")
hap1_name = bam_name.replace(".match.bam", "_H1.bam")
hap2_name = bam_name.replace(".match.bam", "_H2.bam")
out_bamfile_Hap1 = pysam.AlignmentFile(hap1_name, "wb", template = bamfile)
out_bamfile_Hap2 = pysam.AlignmentFile(hap2_name, "wb", template = bamfile)
Hap1_list = []
Hap2_list = []
Hap1_indel_list = []
Hap2_indel_list = []

#以下、snp_listとindel_listが合計して１箇所であるとして解析
#snpの解析
if len(snp_list) == 2:
    snp_pos = snp_pos - 1 #スライスで文字列を取得するために、posを-1

    for read in bamfile:
        seq = read.get_forward_sequence() #get_forward_sequence() returns the original read sequence.
        if read.is_reverse:
            seq = str(Seq(seq).reverse_complement())
        ref_pos = read.get_reference_positions(full_length = True)
        try:
            base_pos = ref_pos.index(snp_pos) #リファレンス上のsnp_pos番目の塩基にマッピングされた塩基はリードでは何番目の塩基か
        except:
            continue
        snp = seq[base_pos]
        if snp == snp_list[0]:
            #print(read.query_name,snp, base_pos,"Hap1")
            Hap1_list.append(read.query_name)
        elif snp == snp_list[1]:
            #print(read.query_name,snp, base_pos,"Hap2")
            Hap2_list.append(read.query_name)

#indelの解析(DFR)
if target_name == "DFR" and len(indel_list) == 2:
    #indelが検出された位置のcigarを確認.cigarの合計がindel_posを上回るまでカウント
    indel_pos = indel_pos[0]
    for read in bamfile:
        cigar_tuples = read.cigartuples
        x = 0
        ref_pos = read.get_reference_positions(full_length = True)
        base_pos = ref_pos.index(indel_pos)#リファレンス上のindel_pos番目はリードでは何番目か
        ins_length = 0
        del_length = 0
        for cigar in cigar_tuples:
            #cigarを先頭からbase_posを上回るまで数え上げ
            #上回るcigarがindelまたはdeletionの場合、indel_lengthとして記録
            x = x + cigar[1]
            if x > base_pos:
                if x[0] == 1:
                    #insertion
                    ins_length = x[1]
                elif x[0] == 2:
                    #deletion
                    del_length = x[1]
        #base_posに応じて、検出されるべきindel_lengthを設定
        if indel_pos == 180: #DFRの180番目には140塩基のdeletion
            if del_length > 130 and del_length < 160:
                Hap1_indel_list.append(read.query_name)
            else:
                Hap2_indel_list.append(read.query_name)
        
        elif indel_pos == 466: #DFRの466番目には516塩基のdeletion
            if del_length > 500 and del_length < 530:
                Hap1_indel_list.append(read.query_name)
            else:
                Hap2_indel_list.append(read.query_name)

        elif indel_pos == 488: #DFRの488番目には130塩基のdeletion
            if del_length > 120 and del_length < 140:
                Hap1_indel_list.append(read.query_name)
            else:
                Hap2_indel_list.append(read.query_name)

        elif indel_pos == 1205: #DFRの1205番目には3塩基のinsertion
            if ins_length > 2 and ins_length < 5:
                Hap1_indel_list.append(read.query_name)
            else:
                Hap2_indel_list.append(read.query_name)

        elif indel_pos == 1992: #DFRの1992番目には1塩基のdeletion
            if del_length < 10:
                Hap1_indel_list.append(read.query_name)
            else:
                Hap2_indel_list.append(read.query_name)

#indelの解析(ANS)
if target_name == "ANS" and len(indel_list) == 2:
    #indelが検出された位置のcigarを確認.cigarの合計がindel_posを上回るまでカウント
    indel_pos = indel_pos[0]
    for read in bamfile:
        cigar_tuples = read.cigartuples
        x = 0
        ref_pos = read.get_reference_positions(full_length = True)
        base_pos = ref_pos.index(indel_pos)#リファレンス上のindel_pos番目はリードでは何番目か
        ins_length = 0
        del_length = 0
        for cigar in cigar_tuples:
            #cigarを先頭からbase_posを上回るまで数え上げ
            #上回るcigarがindelまたはdeletionの場合、indel_lengthとして記録
            x = x + cigar[1]
            if x > base_pos:
                if x[0] == 1:
                    #insertion
                    ins_length = x[1]
                elif x[0] == 2:
                    #deletion
                    del_length = x[1]
        #base_posに応じて、検出されるべきindel_lengthを設定
        if indel_pos == 180: #ANSの1537番目には4塩基のinsertion
            if del_length > 2 and del_length < 6:
                Hap1_indel_list.append(read.query_name)
            else:
                Hap2_indel_list.append(read.query_name)
        
        elif indel_pos == 466: #ANSの1635番目には4塩基のinsertion
            if del_length > 2 and del_length < 6:
                Hap1_indel_list.append(read.query_name)
            else:
                Hap2_indel_list.append(read.query_name)

#Hap1,Hap2,Hap1_indel,Hap2_indelのphasing
if len(snp_list) == 2 and len(indel_list) == 2:
    H1_H1 = set(Hap1_list) & set(Hap1_indel_list)
    H1_H2 = set(Hap1_list) & set(Hap2_indel_list)
    H2_H1 = set(Hap2_list) & set(Hap1_indel_list)
    H2_H2 = set(Hap2_list) & set(Hap2_indel_list)
    H1H1_H2H2 = len(H1_H1) + len(H2_H2)
    H1H2_H2H1 = len(H1_H2) + len(H2_H1)
    if H1H1_H2H2 > H1H2_H2H1:
        Hap1_list = list(H1_H1)
        Hap2_list = list(H2_H2)
    else:
        Hap1_list = list(H1_H2)
        Hap2_list = list(H2_H1)
#print(Hap1_list, Hap2_list)

#indelしか検出されなかった場合
if len(snp_list) == 0:
    Hap1_list = Hap1_indel_list
    Hap2_list = Hap2_indel_list


#hap1とhap2のリード数と比率に関して条件を設定
if (len(Hap1_list) > 2 and len(Hap2_list) > 2):
    bamfile = pysam.AlignmentFile(bam_name, "rb")
    for read in bamfile:
        if read.query_name in Hap1_list:
            out_bamfile_Hap1.write(read)
        elif read.query_name in Hap2_list:
            out_bamfile_Hap2.write(read)

else:
    #print("homo")
    sys.stdout.write("homo")



bamfile.close()
out_bamfile_Hap1.close()
out_bamfile_Hap2.close()
