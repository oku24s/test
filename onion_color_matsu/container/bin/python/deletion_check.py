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
    target_indel = [180, 466, 488, 1205, 1992]
    target_list = target_indel
elif target_name == "ANS":
    target_indel = [1537, 1635]
    target_list = target_indel
else:
    sys.exit("target gene is required !")

#vcffileでターゲットの変異があった場合にリストに追加
indel_list = []
for record in vcffile:
    sample = record.samples
    if record.POS in target_indel:
        indel_list.append(record.POS)

#bamファイルのインポート
bamfile = pysam.AlignmentFile(bam_name, "rb")


#indel_listが両方とも空の場合、以下の動作は実施せず、プログラムを強制終了
if len(indel_list) == 0:
    #print("homo")
    sys.stdout.write("FALSE")
    sys.exit("indel not found")

#indel_listが2個より多い場合（ターゲットが複数カ所発生した場合）、フェージングさせる？今後アップデートする。：削除

#以下、snp_listとindel_listが合計して１箇所であるとして解析
#indelの解析(DFR)
if target_name == "DFR":
    indel180 = 0
    indel466 = 0
    indel488 = 0
    indel1205 = 0
    indel1992 = 0

    for pos in indel_list:

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
                #cigar"4M3I100D":(0,4), (1,3),(2,100)
                #cigarを先頭からbase_posを上回るまで数え上げ
                #上回るcigarがindelまたはdeletionの場合、indel_lengthとして記録
                x = x + cigar[1]
                if x > base_pos:
                    if cigar[0] == 1:#変数xからcigarに変更
                        #insertion
                        ins_length = cigar[1]
                    elif cigar[0] == 2:
                        #deletion
                        del_length = cigar[1]
            #base_posに応じて、検出されるべきindel_lengthを設定
            if indel_pos == 180 and (del_length > 130 and del_length < 160):
                indel180 = indel180 + 1
        
            elif indel_pos == 466 and (del_length > 500 and del_length < 530):
                indel466 = indel466 + 1

            elif indel_pos == 488 and (del_length > 120 and del_length < 140):
                indel488 = indel488 + 1

            elif indel_pos == 1205 and (ins_length > 2 and ins_length < 5):
                indel1205 = indel1205 + 1

            elif indel_pos == 1992 and del_length < 10:
                indel1992 = indel1992 + 1

    #bamファイル中の総リード数をカウントする
    sum_read = bamfile.mapped
    #検出閾値リード数の設定
    if indel180/sum_read > 0.5 or indel466/sum_read > 0.5 or indel488/sum_read > 0.5 or indel1205/sum_read > 0.5 or indel1992/sum_read > 0.5:
        sys.stdout.write("TRUE")
        sys.exit()


#indelの解析(ANS)
if target_name == "ANS":
    indel1537 = 0
    indel1635 = 0
    
    for pos in indel_list:

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
                #cigar"4M3I100D":(0,4), (1,3),(2,100)
                #cigarを先頭からbase_posを上回るまで数え上げ
                #上回るcigarがindelまたはdeletionの場合、indel_lengthとして記録
                x = x + cigar[1]
                if x > base_pos:
                    if cigar[0] == 1:#変数xからcigarに変更
                        #insertion
                        ins_length = cigar[1]
                    elif cigar[0] == 2:
                        #deletion
                        del_length = cigar[1]
            #base_posに応じて、検出されるべきindel_lengthを設定
            if indel_pos == 1537 and (del_length > 2 and del_length < 6):
                indel1537 = indel1537 + 1
        
            elif indel_pos == 1635 and (del_length > 2 and del_length < 6):
                indel1635 = indel1635 + 1


    #bamファイル中の総リード数をカウントする
    sum_read = bamfile.mapped
    #検出閾値リード数の設定
    if indel1537/sum_read > 0.5 or indel1635/sum_read > 0.5:
        sys.stdout.write("TRUE")
        sys.exit()

bamfile.close()
