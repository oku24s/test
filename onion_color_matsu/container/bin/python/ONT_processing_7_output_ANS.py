#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import vcf
import csv
import sys

#ANS_h1の処理
#ANS_h1_accepted_posから変異として検出する箇所のlistをインポート
#ANS_Allele_dict_h1.pyからアリルの辞書をインポート
sys.path.append("/usr/local/lib/Allele_dict")
import ANS_h1_accepted_pos
import ANS_Allele_dict_h1_active
import ANS_Allele_dict_h1_inactive
dic_active = ANS_Allele_dict_h1_active.dic
dic_inactive = ANS_Allele_dict_h1_inactive.dic
ANS_h1_accepted_pos = ANS_h1_accepted_pos.LIST

args = sys.argv
annotation_list_h1 = args[1] + "/ANS_h1_amplicon.ann.list"
print(annotation_list_h1)

with open(annotation_list_h1) as f:
    
    active_sample = list()
    inactive_sample = list()
    pos_all = list()

    for sample in f:
        input = sample.replace("\n","")
        vcff = vcf.Reader(filename = input)
        input = input.split("/")[-1]
        sample_name_bar = input.split("_")[0]
        sample_name_idx = input.split("_")[1]
        sample_name_type = input.split("_")[2]
        sample_name = sample_name_bar + "_" + sample_name_idx + "_" + sample_name_type

        colour = list()
        pos_list = list()
        
        for record in vcff:
            #色の判定とアリルの判定
            ann = record.INFO["ANN"][0].split("|")
            gt = record.genotype(sample_name)["GT"]
            if gt =="1/1" or gt =="0/1" or gt == "1":
                pos = record.POS
                #if pos in ANS_h1_accepted_pos:
                pos_list.append(pos)
                if ann[2] == "MODERATE" and (gt == "1/1" or gt == "0/1" or gt == "1"):
                    if (ann[10] == "p.Ser188Leu" or ann[10] == "p.Gly229Arg"):
                        colour.append([record.POS, ann[1], ann[10]])
                elif ann[2] == "HIGH" and (gt == "1/1" or gt == "0/1" or gt == "1"):
                    colour.append([record.POS, ann[1], ann[10]])


        
        pos_all.append(pos_list)
        pos_tup = tuple(pos_list)

        if len(colour) == 0:
            print(input.split(".")[0])
            print("は機能型です")
            print(colour)
            if len(pos_tup) == 0:
                Allele = "ANS_h1_original"
            else:
                if pos_tup in dic_active:
                    val = dic_active[pos_tup]
                    Allele = val
                else:
                    Allele = "未知のアリルです"
            active_sample.append([input.split(".")[0], Allele, "active", pos_tup, colour])
        else:
            print(input.split(".")[0])
            print("は非機能型です")
            print(colour)
        
            if len(pos_tup) == 0:
                Allele = "ANS_h1_original"
            else:
                if pos_tup in dic_inactive:
                    val = dic_inactive[pos_tup]
                    Allele = val
                else:
                    Allele = "未知のアリルです"
            inactive_sample.append([input.split(".")[0], Allele, "inactive", pos_tup, colour])

#ANS_Lの処理
import ANS_L_accepted_pos
import ANS_Allele_dict_L_active
import ANS_Allele_dict_L_inactive
dic_active = ANS_Allele_dict_L_active.dic
dic_inactive = ANS_Allele_dict_L_inactive.dic
ANS_L_accepted_pos = ANS_L_accepted_pos.LIST

args = sys.argv
annotation_list_L = args[1] + "/ANS_L_amplicon.ann.list"
print(annotation_list_L)
with open(annotation_list_L) as f:
    
    for sample in f:
        input = sample.replace("\n","")
        vcff = vcf.Reader(filename = input)
        input = input.split("/")[-1]
        sample_name_bar = input.split("_")[0]
        sample_name_idx = input.split("_")[1]
        sample_name_type = input.split("_")[2]
        sample_name = sample_name_bar + "_" + sample_name_idx + "_" + sample_name_type

        colour = list()
        pos_list = list()
        
        for record in vcff:
            #色の判定
            ann = record.INFO["ANN"][0].split("|")
            gt = record.genotype(sample_name)["GT"]
            #アリルの判定
            if gt =="1/1" or gt =="0/1" or gt == "1":
                pos = record.POS
                #if pos in ANS_L_accepted_pos:
                pos_list.append(pos)
                if ann[2] == "MODERATE" and (gt == "1/1" or gt == "0/1" or gt == "1"):
                    if (ann[10] == "p.Ser188Leu" or ann[10] == "p.Gly229Arg"):
                        colour.append([record.POS, ann[1], ann[10]])
                elif ann[2] == "HIGH" and (gt == "1/1" or gt == "0/1" or gt == "1"):
                    colour.append([record.POS, ann[1], ann[10]])


        
        pos_all.append(pos_list)
        pos_tup = tuple(pos_list)

        if len(colour) == 0:
            print(input.split(".")[0])
            print("は機能型です")
            print(colour)
            if len(pos_tup) == 0:
                Allele = "ANS_L_original"
            else:
                if pos_tup in dic_active:
                    val = dic_active[pos_tup]
                    Allele = val
                else:
                    Allele = "未知のアリルです"
            active_sample.append([input.split(".")[0], Allele, "active", pos_tup, colour])
        else:
            print(input.split(".")[0])
            print("は非機能型です")
            print(colour)
        
            if len(pos_tup) == 0:
                Allele = "ANS_L_original"
            else:
                if pos_tup in dic_inactive:
                    val = dic_inactive[pos_tup]
                    Allele = val
                else:
                    Allele = "未知のアリルです"
            inactive_sample.append([input.split(".")[0], Allele, "inactive", pos_tup, colour])



active = args[1] + "/active_sample.tsv"
inactive = args[1] + "/inactive_sample.tsv"

with open(active, mode="w", newline="", encoding="utf-8") as fo:
    tsv_writer = csv.writer(fo, delimiter='\t')
    tsv_writer.writerows(active_sample)


with open(inactive, mode="w", newline="", encoding="utf-8") as fo:
    tsv_writer = csv.writer(fo, delimiter='\t')
    tsv_writer.writerows(inactive_sample)

