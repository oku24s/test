#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import vcf
import csv
import sys

#DFR_Allele_dict.pyからアリルの辞書をインポート
sys.path.append("/usr/local/lib/Allele_dict")
import DFR_Allele_dict_active
import DFR_Allele_dict_inactive
import DFR_accepted_pos
dic_active = DFR_Allele_dict_active.dic
dic_inactive = DFR_Allele_dict_inactive.dic
DFR_accepted_pos = DFR_accepted_pos.LIST

args = sys.argv
annotation_list = args[1] + "/DFR_OGBv3.0_Scaffold003_T1280A.ann.list"
print(annotation_list)

with open(annotation_list) as f:
    
    active_sample = list()
    inactive_sample = list()
    pos_all = list()

    for sample in f:
        input = sample.replace("\n","")
        vcff = vcf.Reader(filename = input)
        input = input.split("/")[-1]
        sample_name2 = input.split(".")[0]
        sample_name = sample_name2.replace("_snp_indel_SV_ann","")
        print(sample_name)

        colour = list()
        pos_list = list()
        
        for record in vcff:
            #色の判定
            ann = record.INFO["ANN"][0].split("|")
            gt = record.genotype(sample_name)["GT"]
            pos = record.POS
            if pos in DFR_accepted_pos:
                #p.Val83Ile,p.Phe271Ile,p.Gln336Gluは色に影響を与えない変異として許容
                #得られたデータからはp.Asn185Tyrも色に影響を与えないと推測された
                if (ann[2] == "MODERATE" or ann[2] == "HIGH") and (gt == "1/1" or gt == "0/1" or gt == "1"):
                    if ann[10] != "p.Val83Ile" and ann[10] != "p.Phe271Ile" and ann[10] != "p.Gln336Glu" and ann[10] != "p.Asn185Tyr" and ann[10] != "p.His234Arg":
                        colour.append([record.POS, ann[1], ann[10]])

                #POS=180にある140塩基のdeletionは非機能型になることが報告されている
                #175-185の位置にdeletionが検出された場合は非機能型とカウントする
                if (175 < pos < 185) and (gt == "1/1" or gt == "0/1" or gt == "1"):
                    if "SVTYPE" in record.INFO:
                        SVtype = record.INFO["SVTYPE"]
                        if SVtype == "DEL":
                            colour.append([record.POS, SVtype, ann[1]])

                #splice_region_variantをもつ場合は非機能型アリルとしてカウントする

                #アリルの判定
                if gt =="1/1" or gt =="0/1" or gt == "1":
                    #pos = record.POS
                    pos_list.append(pos)
        
        pos_all.append(pos_list)
        pos_tup = tuple(pos_list)

        if len(colour) == 0:
            print(input.split(".")[0])
            print("は機能型です")
            print(colour)
            if len(pos_tup) == 0:
                Allele = "DFR_Allele_original"
            else:
                if pos_tup in dic_active:
                    val = dic_active[pos_tup]
                    Allele = val
                else:
                    Allele = "未知のアリルです"
            #active_sample.append([input.split(".")[0], Allele, "active", pos_tup, colour])
            active_sample.append([sample_name, Allele, "active", pos_tup, colour])
        else:
            print(input.split(".")[0])
            print("は非機能型です")
            print(colour)
        
            if len(pos_tup) == 0:
                Allele = "DFR_Allele_original"
            else:
                if pos_tup in dic_inactive:
                    val = dic_inactive[pos_tup]
                    Allele = val
                else:
                    Allele = "未知のアリルです"
            #inactive_sample.append([input.split(".")[0], Allele, "inactive", pos_tup, colour])
            inactive_sample.append([sample_name, Allele, "inactive", pos_tup, colour])


active = args[1] + "/active_sample.tsv"
inactive = args[1] + "/inactive_sample.tsv"

with open(active, mode="w", newline="", encoding="utf-8") as fo:
    tsv_writer = csv.writer(fo, delimiter='\t')
    tsv_writer.writerows(active_sample)


with open(inactive, mode="w", newline="", encoding="utf-8") as fo:
    tsv_writer = csv.writer(fo, delimiter='\t')
    tsv_writer.writerows(inactive_sample)

