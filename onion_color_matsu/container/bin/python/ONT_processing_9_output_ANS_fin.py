#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import sys

#result.csvとsample_file.csvを読み込む
args = sys.argv
result_file = args[1] + "/result.csv"
#sample_file.csvを指定
sample_file = args[2]

print(result_file)

result_csv = open(result_file, "r")
sample_csv = open(sample_file, "r", encoding="cp932")

result = csv.reader(result_csv, delimiter=",")
sample = csv.reader(sample_csv, delimiter=",")

fin_list = []
sample_list = []
result_list = []

for x in sample:
    sample_list.append([x[1],x[2]])

for y in result:
    result_list.append(y)

for x in sample_list:
    alleles = [x[0],x[1]]
    active = 0
    for y in result_list:
        if x[0] == y[0]:
            alleles.append(y[1])
            active = active + int(y[3])

    if active > 0:
        alleles.append("active")
    else:
        alleles.append("inactive")

    if len(alleles) == 3:
        del alleles[-1]
        alleles.append("NA")
        alleles.append("NA")
        alleles.append("NA")

    if len(alleles) == 4:
        alleles.insert(2, alleles[2])

    fin_list.append(alleles)




output = args[1] + "/output.tsv"

result_csv.close()
sample_csv.close()

with open(output, mode="w", newline="", encoding="utf8") as fo:
    tsv_writer = csv.writer(fo, delimiter='\t')
    tsv_writer.writerows(fin_list)

with open(output, mode="w", newline="", encoding="utf8") as fo:
    tsv_writer = csv.writer(fo, delimiter='\t')
    tsv_writer.writerows(fin_list)
