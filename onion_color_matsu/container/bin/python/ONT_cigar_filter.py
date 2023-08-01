#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#cigarを参照してbamファイルをフィルタリングする. アダプタートリミング済みのfastqから作成したbamファイルのみ処理可能.

import pysam
import argparse

#inputを指定
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
args = parser.parse_args()
print(args, input)

if args.input == "None":
    print("inputのbamファイルを指定してください")
    sys.exit()

in_bam = args.input
in_bamfile = pysam.AlignmentFile(in_bam, "rb")

#outputを作成
out_name = in_bam.replace("sorted.bam","")
out_match = out_name + "match.bam"
out_unmatch = out_name + "unmatch.bam"
out_bamfile_match = pysam.AlignmentFile(out_match, "wb", template = in_bamfile)
out_bamfile_unmatch = pysam.AlignmentFile(out_unmatch, "wb", template = in_bamfile)

#cigarをチェックして、末端にclipが含まれる配列を除去
for read in in_bamfile:
    cigar = read.cigartuples
    if cigar is None:
        continue
    X = len(cigar) - 1 
    cigar_head = cigar[0]
    cigar_tail = cigar[X]
    #4はsoft_clip, 5はhard clipを示す
    #cigar_head[1]はcigarの長さ.cigar_head[1] < 50は50塩基以上のclipではない場合、許容する.
    if (((cigar_head[0] != 4 and cigar_head[0] != 5) or cigar_head[1] < 50) and ((cigar_tail[0] != 4 and cigar_tail[0] != 5) or cigar_tail[1] < 50)):
        out_bamfile_match.write(read)
    else:
        out_bamfile_unmatch.write(read)

in_bamfile.close()
out_bamfile_match.close()
out_bamfile_unmatch.close()
