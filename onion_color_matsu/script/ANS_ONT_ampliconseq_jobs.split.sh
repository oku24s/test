#!/bin/sh

#$ -j y


indir=$1
outdir=$2
split_num=$3
sample_list=$4
debarcoding=$5
ANS_h1_threshold=$6
ANS_L_threshold=$7
medaka_model=$8
deepvariant_model=$9
program_path=${10}

REF_list=/usr/local/lib/ANS_ref/ANS_ref.list
min=1600
max=2500

#debarcodingとadapter trimming
#ANSの解析に向けて、リード長のフィルタリングを1.6kb - 2.5kbとした
singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_0_demultiplex_trimming.sh -i ${indir} -o ${outdir} -m ${min} -M ${max} -s ${debarcoding}

#分割に使用する中間ファイル"split_num.txt"を作成
sample_count=`wc -l ${outdir}/analysis/fastq_fil.txt | cut -f 1 -d " "`
repeat=$((${sample_count} / ${split_num}))
mod=$((${sample_count} % ${split_num}))

for x in `seq 1 ${repeat}`
do
        seq 1 ${split_num} >> ${outdir}/analysis/split_num.txt
done

seq 1 ${mod} >> ${outdir}/analysis/split_num.txt

#split_num.txtをfastq_fil.txtにpaste
paste ${outdir}/analysis/fastq_fil.txt ${outdir}/analysis/split_num.txt > ${outdir}/analysis/sample_split.txt

#sample_split.txtを参照してfastqをディレクトリに分割
for x in `seq 1 ${split_num}`
do
        #行先のディレクトリを作成
        mkdir ${outdir}/${x}
        mkdir ${outdir}/${x}/analysis

        cat ${outdir}/analysis/sample_split.txt | awk '{if($2 == '${x}') {print $1}}' > ${outdir}/analysis/sample_split_list.txt
        for sample in `cat ${outdir}/analysis/sample_split_list.txt`
        do
                sample_name=`basename ${sample}`
                mv ${sample} ${outdir}/${x}/analysis
                echo "${outdir}/${x}/analysis/${sample_name}" >> ${outdir}/${x}/analysis/fastq_fil.txt
        done
done


#リファレンスとの適合性チェック
for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_1_reference_check.ANS.sh -o ${outdir}/${x} -r ${REF_list} -f ${ANS_h1_threshold} -t ${ANS_L_threshold} &
done
wait

#1回目の変異検出とハプロタイプ分割
for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_2_haplotype_split.ANS.sh -o ${outdir}/${x} -r ${REF_list} &
done
wait

#分割したリードを再度マッピング
for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_3_minimap2.sh -o ${outdir}/${x} -r ${REF_list} &
done
wait

#変異検出
for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_4_bcftools_sniffles.sh -o ${outdir}/${x} -r ${REF_list} &
done
wait

for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/medaka.sif ONT_processing_4_medaka.sh -o ${outdir}/${x} -r ${REF_list} -m ${medaka_model} &
done
wait

for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/pepper_margin_deepvariant.sif ONT_processing_4_deepvariant.sh -o ${outdir}/${x} -r ${REF_list} -m ${deepvariant_model} &
done
wait
#
##結果の統合
for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_5_merge.sh -o ${outdir}/${x} -r ${REF_list} &
done
wait

#snpEff
for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/snpeff_onion.sif ANS_processing_6_annotation.sh -o ${outdir}/${x} -r ${REF_list} &
done
wait

#分割した各ディレクトリから、結果を移動
singularity exec ${program_path}/container/sif/onion_color.sif ANS_ONT_processing_6_merge_files.sh -o ${outdir} -r ${REF_list} -n ${split_num}

#判定スクリプト
singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_7_output_ANS.py ${outdir}
singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_8_output_ANS_fin.sh -o ${outdir} -r ${REF_list}
singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_9_output_ANS_fin.py ${outdir} ${sample_list}
