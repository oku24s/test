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

REF_list=/usr/local/lib/ANS_ref/ANS_ref.list
min=1600
max=2500

#debarcodingとadapter trimming
#ANSの解析に向けて、リード長のフィルタリングを1.6kb - 2.5kbとした
singularity exec ../container/onion_color.sif ONT_processing_0_demultiplex_trimming.sh -i ${indir} -o ${outdir} -m ${min} -M ${max} -s ${debarcoding}

#リファレンスとの適合性チェック
singularity exec ../container/onion_color.sif ONT_processing_1_reference_check.ANS.sh -o ${outdir}/${x} -r ${REF_list} -f ${ANS_h1_threshold} -t ${ANS_L_threshold}

#1回目の変異検出とハプロタイプ分割
singularity exec ../container/onion_color.sif ONT_processing_2_haplotype_split.DFR.sh -o ${outdir}/${x} -r ${REF_list}

#分割したリードを再度マッピング
singularity exec ../container/onion_color.sif ONT_processing_3_minimap2.sh -o ${outdir}/${x} -r ${REF_list}

#変異検出
singularity exec ../container/onion_color.sif ONT_processing_4_bcftools_sniffles.sh -o ${outdir}/${x} -r ${REF_list}

singularity exec ../container/medaka.sif ONT_processing_4_medaka.sh -o ${outdir}/${x} -r ${REF_list} -m ${medaka_model}

singularity exec ../pepper_margin_deepvariant.sif ONT_processing_4_deepvariant.sh -o ${outdir}/${x} -r ${REF_list} -m ${deepvariant_model}

##結果の統合
singularity exec ../container/onion_color.sif ONT_processing_5_merge.sh -o ${outdir}/${x} -r ${REF_list}

#snpEff
singularity exec ../container/snpeff_onion.sif DFR_processing_6_annotation.sh -o ${outdir}/${x} -r ${REF_list}

#分割した各ディレクトリから、結果を移動
singularity exec ../container/onion_color.sif ANS_ONT_processing_6_ddmerge_files.sh -o ${outdir} -r ${REF_list} -n ${split_num}

#判定スクリプト
singularity exec ../container/onion_color.sif ONT_processing_7_output_ANS.py ${outdir}
singularity exec ../container/onion_color.sif ONT_processing_8_output_ANS_fin.sh -o ${outdir} -r ${REF_list}
singularity exec ../container/onion_color.sif ONT_processing_9_output_ANS_fin.py ${outdir} ${sample_list}
