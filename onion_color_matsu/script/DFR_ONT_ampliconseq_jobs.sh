#!/bin/sh

#$ -j y


indir=$1
outdir=$2
split_num=$3
sample_list=$4
debarcoding=$5
DFR_Full_threshold=$6
DFR_TRN_threshold=$7
medaka_model=$8
deepvariant_model=$9

REF_list=/usr/local/lib/DFR_ref/DFR_ref.list
min=1300
max=2800

#debarcodingとadapter trimming
#DFRの解析に向けて、リード長のフィルタリングを1.3kb - 2.8kbとした
singularity exec ../container/onion_color.sif ONT_processing_0_demultiplex_trimming.sh -i ${indir} -o ${outdir} -m ${min} -M ${max} -s ${debarcoding}


#リファレンスとの適合性チェック
singularity exec ../container/onion_color.sif ONT_processing_1_reference_check.DFR.sh -o ${outdir}/${x} -r ${REF_list} -f ${DFR_Full_threshold} -t ${DFR_TRN_threshold}


#1回目の変異検出とハプロタイプ分割
singularity exec ../container/onion_color.sif ONT_processing_2_haplotype_split.DFR.sh -o ${outdir}/${x} -r ${REF_list}


#分割したリードを再度マッピング
singularity exec ../container/onion_color.sif ONT_processing_3_minimap2.sh -o ${outdir}/${x} -r ${REF_list}


#変異検出
singularity exec ../container/onion_color.sif ONT_processing_4_bcftools_sniffles.sh -o ${outdir}/${x} -r ${REF_list}

singularity exec ../container/medaka.sif ONT_processing_4_medaka.sh -o ${outdir}/${x} -r ${REF_list} -m ${medaka_model}

singularity exec ../container/pepper_margin_deepvariant.sif ONT_processing_4_deepvariant.sh -o ${outdir}/${x} -r ${REF_list} -m ${deepvariant_model}


#変異検出した結果の統合
singularity exec ../container/onion_color.sif ONT_processing_5_merge.sh -o ${outdir}/${x} -r ${REF_list}


#変異の影響をアノテーション
singularity exec ../container/snpeff_onion.sif DFR_processing_6_annotation.sh -o ${outdir}/${x} -r ${REF_list}


#分割した各ディレクトリから、結果を移動(singularity内で実行する必要あり）
singularity exec ../container/onion_color.sif DFR_ONT_processing_6_merge_files.sh -o ${outdir} -r ${REF_list} -n ${split_num}


#判定スクリプト
singularity exec ../container/onion_color.sif ONT_processing_7_output_DFR.py ${outdir}
singularity exec ../container/onion_color.sif ONT_processing_8_output_DFR_fin.sh -o ${outdir} -r ${REF_list}
singularity exec ../container/onion_color.sif ONT_processing_9_output_DFR_fin.py ${outdir} ${sample_list}
