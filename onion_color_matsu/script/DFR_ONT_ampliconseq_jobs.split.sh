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
program_path=${10}

REF_list=/usr/local/lib/DFR_ref/DFR_ref.list
min=1300
max=2800

#debarcodingとadapter trimming
#DFRの解析に向けて、リード長のフィルタリングを1.3kb - 2.8kbとした
singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_0_demultiplex_trimming.sh -i ${indir} -o ${outdir} -m ${min} -M ${max} -s ${debarcoding}


#分割に使用する中間ファイル"split_num.txt"を作成
#"split_num.txt"にはサンプルごとに並列処理に向けてサンプルを割り振る行先の番号が記載される
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

#sample_split.txtを参照してfastqを各ディレクトリに割り振る
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


#以下は並列処理のために割り振られたディレクトリごとに実行される
#リファレンスとの適合性チェック
for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_1_reference_check.DFR.sh -o ${outdir}/${x} -r ${REF_list} -f ${DFR_Full_threshold} -t ${DFR_TRN_threshold} &
done
wait

for x in `seq 1 ${split_num}`
do
	#1回目の変異検出とハプロタイプ分割
	singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_2_haplotype_split.DFR.sh -o ${outdir}/${x} -r ${REF_list} &
done
wait

for x in `seq 1 ${split_num}`
do
	#分割したリードを再度マッピング
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
	singularity exec ${program_path}/container/sif/pepper_margin_deepvariant.sif ONT_processing_4_deepvariant.sh -o ${outdir}/${x} -r ${REF_list} -m ${deepvariant_model} &
done
wait

for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/medaka.sif ONT_processing_4_medaka.sh -o ${outdir}/${x} -r ${REF_list} -m ${medaka_model} &
done
wait


#変異検出した結果の統合
for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_5_merge.sh -o ${outdir}/${x} -r ${REF_list} &
done
wait

#snpEffで変異の影響をアノテーション
for x in `seq 1 ${split_num}`
do
	singularity exec ${program_path}/container/sif/snpeff_onion.sif DFR_processing_6_annotation.sh -o ${outdir}/${x} -r ${REF_list} &
done
wait

#分割した各ディレクトリから、結果を移動(singularity内で実行する必要あり）
singularity exec ${program_path}/container/sif/onion_color.sif DFR_ONT_processing_6_merge_files.sh -o ${outdir} -r ${REF_list} -n ${split_num}

#判定スクリプト
singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_7_output_DFR.py ${outdir}
#DPを出力する場合は以下を利用
#singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_7_output_DFR_dp.py ${outdir}
singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_8_output_DFR_fin.sh -o ${outdir} -r ${REF_list}
singularity exec ${program_path}/container/sif/onion_color.sif ONT_processing_9_output_DFR_fin.py ${outdir} ${sample_list}
