#!/bin/bash

##singularity shell ONT_ampliconseq.sifで使用

usage (){
	cat <<EOM
Usage: $(basename "$0") [OPTION]...
	-h		Display this help
	-o[str]		full path of output directory (~/analysis_name)
	-r[str]		full path of reference list (~/reference.list)
EOM

	exit 2
}

while getopts ":o:r:h" optKey
do
	case ${optKey} in
		o)
			outdir=${OPTARG}
			;;
		r)
			REF_list=${OPTARG}
			;;
		'-h'|'--help'|*)
			usage
			;;
	esac
done

char1=`echo ${outdir} | cut -c 1`
if [ -z "${outdir}" ] || [ "${char1}" = '-' ]
then
	echo "output directory is necessary!!"
	exit
else
	echo "-o (output directory) = ${outdir}"
fi

char2=`echo ${REF_list} | cut -c 1`
if [ -z "${REF_list}" ] || [ "${char2}" = '-' ]
then
	echo "reference list is necessary!!"
	exit
else
	echo "-r (reference list) = ${REF_list}"
fi


#各リファレンスの適合性を検討
mkdir ${outdir}/ref_check


#結果を出力するディレクトリを作る
for REF in `cat ${REF_list}`
do
	Name=`basename ${REF}`
	Name2=`echo ${Name} | sed s/\.[^\.]*$//g`
	mkdir ${outdir}/${Name2}
done

#unmatch_readを置くディレクトリを作る
mkdir ${outdir}/unmatch_read

#fastqをunmatchディレクトリにコピーして、unmatch.fastqに名称を変更
for x in `cat ${outdir}/analysis/fastq_fil.txt`
do
	#rename=`basename ${x} | sed s/\.[^\.]*$//g`
	rename=`basename ${x} | sed s/\.[^\.]*$//g | sed s/_fil//g`
	cp ${x} ${outdir}/unmatch_read/${rename}.unmatch.fastq
	echo "${outdir}/unmatch_read/${rename}.unmatch.fastq" >> ${outdir}/unmatch_read/unmatch_read.list
done



#最適なリファレンスを判定して、bamファイルを格納する
for x in `cat ${outdir}/unmatch_read/unmatch_read.list`
do
	for REF in `cat ${REF_list}`
	do
		Name=`basename ${REF}`
		Name2=`echo ${Name} | sed s/\.[^\.]*$//g`
		#mapping
		w=`echo ${x} | cut -f 1 -d "."`
		sample_name=`basename ${w}`
		minimap2 -R "@RG\tID:X\tLB:Y\tSM:${w}\tPL:ONT" -ax map-ont --MD ${REF} ${x} |\
		samtools sort -O BAM > ${outdir}/${Name2}/${sample_name}.sorted.bam

		samtools index ${outdir}/${Name2}/${sample_name}.sorted.bam

		#bamファイルのcigarから適合性を判定(.match.bamが出力される)
		ONT_cigar_filter.py -i ${outdir}/${Name2}/${sample_name}.sorted.bam
		#2リード以上残らなかった場合はファイルを削除
		if [ `samtools view -F 0x04 -c ${outdir}/${Name2}/${sample_name}.match.bam` -lt 2 ]; then
			rm ${outdir}/${Name2}/${sample_name}.match.bam
		else
			echo "${outdir}/${Name2}/${sample_name}.match.bam" >> ${outdir}/${Name2}/${Name2}.match_bam.list

		fi

		#処理後のsorted.bamを削除
		rm ${outdir}/${Name2}/${sample_name}.sorted.bam 
		rm ${outdir}/${Name2}/${sample_name}.sorted.bam.bai

		#unmatch.bamからfastqに変換
		samtools fastq ${outdir}/${Name2}/${sample_name}.unmatch.bam > ${outdir}/unmatch_read/${sample_name}.unmatch.fastq
		
		rm ${outdir}/${Name2}/${sample_name}.unmatch.bam
	done
done

#unmatch_read.listを削除
rm ${outdir}/unmatch_read/unmatch_read.list
