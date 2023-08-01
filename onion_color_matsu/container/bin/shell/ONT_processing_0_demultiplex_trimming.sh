#!/bin/bash

usage (){
	cat <<EOM
Usage: $(basename "$0") [OPTION]...
	-h		Display this help
	-i [str] 	full path of input directory (~/fastq_pass)
	-o [str] 	full path of output directory (~/analysis_name)
	-m [num] 	min length of reads 
	-M [num] 	max length of reads
	-s [str]	strategy for analyzing (select from 'old', 'old_complete', or 'v4')
EOM

	exit 2
}

while getopts ":i:o:m:M:s:h" optKey
do
	case ${optKey} in
		i)
			indir=${OPTARG}
			;;
		o)
			outdir=${OPTARG}
			;;
		m)
			len_min=${OPTARG}
			;;
		M)
			len_max=${OPTARG}
			;;
		s)
			strategy=${OPTARG}
			;;
		'-h'|'--help'|*)
			usage
			;;
	esac
done

if [ ${indir: -1} = '/' ]
then
	indir=${indir/%?/}
fi


if [ ${outdir: -1} = '/' ]
then
	outdir=${outdir/%?/}
fi


char1=`echo ${indir} | cut -c 1`
if [ -z "${indir}" ] || [ "${char1}" = '-' ]
then
	echo "input directory is necessary!!"
	exit
else
	echo "-i (input directory) = ${indir}"
fi

char2=`echo ${outdir} | cut -c 1`
if [ -z "${outdir}" ] || [ "${char2}" = '-' ]
then
	echo "output directory is necessary!!"
	exit
else
	echo "-o (output directory) = ${outdir}"
fi

char3=`echo ${len_min} | cut -c 1`
if [ -z "${len_min}" ] || [ "${char3}" = '-' ]
then
	echo "Min length of reads is missing. -m = 1"
	len_min=1
else
	echo "-m (min length of reads) = ${len_min}"
fi

char4=`echo ${len_max} | cut -c 1`
if [ -z "${len_max}" ] || [ "${char4}" = '-' ]
then
	echo "Max length of reads is missing. -M = 100000000"
	len_max=100000000
else
	echo "-M (max length of reads) = ${len_max}"
fi

char5=`echo ${strategy} | cut -c 1`
if [ -z "${strategy}" ] || [ "${char5}" = '-' ]
then
	echo "debarcoding strategy is empty. This script is running by default setting (old)!"
	strategy=old
else
	echo "-s (debarcoding strategy) = ${strategy}"
fi
#
#
mkdir ${outdir}/analysis

##バーコードの振り分け
for i in `seq -w 1 96`
do
	cat ${indir}/barcode${i}/*fastq.gz > ${outdir}/analysis/barcode${i}.me.fastq.gz 2>/dev/null
	cat ${indir}/barcode${i}/*fastq > ${outdir}/analysis/barcode${i}.merge.fastq 2>/dev/null
	bgzip ${outdir}/analysis/barcode${i}.merge.fastq 2>/dev/null
	cat ${outdir}/analysis/barcode${i}.merge.fastq.gz ${outdir}/analysis/barcode${i}.me.fastq.gz > ${outdir}/analysis/barcode${i}.fastq.gz 2>/dev/null
	rm ${outdir}/analysis/barcode${i}.me.fastq.gz 2>/dev/null
	rm ${outdir}/analysis/barcode${i}.merge.fastq.gz 2>/dev/null

	if [ ${strategy} = 'old' ]
	then
		ONT_demultiplex_trimming_old.py -i ${outdir}/analysis/barcode${i}.fastq.gz
	elif [ ${strategy} = 'v4' ]
	then
		ONT_demultiplex_trimming.v4.py -i ${outdir}/analysis/barcode${i}.fastq.gz
	else
		ONT_demultiplex_trimming_old_complete.py -i ${outdir}/analysis/barcode${i}.fastq.gz
	fi

	rm -f ${outdir}/analysis/barcode${i}.fastq.gz
done

ls ${outdir}/analysis/barcode*_index*.fastq > ${outdir}/analysis/fastq.list

#fastqからreadの長さが${len_min}以上${len_max}以下のみにフィルタリング後、readが5以上検出されたもののみ解析に使用
for x in `cat ${outdir}/analysis/fastq.list`
do
	w=`echo ${x} | sed s/.fastq//g`
	seqkit seq --remove-gaps -m ${len_min} -M ${len_max} ${w}.fastq > ${w}_fil.fastq
	if [ `wc -l ${w}_fil.fastq | awk '{print $1}'` -gt 19 ]
	then 
		echo "${w}_fil.fastq" >> ${outdir}/analysis/fastq_fil.txt
	else
		rm ${w}_fil.fastq -f
	fi

	rm ${w}.fastq -f

done

