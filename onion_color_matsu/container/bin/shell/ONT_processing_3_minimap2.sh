#!/bin/bash
##singularity shell ONT_ampliconseq.sifで使用

usage (){
	cat <<EOM
Usage: $(basename "$0")[OPTION]...
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


#リファレンスごとにH1,H2,homo.fastqを再度mapping

for REF in `cat ${REF_list}`
do
	Name=`basename ${REF}`
	Name2=`echo ${Name} | sed s/\.[^\.]*$//g`
	#fastqのリストを作る.fastqをフィルタリングするなら、ここで行う.
	ls ${outdir}/${Name2}/*_H1.fastq > ${outdir}/${Name2}/fastq_fil_hap.list
	ls ${outdir}/${Name2}/*_H2.fastq >> ${outdir}/${Name2}/fastq_fil_hap.list
	ls ${outdir}/${Name2}/*_homo.fastq >> ${outdir}/${Name2}/fastq_fil_hap.list

	for x in `cat ${outdir}/${Name2}/fastq_fil_hap.list`
	do
		#mapping
		w=`echo ${x} | sed s/.fastq//g`
		sample_name=`basename ${w}`
		minimap2 -R "@RG\tID:X\tLB:Y\tSM:${sample_name}\tPL:ONT" -ax map-ont --MD ${REF} ${x} |\
		samtools sort -O BAM > ${outdir}/${Name2}/${sample_name}.sorted.bam

		samtools index ${outdir}/${Name2}/${sample_name}.sorted.bam

		echo "${outdir}/${Name2}/${sample_name}.sorted.bam" >> ${outdir}/${Name2}/${Name2}_hap_bam.list
	done

done

