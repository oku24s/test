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



for REF in `cat ${REF_list}`
do
	Name=`basename ${REF}`
	Name2=`echo ${Name} | sed s/\.[^\.]*$//g`
	Dir=${outdir}/${Name2}


	#passしたハプロタイプを処理
	for x in `cat ${Dir}/${Name2}_hap_bam.list`
	do
		w=`echo ${x} | sed s/.sorted.bam//g`
		samtools index ${x}

		#SNP call bcftools
		bcftools mpileup -Ou -f ${REF} -a DP,AD,INFO/AD ${x} | bcftools call -v -m -Ov -f GQ -o ${w}_bcftools_raw.vcf &
		
		#SV call
		sniffles -s 5 -m ${x} -v ${w}.sv.raw.vcf &

		wait
		#bcftools
		echo "${outdir}/analysis/`basename ${w}.bam` `basename ${w}.bam`" > ${Dir}/sample_name.txt
		bcftools reheader -s ${Dir}/sample_name.txt ${w}_bcftools_raw.vcf -o ${w}_bcftools.vcf
	
		##SV callでIMPRECISEを除去,それぞれの結果をgz+index
		echo "${w}.sorted.bam `basename ${w}.bam`" | awk '{ if($2 != 2){print $0} }' > ${Dir}/sample_name.txt
		bcftools reheader -s ${Dir}/sample_name.txt ${w}.sv.raw.vcf -o ${w}.sv.raw2.vcf
		grep -v "IMPRECISE" ${w}.sv.raw2.vcf > ${w}.sv.vcf
		
		rm ${w}.sv.raw.vcf ${w}.sv.raw2.vcf
		
		bgzip ${w}_bcftools.vcf -f 
		bcftools index ${w}_bcftools.vcf.gz
		
		bgzip ${w}.sv.vcf -f
		bcftools index ${w}.sv.vcf.gz
	
	done
done
