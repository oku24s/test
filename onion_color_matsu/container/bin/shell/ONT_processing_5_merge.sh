#!/bin/bash
##singularity shell ONT_ampliconseq.sifで使用

usage (){
        cat <<EOM
Usage: $(basename "$0") [OPTION]...
        -h              Display this help
        -o[str]         full path of output directory (~/analysis_name)
        -r[str]         full path of reference list (~/reference.list)
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


#passしたハプロタイプを処理
for REF in `cat ${REF_list}`
do
	Name=`basename ${REF}`
	Name2=`echo ${Name} | sed s/\.[^\.]*$//g`
	Dir=${outdir}/${Name2}
	rm ${Dir}/${Name2}_merge.list
	rm ${Dir}/snp_indel_sv.list

	for x in `cat ${Dir}/${Name2}_hap_bam.list`
	do
		w=`echo ${x} | sed s/.sorted.bam//g`
	
		# normalize indels
		bcftools norm --fasta-ref ${REF} -m - -o ${w}_bcftools_norm.vcf ${w}_bcftools.vcf.gz
		bcftools norm --fasta-ref ${REF} -m - -o ${w}_medaka_norm.vcf ${w}_medaka.vcf.gz
		bcftools norm --fasta-ref ${REF} -m - -o ${w}_deepvariant_norm.vcf ${w}_pepper_variant.vcf.gz
	
		#deepvariantのみサンプル名の変更が要らないため、bgzip + index
		bgzip ${w}_deepvariant_norm.vcf -f
		bcftools index ${w}_deepvariant_norm.vcf.gz


		y=${w}
	        echo ${y}
		h=`basename ${y}`
	        echo "`basename ${y}`.bam `basename ${w}`" > ${Dir}/sample_name.txt
	      	bcftools reheader -s ${Dir}/sample_name.txt -o ${w}_bcftools_norm_rename.vcf ${w}_bcftools_norm.vcf
		bgzip ${w}_bcftools_norm_rename.vcf -f
		bcftools index ${w}_bcftools_norm_rename.vcf.gz

		

	      	bcftools reheader -s ${Dir}/sample_name.txt -o ${w}_medaka_norm_rename.vcf ${w}_medaka_norm.vcf
		bgzip ${w}_medaka_norm_rename.vcf -f
		bcftools index ${w}_medaka_norm_rename.vcf.gz
		

		#rename前のvcfを削除
		rm ${w}_bcftools_norm.vcf ${w}_medaka_norm.vcf
	
		#medaka, bcftools, deepvariantをmerge. 2つ以上のソフトでコールされた変異を採用
		bcftools isec -n +2 ${w}_deepvariant_norm.vcf.gz ${w}_medaka_norm_rename.vcf.gz ${w}_bcftools_norm_rename.vcf.gz -Ov -o ${w}_snp_indel_norm.list
		zcat ${w}_bcftools_norm_rename.vcf.gz | grep "#" > ${w}_snp_indel_norm.vcf
		cat ${w}_snp_indel_norm.list | while read line
		do
			z=`echo ${line} | cut -d " " -f 5`
			pos=`echo ${line} | cut -d " " -f 2`
			echo ${z}
			if [ ${z} = 001 ] || [ ${z} = 101 ] || [ ${z} = 111 ] || [ ${z} = 011 ];then
				zcat ${w}_bcftools_norm_rename.vcf.gz | grep -v "#" | awk -v foo=${pos} '{if($2 == foo) { print $0}}'  >> ${w}_snp_indel_norm.vcf
			
			elif [ ${z} = 010 ] || [ ${z} = 110 ];then
				echo ${z}
				zcat ${w}_medaka_norm_rename.vcf.gz | grep -v "#" | awk -v foo=${pos} '{if($2 == foo) {print $0}}' >> ${w}_snp_indel_norm.vcf
			
			else
				zcat ${w}_deepvariant_norm.vcf.gz | grep -v "#" | awk -v foo=${pos} '{if($2 == foo) {print $0}}' >> ${w}_snp_indel_norm.vcf
			fi
		done

		#ANS_Lのpos:1252のdeletionはmedakaでのみ安定して検出できるため、例外的にmedakaの出力を転記する
		if `grep 1252 ${w}_snp_indel_norm.vcf | wc -l` == 0; then
			zcat ${w}_medaka_norm_rename.vcf.gz | grep -v "#" | awk '{if($2 == 1252 && $1 == "ANS-L_amplicon") {print $0}}' >> ${w}_snp_indel_norm.vcf
		fi

		#DFRのpos:1280のsnpはbcftoolsでのみ安定して検出できるため、例外的にbcftoolsの出力を転記する
		if `grep 1280 ${w}_snp_indel_norm.vcf | wc -l` == 0; then
			zcat ${w}_bcftools_norm_rename.vcf.gz | grep -v "#" | awk '{if($2 == 1280 && $1 == "DFR_OGBv30_Scaffold003_T1280A") {print $0}}' >> ${w}_snp_indel_norm.vcf
		fi


		cp ${w}_snp_indel_norm.vcf ${w}_snp_indel_SV.vcf
		zcat ${w}.sv.vcf.gz | grep -v "#" >> ${w}_snp_indel_SV.vcf
	
	

		echo "${w}_snp_indel_SV.vcf" >> ${Dir}/snp_indel_sv.list
	
	done

done
