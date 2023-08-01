#!/bin/bash
##singularity shell ONT_ampliconseq.sifで使用

usage (){
        cat <<EOM
Usage: $(basename "$0") [OPTION]...
        -h              Display this help
        -o[str]         full path of output directory (~/analysis_name)
        -r[str]         full path of reference list (~/reference.list)
	-n[int]		threads (= split number)
EOM

        exit 2
}

while getopts ":o:r:n:h" optKey
do
        case ${optKey} in
                o)
                        outdir=${OPTARG}
                        ;;
                r)
                        REF_list=${OPTARG}
                        ;;
		n)
			split_num=${OPTARG}
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



for x in `cat ${REF_list}`
do
        #行先のディレクトリを準備
        ref_name=`basename ${x}| sed s/.fa//g`
        mkdir ${outdir}/${ref_name}

        #アノテーション後のvcfファイルとbamファイルを移動
        for y in `seq 1 ${split_num}`
        do
                mv ${outdir}/${y}/${ref_name}/*.vcf* ${outdir}/${ref_name}
                mv ${outdir}/${y}/${ref_name}/*.bam ${outdir}/${ref_name}
                mv ${outdir}/${y}/${ref_name}/*.fastq ${outdir}/${ref_name}
                cat ${outdir}/${y}/${ref_name}/snp_indel_sv.list >> ${outdir}/${ref_name}/snp_indel_sv.list
        done

	ls ${outdir}/${ref_name}/*match.bam > ${outdir}/${ref_name}/${ref_name}.match_bam.list
done

#DFR_OGBv3.0_Scaffold003_T1280A.ann.listの中身を統合
ls ${outdir}/DFR_OGBv3.0_Scaffold003_T1280A/*ann.vcf > ${outdir}/DFR_OGBv3.0_Scaffold003_T1280A.ann.list

for y in `seq 1 ${split_num}`
do
	rm -r ${outdir}/${y}
done


#DFR_TRN,DTP,LTRをもつサンプルのリストをつくる
cat ${outdir}/DFR-TRN_amplicon/DFR-TRN_amplicon.match_bam.list >> ${outdir}/TRN_DTP_LTR.list
cat ${outdir}/DFR-ADTP_amplicon/DFR-ADTP_amplicon.match_bam.list >> ${outdir}/TRN_DTP_LTR.list
cat ${outdir}/DFR-ALTR_amplicon/DFR-ALTR_amplicon.match_bam.list >> ${outdir}/TRN_DTP_LTR.list


for x in `cat ${outdir}/analysis/fastq.list`
do
	sample_name=`basename ${x} | sed "s/.fastq//g"`
	echo ${sample_name}
	Full_count=`cat ${outdir}/DFR_OGBv3.0_Scaffold003_T1280A.ann.list | grep ${sample_name} | wc -l`
	TRN_DTP_LTR_count=`cat ${outdir}/TRN_DTP_LTR.list | grep ${sample_name} | wc -l`
	
	#Fullのアリルが２種類とTRNが検出された場合、
	if [ ${Full_count} -ge 2 ] && [ ${TRN_DTP_LTR_count} -ge 1 ]; then
		Full_hap1=`samtools view -F 0x04 -c ${outdir}/DFR_OGBv3.0_Scaffold003_T1280A/${sample_name}_H1.sorted.bam`
		Full_hap2=`samtools view -F 0x04 -c ${outdir}/DFR_OGBv3.0_Scaffold003_T1280A/${sample_name}_H2.sorted.bam`
		TRN=`samtools view -F 0x04 -c ${outdir}/DFR-TRN_amplicon/${sample_name}.match.bam`
		
		#リード数が最下位のアリルを削除する
		##Full_hap1がFull_hap2より少なく、Full_hap1がTRNより少ない場合
		if [ ${Full_hap1} -lt ${Full_hap2} ] && [ ${Full_hap1} -lt ${TRN} ];then
			#Full_hap1をDFR_OGBv3.0_Scaffold003_T1280A.ann.listから削除
			sed -i -e "/${sample_name}_H1/d" ${outdir}/DFR_OGBv3.0_Scaffold003_T1280A.ann.list
		
		##Full_hap2がFull_hap1より少なく、Full_hap2がTRNより少ない場合
		elif [ ${Full_hap2} -lt ${Full_hap1} ] && [ ${Full_hap2} -lt ${TRN} ];then
			#Full_hap2をDFR_OGBv3.0_Scaffold003_T1280A.ann.listから削除
			sed -i -e "/${sample_name}_H2/d" ${outdir}/DFR_OGBv3.0_Scaffold003_T1280A.ann.list
		
		##TRNがFull_hap1より少なく、TRNがFull_hap2より少ない場合
		else
			#TRNをTRN_DTP_LTR.listから削除
			sed -i -e "/${sample_name}/d" ${outdir}/TRN_DTP_LTR.list

		fi
	fi
done
