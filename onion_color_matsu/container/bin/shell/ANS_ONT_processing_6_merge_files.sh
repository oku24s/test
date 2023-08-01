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
                mv ${outdir}/${y}/${ref_name}/*.bam* ${outdir}/${ref_name}
                mv ${outdir}/${y}/${ref_name}/*.fastq ${outdir}/${ref_name}
        done
	#.ann.listの中身を統合
	ls ${outdir}/${ref_name}/*ann.vcf > ${outdir}/${ref_name}.ann.list

	#.match_bam.listの統合
	ls ${outdir}/${x}/*match.bam > ${outdir}/${ref_name}.match_bam.list
done

for x in `seq 1 ${split_num}`
do
	rm -r ${outdir}/${x}
done




for x in `cat ${outdir}/analysis/fastq.list`
do
	sample_name=`basename ${x} | sed "s/.fastq//g"`
	echo ${sample_name}
	h1_count=`cat ${outdir}/ANS_h1_amplicon.ann.list | grep ${sample_name} | wc -l`
	L_count=`cat ${outdir}/ANS_L_amplicon.ann.list | grep ${sample_name} | wc -l`

	if [ ${h1_count} -ge 2 ] && [ ${L_count} -ge 2 ]; then
		h1_hap1=`samtools view -F 0x04 -c ${outdir}/ANS_h1_amplicon/${sample_name}_H1.sorted.bam`
		h1_hap2=`samtools view -F 0x04 -c ${outdir}/ANS_h1_amplicon/${sample_name}_H2.sorted.bam`
		if [ ${h1_hap1} -ge ${h1_hap2} ];then
			#h1_hap2をANS_h1_amplicon.ann.listから削除
			sed -i -e "/${sample_name}_H2/d" ${outdir}/ANS_h1_amplicon.ann.list
			echo ${sample_name}_h1_H2
		else
			#h1_hap2をANS_h1_amplicon.ann.listから削除
			sed -i -e "/${sample_name}_H1/d" ${outdir}/ANS_h1_amplicon.ann.list
			echo ${sample_name}_h1_H1
		fi


		L_hap1=`samtools view -F 0x04 -c ${outdir}/ANS_L_amplicon/${sample_name}_H1.sorted.bam`
		L_hap2=`samtools view -F 0x04 -c ${outdir}/ANS_L_amplicon/${sample_name}_H2.sorted.bam`
		if [ ${L_hap1} -ge ${L_hap2} ];then
			#h1_hap2をANS_L_amplicon.ann.listから削除
			sed -i -e "/${sample_name}_H2/d" ${outdir}/ANS_L_amplicon.ann.list
			echo ${sample_name}_L_H2
		else
			#h1_hap2をANS_L_amplicon.ann.listから削除
			sed -i -e "/${sample_name}_H1/d" ${outdir}/ANS_L_amplicon.ann.list
			echo ${sample_name}_L_H1
		fi

	elif [ ${h1_count} -ge 2 ] && [ ${L_count} -eq 1 ]; then
		h1_hap1=`samtools view -F 0x04 -c ${outdir}/ANS_h1_amplicon/${sample_name}_H1.sorted.bam`
		h1_hap2=`samtools view -F 0x04 -c ${outdir}/ANS_h1_amplicon/${sample_name}_H2.sorted.bam`
		if [ ${h1_hap1} -ge ${h1_hap2} ];then
			#h1_hap2をANS_h1_amplicon.ann.listから削除
			sed -i -e "/${sample_name}_H2/d" ${outdir}/ANS_h1_amplicon.ann.list
			echo ${sample_name}_h1_H2
		else
			#h1_hap2をANS_h1_amplicon.ann.listから削除
			sed -i -e "/${sample_name}_H1/d" ${outdir}/ANS_h1_amplicon.ann.list
			echo ${sample_name}_h1_H1
		fi
	

	elif [ ${h1_count} -eq 1 ] && [ ${L_count} -ge 2 ]; then
		L_hap1=`samtools view -F 0x04 -c ${outdir}/ANS_L_amplicon/${sample_name}_H1.sorted.bam`
		L_hap2=`samtools view -F 0x04 -c ${outdir}/ANS_L_amplicon/${sample_name}_H2.sorted.bam`
		if [ ${L_hap1} -ge ${L_hap2} ];then
			#h1_hap2をANS_L_amplicon.ann.listから削除
			sed -i -e "/${sample_name}_H2/d" ${outdir}/ANS_L_amplicon.ann.list
			echo ${sample_name}_L_H2
		else
			#h1_hap2をANS_L_amplicon.ann.listから削除
			sed -i -e "/${sample_name}_H1/d" ${outdir}/ANS_L_amplicon.ann.list
			echo ${sample_name}_L_H1
		fi
	fi
	echo ${h1_count} ${L_count}
	h1_count=0
	L_count=0
done
