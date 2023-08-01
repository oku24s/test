#!/bin/bash
##singularity shell deepvariant_pepper.sifで使用

usage (){
        cat <<EOM
Usage: $(basename "$0") [OPTION]...
        -h              Display this help
        -o[str]         full path of output directory (~/analysis_name)
        -r[str]         full path of reference list (~/reference.list)
        -m[str]         deepvariant model ({ont_r9_guppy5_sup, ont_r10_q20})
EOM

        exit 2
}

while getopts ":o:r:m:h" optKey
do
        case ${optKey} in
                o)
                        outdir=${OPTARG}
                        ;;
                r)
                        REF_list=${OPTARG}
                        ;;
                m)
                        model=${OPTARG}
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
        echo "Reference list is necessary!!"
        exit
else
        echo "-r (reference list) = ${REF_list}"
fi

char3=`echo ${model} | cut -c 1`
if [ -z "${model}" ] || [ "${char3}" = '-' ]
then
        echo "model is set as ont_r9_guppy5_sup"
        model=ont_r9_guppy5_sup
else
        echo "-m (deepvariant model) = ${model}"
fi


for REF in `cat ${REF_list}`
do
        Name=`basename ${REF}`
        Name2=`echo ${Name} | sed s/.fa//g`
        Dir=${outdir}/${Name2}


        for x in `cat ${Dir}/${Name2}_hap_bam.list`
        do
                w=`echo ${x} | sed s/.sorted.bam//g`
                samtools index ${x}

                z=`basename ${w}`
                echo "Sample ${z}" > ${Dir}/sample_name.txt
		
		#pepper_margin_deepvariantを実行
                run_pepper_margin_deepvariant call_variant -f ${REF} -b ${x} -o ${Dir} -t 1 --${model}
                if [ -e ${Dir}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz ]; then
                        echo "reheader!! ${w}"
                        bcftools reheader -s ${Dir}/sample_name.txt ${Dir}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz -o ${Dir}/${z}_pepper_variant.vcf.gz

                        bcftools index ${Dir}/${z}_pepper_variant.vcf.gz

                        rm ${Dir}/PEPPER_VARIANT_FULL.vcf.gz ${Dir}/PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz ${Dir}/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz
                        rm ${Dir}/PEPPER_VARIANT_FULL.vcf.gz.tbi ${Dir}/PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz.tbi ${Dir}/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz.tbi
                        rm ${Dir}/PHASED.PEPPER_MARGIN.chunks.csv ${Dir}/PHASED.PEPPER_MARGIN.haplotagged.bam ${Dir}/PHASED.PEPPER_MARGIN.haplotagged.bam.bai
                        rm -r ${Dir}/logs ${Dir}/intermediate_files ${Dir}/dv_intermediate_outputs
                        rm ${Dir}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz ${Dir}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz.tbi ${Dir}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.visual_report.html
                else
                        bcftools reheader -s ${Dir}/sample_name.txt ${Dir}/intermediate_files/PEPPER_VARIANT_FULL.vcf.gz -o ${Dir}/${z}_pepper_variant.vcf.gz
                        bcftools index ${Dir}/${z}_pepper_variant.vcf.gz

                        rm ${Dir}/PEPPER_VARIANT_FULL.vcf.gz ${Dir}/PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz ${Dir}/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz
                        rm ${Dir}/PEPPER_VARIANT_FULL.vcf.gz.tbi ${Dir}/PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz.tbi ${Dir}/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz.tbi
                        rm ${Dir}/PHASED.PEPPER_MARGIN.chunks.csv ${Dir}/PHASED.PEPPER_MARGIN.haplotagged.bam ${Dir}/PHASED.PEPPER_MARGIN.haplotagged.bam.bai
                        rm -r ${Dir}/logs ${Dir}/intermediate_files ${Dir}/dv_intermediate_outputs
                        rm ${Dir}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz ${Dir}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz.tbi ${Dir}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.visual_report.html
                fi

        done

done
