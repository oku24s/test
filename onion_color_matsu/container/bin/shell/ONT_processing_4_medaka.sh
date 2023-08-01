#!/bin/bash
##singularity shell medaka.sifで使用

usage (){
        cat <<EOM
Usage: $(basename "$0")[OPTION]...
        -h              Display this help
        -o[str]         full path of output directory (~/analuysis_name)
        -r[str]         fullpath of reference list (~/reference.list)
        -m[str]         medaka model
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
                        medaka_model=${OPTARG}
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


        for x in `cat ${Dir}/${Name2}_hap_bam.list`
        do
                w=`echo ${x} | sed s/.sorted.bam//g`
                samtools index ${x}
                medaka consensus ${x} ${w}_medaka.hdf --model ${medaka_model}
                medaka variant ${REF} ${w}_medaka.hdf ${w}_medaka.raw.vcf

                echo "SAMPLE `basename ${w}`" > ${Dir}/sample_name.txt
                bcftools reheader -s ${Dir}/sample_name.txt ${w}_medaka.raw.vcf -o ${w}_medaka.vcf

                bgzip ${w}_medaka.vcf -f
                bcftools index ${w}_medaka.vcf.gz

                rm ${w}_medaka.raw.vcf ${w}_medaka.hdf


        done

done
