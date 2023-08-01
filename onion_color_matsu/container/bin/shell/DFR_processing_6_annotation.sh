#!/bin/bash
##singularity shell snpEff_onion_DFR_ANS.sifで使用

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

#DFRのsnpEffアノテーション
for REF in `cat ${REF_list}`
do

        Name=`basename ${REF}`
        if [ ${Name} = "DFR_OGBv3.0_Scaffold003_T1280A.fa" ]; then

                Name2=`echo ${Name} | sed s/\.[^\.]*$//g`
                Dir=${outdir}/${Name2}

                for Sample in `cat ${Dir}/snp_indel_sv.list`
                do

                        Sample_name=`basename ${Sample}`
                        Sample_out=`echo ${Sample_name} | sed s/.vcf//g`

                        echo ${Sample} ${Sample_name} ${Sample_out}
                        cat ${Sample} | sed s/19-WGS-ACE.3.0_Scaffold003/19-WGS-ACE30_Scaffold003/g > ${Dir}/${Sample_out}.rename.vcf
                        snpEff Onion_DFR ${Dir}/${Sample_out}.rename.vcf > ${Dir}/${Sample_out}_ann.vcf
                        echo "${Dir}/${Sample_out}_ann.vcf" >> ${outdir}/DFR_OGBv3.0_Scaffold003_T1280A.ann.list
                done
        fi
done
