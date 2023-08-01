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

awk -v 'OFS=\t' '{ print $1,$2,$3,1}' < ${outdir}/active_sample.tsv | sed "s/\t/,/g" | sed "s/_H1//g" | sed "s/_H2//g" | sed "s/_homo//g" >> ${outdir}/result.raw.csv
awk -v 'OFS=\t' '{ print $1,$2,$3,0}' < ${outdir}/inactive_sample.tsv | sed "s/\t/,/g" | sed "s/_H1//g" | sed "s/_H2//g" | sed "s/_homo//g" >> ${outdir}/result.raw.csv


for y in `cat ${outdir}/TRN_DTP_LTR.list`
do
	z=`basename ${y} | sed "s/.match.bam//g"`
	echo "${z},TRN,inactive,0" >> ${outdir}/result.raw.csv
done



#result.csvの重複削除
sort -u ${outdir}/result.raw.csv | uniq -u > ${outdir}/result.csv
