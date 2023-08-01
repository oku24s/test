#!/bin/bash

##singularity shell ONT_ampliconseq.sifで使用

usage (){
	cat <<EOM
Usage: $(basename "$0") [OPTION]...
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


for REF in `cat ${REF_list}`
do
	Name=`basename ${REF}`
	Name2=`echo ${Name} | sed s/\.[^\.]*$//g`
	for x in `cat ${outdir}/${Name2}/${Name2}.match_bam.list`
	do
		w=`echo ${x} | sed s/.match.bam//g`
		#.match.bamのリード数が5未満の場合、解析を進めない
		read_count=`samtools view -F 0x04 -c ${x}`
		if [ ${read_count} -lt 5 ];then
			echo "the number of reads is ${read_count}"
			continue
		fi

		#bcftoolsで1回目の変異コール
		samtools index ${x}

		bcftools mpileup -Ou -f ${REF} -a DP,AD,INFO/AD ${x} | bcftools call -v -m -f GQ -o ${w}_bcftools.raw.vcf
		#bcftools index ${w}_bcftools.raw.vcf.gz -f
		rename1=`grep "#" ${w}_bcftools.raw.vcf | grep -v "##" | awk '{print $10}'`
		rename2=`basename ${rename1}`
		echo ${rename1} ${rename2} > ${outdir}/${Name2}/sample_name.txt
		bcftools reheader -s ${outdir}/${Name2}/sample_name.txt ${w}_bcftools.raw.vcf -o ${w}_bcftools.vcf
		bgzip ${w}_bcftools.raw.vcf
		bcftools index ${w}_bcftools.raw.vcf.gz
		
		sniffles -s 5 -m ${x} -v ${w}.sv.raw.vcf
		echo "${x} `basename ${x} | sed s/\.match\.bam//g`" > ${outdir}/${Name2}/sample_name.txt
		bcftools reheader -s ${outdir}/${Name2}/sample_name.txt ${w}.sv.raw.vcf -o ${w}.sv.raw2.vcf
		grep -v "IMPRECISE" ${w}.sv.raw2.vcf > ${w}.sv.vcf
		bgzip ${w}.sv.vcf
		bcftools index ${w}.sv.vcf.gz
		
		bgzip ${w}_bcftools.vcf
		bcftools index ${w}_bcftools.vcf.gz

		bcftools concat ${w}_bcftools.vcf.gz ${w}.sv.vcf.gz -o ${w}_bcftools.sv.vcf.gz -Oz
		bcftools index ${w}_bcftools.sv.vcf
		bgzip ${w}_bcftools.vcf -f

		#要検討箇所　GQのフィルタリングする？リードが少ない場合、ヘテロについてGQが十分に高く出るか？

		#phasingと split bam
		whatshap phase -o ${w}.phased.vcf --reference ${REF} ${w}_bcftools.raw.vcf.gz ${w}.match.bam
		bgzip ${w}.phased.vcf -f
		bcftools index ${w}.phased.vcf.gz

		whatshap haplotag --reference ${REF} --output-haplotag-list ${w}.haplotagged.list\
			-o ${w}.haplotagged.bam ${w}.phased.vcf.gz ${w}.match.bam --ignore-read-groups

		samtools index ${w}.haplotagged.bam

		whatshap split --output-h1 ${w}_H1.bam --output-h2 ${w}_H2.bam ${w}.haplotagged.bam ${w}.haplotagged.list

		#要検討箇所　bamからfastqに戻す。Hap1とHap2が両方とも十分にあればsorted.bamは削除。Hap1とHap2が両方ともなければホモとしてmatch.bamを使う
		#Hap1だけ or Hap2だけが十分なデータ量となった場合はどうする？
		#⇒とりあえず、Hap1.bam or Hap2.bamとして出力されたものをホモとして扱うことにして以下を記載。Hap1.bam or Hap2.bamよりもmatch.bamを使ったほうがいいかも？

		Hap1=`samtools view -c ${w}_H1.bam`
		Hap2=`samtools view -c ${w}_H2.bam`

		result=""

		if [ ${Hap1} -ge 2 ] && [ ${Hap2} -ge 2 ];then
			samtools fastq ${w}_H1.bam > ${w}_H1.fastq
			samtools fastq ${w}_H2.bam > ${w}_H2.fastq
			if [ `wc -l ${w}_H1.fastq` -gt 1 ];then
                        	echo ${w}_H1.fastq >> ${outdir}/analysis/fastq_split.list
			fi
			if [ `wc -l ${w}_H2.fastq` ];then
                        	echo ${w}_H2.fastq >> ${outdir}/analysis/fastq_split.list
			fi

		elif [ ${Hap1} -ge 2 ] && [ ${Hap2} -lt 2 ];then
			#ホモと判定する前にhaplotype_splitter.pyを実行
			gunzip ${w}_bcftools.sv.vcf.gz
			result=`haplotype_splitter.py -v ${w}_bcftools.sv.vcf -b ${x} -t DFR`
			echo ${result}
			if [ "${result}" = "homo" ] || [ "${result}" = "multi" ];then
				samtools fastq ${w}_H1.bam > ${w}_homo.fastq
				echo ${w}_homo.fastq >> ${outdir}/analysis/fastq_split.list
			else
				samtools fastq ${w}_H1.bam > ${w}_H1.fastq
                        	samtools fastq ${w}_H2.bam > ${w}_H2.fastq
				if [ `wc -l ${w}_H1.fastq` -gt 1 ];then
                        		echo ${w}_H1.fastq >> ${outdir}/analysis/fastq_split.list
				fi
				if [ `wc -l ${w}_H2.fastq` ];then
                        		echo ${w}_H2.fastq >> ${outdir}/analysis/fastq_split.list
				fi
			fi

		elif [ ${Hap1} -lt 2 ] && [ ${Hap2} -ge 2 ];then
			#ホモと判定する前にhaplotype_splitter.pyを実行
			gunzip ${w}_bcftools.sv.vcf.gz
			result=`haplotype_splitter.py -v ${w}_bcftools.sv.vcf -b ${x} -t DFR`
			echo ${result}
			if [ "${result}" = "homo" ] || [ "${result}" = "multi" ];then
				samtools fastq ${w}_H2.bam > ${w}_homo.fastq
				echo ${w}_homo.fastq >> ${outdir}/analysis/fastq_split.list
			else
				samtools fastq ${w}_H1.bam > ${w}_H1.fastq
                        	samtools fastq ${w}_H2.bam > ${w}_H2.fastq
				if [ `wc -l ${w}_H1.fastq` -gt 1 ];then
                        		echo ${w}_H1.fastq >> ${outdir}/analysis/fastq_split.list
				fi
				if [ `wc -l ${w}_H2.fastq` ];then
                        		echo ${w}_H2.fastq >> ${outdir}/analysis/fastq_split.list
				fi
			fi

		elif [ ${Hap1} -lt 2 ] && [ ${Hap2} -lt 2 ];then
			#ホモと判定する前にhaplotype_splitter.pyを実行
			gunzip ${w}_bcftools.sv.vcf.gz
			result=`haplotype_splitter.py -v ${w}_bcftools.sv.vcf -b ${x} -t DFR`
			echo ${result}
			if [ "${result}" = "homo" ] || [ "${result}" = "multi" ];then
				samtools fastq ${w}.match.bam > ${w}_homo.fastq
				echo ${w}_homo.fastq >> ${outdir}/analysis/fastq_split.list
			else
				samtools fastq ${w}_H1.bam > ${w}_H1.fastq
                        	samtools fastq ${w}_H2.bam > ${w}_H2.fastq
				if [ `wc -l ${w}_H1.fastq` -gt 1 ];then
                        		echo ${w}_H1.fastq >> ${outdir}/analysis/fastq_split.list
				fi
				if [ `wc -l ${w}_H2.fastq` ];then
                        		echo ${w}_H2.fastq >> ${outdir}/analysis/fastq_split.list
				fi
			fi
		fi

		#haplotagged.bam, haplotagged.bam.bai, phased.vcf.gz, phased.vcf.gz.csi, H1.bam, H2.bam, haplotagged.list, bcftools.vcf.gz, bcftools.vcf.gz.csi
		#rm ${w}.haplotagged.bam ${w}.haplotagged.bam.bai ${w}.phased.vcf.gz ${w}.phased.vcf.gz.csi ${w}_H1.bam ${w}_H2.bam ${w}.haplotagged.list #${w}_bcftools.vcf.gz ${w}_bcftools.vcf.gz.csi

	done
done
