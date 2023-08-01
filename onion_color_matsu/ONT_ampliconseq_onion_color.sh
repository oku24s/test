#!/bin/sh


echo "解析の対象を選択してください"
echo "タマネギDFRの解析：1"
echo "タマネギANSの解析：2"
read menu

echo "config fileを指定してください"
read conf

#config fileを読み込む
source ${conf}

#programのpathを取得
program_path=`dirname $0`
echo ${program_path}

if [ ${menu} = 1 ]; then
	
	echo "DFRの解析をします"
	qsub -pe threads ${thread} ${program_path}/script/DFR_ONT_ampliconseq_jobs.split.sh ${indir} ${outdir} ${thread} ${sample_list} ${debarcoding} ${DFR_Full_threshold} ${DFR_TRN_threshold} ${medaka_model} ${deepvariant_model} ${program_path}

elif [ ${menu} = 2 ]; then
	
	echo "ANSの解析をします"
	qsub -pe threads ${thread} ${program_path}/script/ANS_ONT_ampliconseq_jobs.split.sh ${indir} ${outdir} ${thread} ${sample_list} ${debarcoding} ${ANS_h1_threshold} ${ANS_L_threshold} ${medaka_model} ${deepvariant_model} ${program_path}

else
	
	echo "1か2を入力してください"

fi
