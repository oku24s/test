#!/bin/bash

program_path=`dirname $0`

mkdir ${program_path}/container/sif

##onion_color.sifの作成------------------------------------------------------------------------------------------------
#defファイルの書き出し
echo "
BootStrap: docker-archive
From: ${program_path}/container/base_images/continuumio3.tar

%post
        python --version
        conda --version
        conda install -c bioconda -c conda-forge whatshap=1.1
        conda install -c bioconda bcftools
        conda install -c bioconda samtools
        conda install -c bioconda minimap2
        conda install -c bioconda seqkit
        conda install htslib
        conda install -c bioconda sniffles
        conda install python-Levenshtein
        conda install -c bioconda biopython
        conda install -c bioconda pysam
	ln -s /opt/conda/lib/libcrypto.so.1.1 /opt/conda/lib/libcrypto.so.1.0.0

%files
        ${program_path}/container/bin/python/haplotype_splitter.py /usr/local/bin
        ${program_path}/container/bin/python/ONT_cigar_filter.py /usr/local/bin
        ${program_path}/container/bin/python/ONT_demultiplex_trimming_old_complete.py /usr/local/bin
        ${program_path}/container/bin/python/ONT_demultiplex_trimming_old.py /usr/local/bin
        ${program_path}/container/bin/python/ONT_demultiplex_trimming.v4.py /usr/local/bin
        ${program_path}/container/bin/python/ONT_processing_7_output_ANS.py /usr/local/bin
        ${program_path}/container/bin/python/ONT_processing_7_output_DFR.py /usr/local/bin
        ${program_path}/container/bin/python/ONT_processing_9_output_ANS_fin.py /usr/local/bin
        ${program_path}/container/bin/python/ONT_processing_9_output_DFR_fin.py /usr/local/bin
        ${program_path}/container/bin/shell/ONT_processing_0_demultiplex_trimming.sh /usr/local/bin
        ${program_path}/container/bin/shell/ONT_processing_1_reference_check.DFR.sh /usr/local/bin
        ${program_path}/container/bin/shell/ONT_processing_1_reference_check.ANS.sh /usr/local/bin
        ${program_path}/container/bin/shell/ONT_processing_2_haplotype_split.ANS.sh /usr/local/bin
        ${program_path}/container/bin/shell/ONT_processing_2_haplotype_split.DFR.sh /usr/local/bin
        ${program_path}/container/bin/shell/ONT_processing_3_minimap2.sh /usr/local/bin
        ${program_path}/container/bin/shell/ONT_processing_4_bcftools_sniffles.sh /usr/local/bin
        ${program_path}/container/bin/shell/ONT_processing_5_merge.sh /usr/local/bin
        ${program_path}/container/bin/shell/ANS_ONT_processing_6_merge_files.sh /usr/local/bin
        ${program_path}/container/bin/shell/DFR_ONT_processing_6_merge_files.sh /usr/local/bin
        ${program_path}/container/bin/shell/ONT_processing_8_output_ANS_fin.sh /usr/local/bin
        ${program_path}/container/bin/shell/ONT_processing_8_output_DFR_fin.sh /usr/local/bin
        ${program_path}/container/lib/ANS_ref /usr/local/lib
        ${program_path}/container/lib/DFR_ref /usr/local/lib
        ${program_path}/container/lib/Allele_dict /usr/local/lib

%environment
        export LC_ALL=C.UTF-8" > ${program_path}/container/def/onion_color.def

#sifファイルの作成
singularity build --fakeroot ${program_path}/container/sif/onion_color.sif ${program_path}/container/def/onion_color.def




##pepper_margin_deepvariantのコンテナ作成-------------------------------------------------------------------------------
#defファイルの書き出し
echo "
#BootStrap: docker
#From: kishwars/pepper_deepvariant:r0.8
BootStrap: docker-archive
From: ${program_path}/container/base_images/pepper_margin_deepvariant.tar

%files
        ${program_path}/container/bin/shell/ONT_processing_4_deepvariant.sh /usr/local/bin
        ${program_path}/container/lib/DFR_ref /usr/local/lib
        ${program_path}/container/lib/ANS_ref /usr/local/lib

%environment
        export LC_ALL=C" > ${program_path}/container/def/pepper_margin_deepvariant.def

#sifファイルの作成
singularity build --fakeroot ${program_path}/container/sif/pepper_margin_deepvariant.sif ${program_path}/container/def/pepper_margin_deepvariant.def


##medakaのコンテナ作成--------------------------------------------------------------------------------------------------
#defファイルの書き出し
echo "
#BootStrap: docker
#From: nanozoo/medaka:1.7.2--aa54076
BootStrap: docker-archive
From: ${program_path}/container/base_images/medaka.tar

%files
        ${program_path}/container/bin/shell/ONT_processing_4_medaka.sh /usr/local/bin
        ${program_path}/container/lib/ANS_ref /usr/local/lib
        ${program_path}/container/lib/DFR_ref /usr/local/lib

%environment
        export LC_ALL=C" > ${program_path}/container/def/medaka.def

#sifファイルの作成
singularity build --fakeroot ${program_path}/container/sif/medaka.sif ${program_path}/container/def/medaka.def



##SnpEffのコンテナ作成----------------------------------------------------------------------------------------------------------------------
#defファイルの書き出し
echo "
#BootStrap: docker
#From: nfcore/snpeff:5.1.R64-1-1
BootStrap: docker-archive
From: ${program_path}/container/base_images/snpeff.tar

%post
        #snpeffでデータベース構築
        snpEff build -gff3 -v Onion_DFR -noCheckProtein -noCheckCds
        snpEff build -gff3 -v Onion_ANS_L -noCheckProtein -noCheckCds
        snpEff build -gff3 -v Onion_ANS_h1 -noCheckProtein -noCheckCds

%files
	#スクリプトをコンテナにコピー
       	${program_path}/container/bin/shell/ANS_processing_6_annotation.sh /usr/local/bin
        ${program_path}/container/bin/shell/DFR_processing_6_annotation.sh /usr/local/bin
	#referenceをコンテナにコピー
        ${program_path}/container/lib/ANS_ref /usr/local/lib
        ${program_path}/container/lib/DFR_ref /usr/local/lib
        ${program_path}/container/lib/Allele_dict /usr/local/lib
	#データベース構築用データをコンテナにコピー
	${program_path}/container/snpEff/Onion_DFR /opt/conda/envs/nf-core-snpeff-5.1/share/snpeff-5.1-2/data
	${program_path}/container/snpEff/Onion_ANS_L /opt/conda/envs/nf-core-snpeff-5.1/share/snpeff-5.1-2/data
	${program_path}/container/snpEff/Onion_ANS_h1 /opt/conda/envs/nf-core-snpeff-5.1/share/snpeff-5.1-2/data
	${program_path}/container/snpEff/snpEff.config /opt/conda/envs/nf-core-snpeff-5.1/share/snpeff-5.1-2

%environment
        export LC_ALL=C.UTF-8" > ${program_path}/container/def/snpeff_onion.def

#sifファイルの作成
singularity build --fakeroot ${program_path}/container/sif/snpeff_onion.sif ${program_path}/container/def/snpeff_onion.def
