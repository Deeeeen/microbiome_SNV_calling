#!/bin/bash

#### !!!! ####
#### you may want to link the software call commands to 
#### where you install the corresponding software
#### !!!! ####

#### software version used ####
#### bowtie2-2.4.3 
#### samtools-1.12
#### bcftools-1.12

#### directory where you want to keep the output data
dir=$1
#### directory where you keep the reference genomes
ref_dir=$2
#### directory where you keep the fastq sequence files
data_dir=$3


START=$(date +%s)
#### TODO: change to read from ${ref_dir}
refs=(NZ_CP022479.1 NZ_CP023819.1 NZ_CP030777.1 NZ_CP048437.1)
#### TODO: change to read from ${data_dir}
samples=(B1_A10 B1_A11 B1_A12 B1_A13 B1_A14 B1_A15 B1_A16 B1_A17 B1_A18 B1_A19 B1_A1 B1_A20 B1_A21 B1_A2 B1_A3 B1_A4 B1_A5 B1_A6 B1_A7 B1_A8 B1_A9 B1_B10 B1_B11 B1_B12 B1_B13 B1_B14 B1_B15 B1_B16 B1_B17 B1_B18 B1_B19 B1_B1 B1_B20 B1_B21 B1_B2 B1_B3 B1_B4 B1_B5 B1_B6 B1_B7 B1_B8 B1_B9)

for ref in ${refs[@]}
do
	echo -e "\n----- ${ref} -----"
	
	echo -e "\n----- Index Fasta file -----"
	bowtie2-build \
		${ref_dir}/${ref}.fasta ${ref_dir}/${ref}.index

	mkdir ${dir}/${ref}
	mkdir ${dir}/${ref}/bam
	mkdir ${dir}/${ref}/vcf

	for sample in ${samples[@]}
	do
		echo -e "\n----- bowtie2 Alignment -----"
		bowtie2 -p 8 \
			-x ${ref_dir}/${ref}.index --no-mixed --very-sensitive \
			--n-ceil 0,0.01 -U ${data_dir}/${sample}.fq.gz -S ${dir}/${ref}/bam/${sample}.sam

		echo -e "\n----- SAM -> BAM, sort BAM -----"
		samtools view \
			-b ${dir}/${ref}/bam/${sample}.sam -o ${dir}/${ref}/bam/${sample}.bam

		samtools sort \
			-M -m 10G -o ${dir}/${ref}/bam/${sample}_sorted.bam ${dir}/${ref}/bam/${sample}.bam

		echo -e "\n----- SNV Calling -----"
		bcftools mpileup \
			-C 50 -m 3 -F 0.0002 -f ${ref_dir}/${ref}.fasta -Ou ${dir}/${ref}/bam/${sample}_sorted.bam | \
		bcftools call -c -v \
		   --ploidy 1 -Ov -o ${dir}/${ref}/vcf/${sample}.vcf

		rm ${dir}/${ref}/bam/${sample}.bam
		rm ${dir}/${ref}/bam/${sample}.sam

		echo -e "\n----- variants filtering -----"
		vcfutils.pl varFilter \
			-d 100 ${dir}/${ref}/vcf/${sample}.vcf > ${dir}/${ref}/vcf/${sample}_filtered_temp.vcf
		
		bgzip ${dir}/${ref}/vcf/${sample}_filtered_temp.vcf

		echo -e "\n----- index vcf -----"
		bcftools index \
			-t ${dir}/${ref}/vcf/${sample}_filtered_temp.vcf.gz \
			-o ${dir}/${ref}/vcf/${sample}_filtered_temp.vcf.gz.tbi

		echo -e "\n----- filter vcfs QUAL>=60 -----"
		/projects/ps-palmer/software/local/src/bcftools-1.12/bcftools view \
			-i '%QUAL>=60' -Oz -o ${dir}/${ref}/vcf/${sample}_filtered.vcf.gz ${dir}/${ref}/vcf/${sample}_filtered_temp.vcf.gz
		
		echo -e "\n----- index vcf -----"
		bcftools index \
			-t ${dir}/${ref}/vcf/${sample}_filtered.vcf.gz \
			-o ${dir}/${ref}/vcf/${sample}_filtered.vcf.gz.tbi
		
		rm ${dir}/${ref}/vcf/${sample}_filtered_temp.vcf.gz 
		rm ${dir}/${ref}/vcf/${sample}_filtered_temp.vcf.gz.tbi
	done

	ls ${dir}/${ref}/vcf/*_filtered.vcf.gz > ${dir}/${ref}/vcf_list

	echo -e "\n----- merge vcfs -----"
	bcftools merge \
		-m all -l ${dir}/${ref}/vcf_list -Oz -o ${dir}/${ref}/BH1206_${ref}.vcf.gz

	echo -e "\n----- index vcf -----"
	bcftools index \
		-t ${dir}/${ref}/BH1206_${ref}.vcf.gz \
		-o ${dir}/${ref}/BH1206_${ref}.vcf.gz.tbi
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "SNV Calling Time elapsed: $(( $END - $START )) seconds"
