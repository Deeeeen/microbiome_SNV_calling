#!/bin/bash

#### !!!! ####
#### you may want to link the software call commands to 
#### where you install the corresponding software
#### !!!! ####

#### software version used ####
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

	for sample in ${samples[@]}
	do
		echo -e "\n----- Alignment Statistics -----"
		samtools stats \
			${dir}/${ref}/bam/${sample}_sorted.bam \
			> ${dir}/${ref}/bam/${sample}_stats

		samtools coverage \
			${dir}/${ref}/bam/${sample}_sorted.bam \
			> ${dir}/${ref}/bam/${sample}_coverage

		samtools depth \
			-aa ${dir}/${ref}/bam/${sample}_sorted.bam \
			> ${dir}/${ref}/bam/${sample}_depth

		echo -e "\n----- SNV Calling Statistics-----"
		bcftools stats \
			${dir}/${ref}/vcf/${sample}_filtered.vcf.gz > ${dir}/${ref}/vcf/${sample}_filtered_stats
		
		bcftools query \
			-f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' \
			${dir}/${ref}/vcf/${sample}_filtered.vcf.gz > ${dir}/${ref}/vcf/${sample}_filtered_GT

	done
	bcftools stats \
			${dir}/${ref}/BH1206_${ref}_filtered.vcf.gz > ${dir}/${ref}/BH1206_${ref}_filtered_stats
done

for sample in ${samples[@]}
do
	echo -e "${sample} "
	zcat ${data_dir}/${sample}.fq.gz | echo $((`wc -l`/4))
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Alignment and SNV Calling Statistics Time elapsed: $(( $END - $START )) seconds"
