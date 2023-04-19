#!/usr/bin/bash

GATK_Package_jarFile_address='/media/ad/9117284696_AD/carrying/NGSguide/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar'
Varscan2_jarFile_address='VarScan.v2.3.9.jar'

Thread=12
Genome_FASTA='/media/ad/9117284696_AD/NGSneeds/Genome_FASTA/hg19/hg19.fa'
GATK_Known_site_db1='/media/ad/9117284696_AD/NGSneeds/Variant_Known_Databases/dbSNP_146_hg19/dbsnp_138.hg19.vcf'

primaryBam_Normal='/media/ad/9117284696_AD/SomaticVariantCalling/Exome_Norm/Exome_Normal.bam'
primaryBam_Tumor='/media/ad/9117284696_AD/SomaticVariantCalling/Exome_Tumor/Exome_Tumor.bam'

Varscan2_Output_path='/media/ad/9117284696_AD/SomaticVariantCalling/Varscan/exome/'


####BAM sorting
samtools sort -@ ${Thread} -o ${primaryBam_Normal}_Sorted ${primaryBam_Normal}
samtools sort -@ ${Thread} -o ${primaryBam_Tumor}_Sorted  ${primaryBam_Tumor} 

rm -f ${primaryBam_Normal}
rm -f ${primaryBam_Tumor}

#java -jar ${GATK_Package_jarFile_address} SortSam -I ${primaryBam_Normal} -O ${primaryBam_Normal}_Sorted -SO coordinate
#java -jar ${GATK_Package_jarFile_address} SortSam -I ${primaryBam_Tumor} -O ${primaryBam_Tumor}_Sorted -SO coordinate


####Read group adding
java -jar ${GATK_Package_jarFile_address} AddOrReplaceReadGroups -I ${primaryBam_Normal}_Sorted -O ${primaryBam_Normal}_Sorted_RG -LB WES-Somatic -PL Illumina -PU Unit1 -SM Normal
java -jar ${GATK_Package_jarFile_address} AddOrReplaceReadGroups -I ${primaryBam_Tumor}_Sorted -O ${primaryBam_Tumor}_Sorted_RG -LB WES-Somatic -PL Illumina -PU Unit1 -SM Tumor

rm -f ${primaryBam_Normal}_Sorted
rm -f ${primaryBam_Tumor}_Sorted


####Mark duplicated reads
java -jar ${GATK_Package_jarFile_address} MarkDuplicates -I ${primaryBam_Normal}_Sorted_RG -O ${primaryBam_Normal}_Sorted_RG_MD  -M ${primaryBam_Normal}_Sorted_RG_MD.stat
java -jar ${GATK_Package_jarFile_address} MarkDuplicates -I ${primaryBam_Tumor}_Sorted_RG -O ${primaryBam_Tumor}_Sorted_RG_MD -M ${primaryBam_Tumor}_Sorted_RG_MD.stat

rm -f ${primaryBam_Normal}_Sorted_RG
rm -f ${primaryBam_Tumor}_Sorted_RG


####Base Quality Score Recalibration
java -jar ${GATK_Package_jarFile_address} BaseRecalibrator -I ${primaryBam_Normal}_Sorted_RG_MD -O ${primaryBam_Normal}_Sorted_RG_MD_bqsr.table -R ${Genome_FASTA} --known-sites ${GATK_Known_site_db1}
java -jar ${GATK_Package_jarFile_address} ApplyBQSR -bqsr ${primaryBam_Normal}_Sorted_RG_MD_bqsr.table -I ${primaryBam_Normal}_Sorted_RG_MD -O ${primaryBam_Normal}_Sorted_RG_MD_BQSR.bam

java -jar ${GATK_Package_jarFile_address} BaseRecalibrator -I ${primaryBam_Tumor}_Sorted_RG_MD -O ${primaryBam_Tumor}_Sorted_RG_MD_bqsr.table -R ${Genome_FASTA} --known-sites ${GATK_Known_site_db1}
java -jar ${GATK_Package_jarFile_address} ApplyBQSR -bqsr ${primaryBam_Tumor}_Sorted_RG_MD_bqsr.table -I ${primaryBam_Tumor}_Sorted_RG_MD -O ${primaryBam_Tumor}_Sorted_RG_MD_BQSR.bam

rm -f ${primaryBam_Normal}_Sorted_RG_MD
rm -f ${primaryBam_Tumor}_Sorted_RG_MD



#Varscan2 somatic variant caller
java -jar ${Varscan2_jarFile_address} somatic <(samtools mpileup --no-BAQ -f ${Genome_FASTA} ${primaryBam_Normal}_Sorted_RG_MD_BQSR.bam ${primaryBam_Tumor}_Sorted_RG_MD_BQSR.bam) ${Varscan2_Output_path} --mpileup 1 --output-vcf

##Should be run by your own
#java -jar VarScan.v2.3.9.jar processSomatic '/media/ad/9117284696_AD/SomaticVariantCalling/Varscan/exome.snp.vcf' 
#java -jar VarScan.v2.3.9.jar processSomatic '/media/ad/9117284696_AD/SomaticVariantCalling/Varscan/exome.indel.vcf' 




#############################################
##Let's get into the broad institute pipeline(GATK4.4)
##Firstly, download af only germline resources using the link below. 
##Keep in mind that this database isn't necessary for the somatic calling but it's recommended.
###  https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37?pli=1
## I also upload the compressed version of the db in my own server for your convinient. Find them on easy databases section at adbioinformatics.net

#Compress af-only database using bgzip
bgzip -@ 12 -c '/media/ad/9117284696_AD/NGSneeds/somatic-b37_af-only-gnomad.raw.sites.vcf' > '/media/ad/9117284696_AD/NGSneeds/somatic-b37_af-only-gnomad.raw.sites.vcf.gz'
#then, index it using tabix tool
tabix -p vcf '/media/ad/9117284696_AD/NGSneeds/somatic-b37_af-only-gnomad.raw.sites.vcf.gz'
#To make af only database compatible with the used reference genome. It needs to add "chr" prefix to chr numbers using command below.
awk -F'\t' '{ if(substr($1, 1, 1) == "#" || substr($1, 1, 3) == "chr") print $0; else print "chr"$0 }' '/media/ad/9117284696_AD/NGSneeds/somatic-b37_af-only-gnomad.raw.sites.vcf' > '/media/ad/9117284696_AD/NGSneeds/somatic-b37_af-only-gnomad.raw.sites_chradded.vcf'

#Creating a panel of normals
#Remember that if you have more than one normal sample, you will need to run this step for each of them one time, then merge them to make a unified PON VCF file
java -jar ${GATK_Package_jarFile_address} Mutect2 --native-pair-hmm-threads 12 -R '/media/ad/9117284696_AD/NGSneeds/Genome_FASTA/hg19/hg19.fa' -I '/media/ad/9117284696_AD/SomaticVariantCalling/Exome_Normal.bam_Sorted_RG_MD_BQSR.bam' -tumor-sample Normal -O '/media/ad/9117284696_AD/SomaticVariantCalling/Mutect2/Exome_Norm_PON.vcf.gz'

#Running Mutect2 Using latest version of GATK
java -jar ${GATK_Package_jarFile_address} Mutect2 --native-pair-hmm-threads 12 -R '/media/ad/9117284696_AD/NGSneeds/Genome_FASTA/hg19/hg19.fa' -I '/media/ad/9117284696_AD/SomaticVariantCalling/Exome_Tumor.bam_Sorted_RG_MD_BQSR.bam' -tumor Tumor -I '/media/ad/9117284696_AD/SomaticVariantCalling/Exome_Normal.bam_Sorted_RG_MD_BQSR.bam' -normal Normal --panel-of-normals '/media/ad/9117284696_AD/SomaticVariantCalling/Mutect2/Exome_Norm_PON.vcf.gz' -O '/media/ad/9117284696_AD/SomaticVariantCalling/Mutect2/finalMutect.vcf.gz' --germline-resource '/media/ad/9117284696_AD/NGSneeds/somatic-b37_af-only-gnomad.raw.sites.vcf.gz' --af-of-alleles-not-in-resource 0.00003 


echo "
Finished!
Presented to you by Milad Eidi
AD Bioinformatics.net
March 2023"
