#Stage_Two

#downloading sample datasets
mkdir -p raw_data 
cd raw_data/

wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

#downloading reference sequence 
mkdir reference
cd reference/
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#unzip reference
gunzip hg19.chr5_12_17.fa.gz

cd ../

#install fastqc and multiqc on the terminal
conda install -c bioconda fastqc multiqc --yes

#create directory for the fastqc output
mkdir -p Fastqc_Reports 

#fastqc analysis 
fastqc raw_data/*.fastq.gz -o Fastqc_Reports

#multiqc analysis
multiqc Fastqc_Reports -o Fastqc_Reports

#downloading fastp
conda install -c bioconda fastp

cd raw_data/
nano trimmed.sh
bash trimmed.sh

#applying multiqc analysis on fastp results
multiqc trimmed_reads -o trimmed_reads 
cd ../

#Read mapping
cd reference/

#Index reference file	
bwa index hg19.chr5_12_17.fa 
cd ../

#perform alignment
mkdir Mapping

bwa mem -t 8 reference/hg19.chr5_12_17.fa \
raw_data/trimmed_reads/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz raw_data/trimmed_reads/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz \
| samtools view -b > Mapping/SLGFSK-N.sam

bwa mem -t 8 reference/hg19.chr5_12_17.fa \
raw_data/trimmed_reads/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz raw_data/trimmed_reads/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz \
| samtools view -b > Mapping/SLGFSK-T.sam

#convert SAM to BAM 
samtools view -S -b SLGFSK-N.sam > SLGFSK-N.bam
samtools view -S -b SLGFSK-T.sam > SLGFSK-T.bam

#sort BAM file
samtools sort SLGFSK-N.bam -o SLGFSK-N.sorted.bam
samtools sort SLGFSK-T.bam -o SLGFSK-T.sorted.bam

#index sorted BAM file
samtools index SLGFSK-N.sorted.bam
samtools index SLGFSK-T.sorted.bam

#mapped reads filtering 
samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/SLGFSK-N.sorted.bam > Mapping/SLGFSK-N.filtered1.bam 
samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/SLGFSK-T.sorted.bam > Mapping/SLGFSK-T.filtered1.bam

#duplicates removal
samtools rmdup SLGFSK-N.sorted.bam  SLGFSK-N.rdup 
samtools rmdup SLGFSK-T.sorted.bam  SLGFSK-T.rdup

#leftalign BAM 
cat Mapping/SLGFSK-N.rdup  | bamleftalign -f reference/hg19.chr5_12_17.fa -m 5 -c > Mapping/SLGFSK-N.leftAlign.bam
cat Mapping/SLGFSK-T.rdup  | bamleftalign -f reference/hg19.chr5_12_17.fa -m 5 -c > Mapping/SLGFSK-T.leftAlign.bam

#recalibrate read mapping qualities
samtools calmd -@ 32 -b Mapping/SLGFSK-N.leftAlign.bam reference/hg19.chr5_12_17.fa > Mapping/SLGFSK-N.recalibrate.bam
samtools calmd -@ 32 -b Mapping/SLGFSK-T.leftAlign.bam reference/hg19.chr5_12_17.fa > Mapping/SLGFSK-T.recalibrate.bam

#Refilter read mapping qualities
bamtools filter -in Mapping/SLGFSK-N.recalibrate.bam -mapQuality "<=254 >" Mapping/SLGFSK-N.refilter.bam
bamtools filter -in Mapping/SLGFSK-T.recalibrate.bam -mapQuality "<=254 >" Mapping/SLGFSK-T.refilter.bam

#Variant calling and classification

wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar	

mkdir Variants 

#convert data to pileup 
samtools mpileup -f reference/hg19.chr5_12_17.fa Mapping/SLGFSK-N.refilter.bam --min-MQ 1 --min-BQ 28 > Variants/SLGFSK-N.pileup
samtools mpileup -f reference/hg19.chr5_12_17.fa Mapping/SLGFSK-T.refilter.bam --min-MQ 1 --min-BQ 28 > Variants/SLGFSK-T.pileup

#call Variants
java -jar VarScan.v2.3.9.jar somatic Variants/SLGFSK-N_231335.pileup \
Variants/SLGFSK-T_231336.pileup Variants/SLGFSK --normal-purity 1  --tumor-purity 0.5 --output-vcf 1 

cd Variants/ 

#change the names of the ouput
mv SLGFSK--normal-purity.snp.vcf SLGFSK.snp.vcf 

#merge vcf
bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
tabix Variants/SLGFSK.snp.vcf.gz
tabix Variants/SLGFSK.indel.vcf.gz
bcftools merge --force-sample Variants/SLGFSK.snp.vcf.gz Variants/SLGFSK.indel.vcf.gz > Variants/SLGFSK.vcf

#Variant Annotation
#Functional Annotation using SnpEff

#download snpEff database
snpEff download hg19

#annotate variants
snpEff hg19 Variants/SLGFSK.vcf > Variants/SLGFSK.ann.vcf
