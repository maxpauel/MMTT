		# STEP 1 - QUALITY CONTROL AND TRIMMING
# Tools: 
#	- FastQC
#	- Trimmomatic
# Input files:
#	- ./Fastq/*.fastq - raw data files
#	- ./adapters.fa - NEB adapter sequences
# Output files:
#	- ./Fastq_trimmed/*.trim.fastq.gz - trimmed fastq files
#	- ./Fastq/*.fastqc.html - quality control files
mkdir ./Fastq_trimmed
path_in='./Fastq'
path_out='./Fastq_trimmed'
for i in $path_in/*.fastq.gz; do {
res_file=${i%.fastq*}.trim.fastq.gz
res_file=${res_file/$path_in/}
res_file=$path_out$res_file
# quality control of raw data
fastqc $i -o ./Fastq;
# adapter trimming and deletion of low quality reads
java -jar /usr/share/java/trimmomatic.jar SE -threads 8 $i $res_file ILLUMINACLIP:./adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35;
} ; done


		# STEP 2 - ALIGNMENT
# Tools:
#	- Hisat2
#	- Samtools
# Input files:
# 	- ./Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz - genome fasta file, source - https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# 	- ./Homo_sapiens.GRCh38.101.gtf - genome annotation file, source -  https://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
#	  - ./Fastq_trimmed/*.trim.fastq.gz - trimmed fastq files
# Output files:
#	- ./bam/*.sorted.bam - sorted bam files

# build genome index
hisat2-build ./Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ./hisat2_index/Homo_sapiens.GRCh38.dna.primary_assembly
# generate splice sites txt
hisat2_extract_splice_sites.py ./Homo_sapiens.GRCh38.101.gtf > ./GRCh38.101.splicesites.txt
# align trimmed fastq
ht_index='./hisat2_index/Homo_sapiens.GRCh38.dna.primary_assembly'
path_in='./Fastq_trimmed/'
path_out='./bam/'
for i in $path_in/*.trim.fastq.gz; do {
res_file=${i%.trim.fastq*}.sam
res_file=${res_file/$path_in/}
res_file=$path_out$res_file
hisat2 -p 8 -x $ht_index -U $i  -S $res_file \
--known-splicesite-infile ./GRCh38.101.splicesites.txt \
--rna-strandness R;
samtools view -@ 8 -bSo ${res_file%.sam*}.bam $res_file;
rm $res_file;
samtools sort -@ 8 ${res_file%.sam*}.bam -o ${res_file%.sam*}.sorted.bam;
rm ${res_file%.sam*}.bam;
} ; done

		# STEP 3 - CALCULATE READ COUNTS 
# Tools:
#	- R, R packages (Rsubread)
# Input files:
#	  - ./bam/*.sorted.bam - sorted bam files
# 	- ./Homo_sapiens.GRCh38.101.gtf - genome annotation file, source -  https://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
# Output files:
#	- ./bam/counts.txt - read count file

ls ./bam/*.sorted.bam >./bam/bam_list.txt
R
library(Rsubread)
files=readLines('./bam/bam_list.txt')
a=featureCounts(files,
isGTFAnnotationFile = TRUE,
annot.ext ='./Homo_sapiens.GRCh38.101.gtf',
GTF.featureType="exon",
GTF.attrType="gene_id",
strandSpecific=2,
countMultiMappingReads=FALSE,
allowMultiOverlap=F,
nthreads=8)
write.table(a[[1]],'./bam/counts.txt',quote=F,sep='\t')
q()
n

		# STEP 4 - DIFFERENTIAL GENE EXPRESSION
# Tools: 
#	- R, R packages (DESeq2)
# Input files:

# Output files:


R
library(DESeq2)
library(rtracklayer)
cts=read.table('./bam/counts.txt',header=T,row.names=1)
design=read.csv('./design.txt',row.names=1)
cts_basal=cts1[,c(TRUE,FALSE)]
design_basal=design1[c(TRUE,FALSE),]
dds <- DESeqDataSetFromMatrix(countData = cts_g,
                              colData = coldata,
                              design = ~ con_group1)
dds1 <- DESeq(dds)
res1_basal=results(dds1,contrast=c('con_group1','Ob','H'))
res2_basal=results(dds1,contrast=c('con_group1','T2d','H'))
res3_basal=results(dds1,contrast=c('con_group1','T2d','Ob'))
res_basal=cbind(res1[,c(2,6)],res2[,c(2,6)],res3[,c(2,6)])
colnames(res_basal)=c('LFC_ObH','Padj_ObH','LFC_T2dH','Padj_T2dH','LFC_T2dOb','Padj_T2dOb')

cts_t2d=cts[,1:12]
cts_h=cts[,13:28]
cts_ob=cts[,29:42]
design_t2d=design[1:12,]
design_h=design[13:28,]
design_ob=design[29:42,]
coldata=design_t2d
dds <- DESeqDataSetFromMatrix(countData = cts_t2d,
                              colData = coldata,
                              design = ~ person+group)
dds1 <- DESeq(dds)
res1=results(dds1,contrast=c('group','Post','Pre'))
coldata=design_h
dds <- DESeqDataSetFromMatrix(countData = cts_h,
                              colData = coldata,
                              design = ~ person+group)
dds2 <- DESeq(dds)
res2=results(dds2,contrast=c('group','Post','Pre'))
coldata=design_ob
dds <- DESeqDataSetFromMatrix(countData = cts_ob,
                              colData = coldata,
                              design = ~ person+group)
dds3 <- DESeq(dds)
res3=results(dds3,contrast=c('group','Post','Pre'))
resMMT=cbind(res1[,c(2,6)],res2[,c(2,6)],res3[,c(2,6)])
colnames(resMMT)=c('LFC_T2d','Padj_T2d','LFC_H','Padj_H','LFC_Ob','Padj_Ob')
results=cbind(res_basal,resMMT)

gtf=readGFF('./Homo_sapiens.GRCh38.101.gtf')
gtf=gtf[gtf$type=='gene',c(9,11,13)]
DEG=merge(gtf,results,by.x='gene_id',by.y='row.names')
write.table(DEG,'./bam/DEG_analysis.txt',quote=F,sep='\t',row.names=F)
q()
n

		# STEP 5 - TPM CALCULATION
path_in='./Fastq_trimmed/'
path_out='./TPM'
for i in $path_in/*.trim.fastq.gz; do {
mkdir ${i%.trim.fastq*}
kallisto quant \
--rf-stranded \
-b 0 \
-i /media/maxpauel/01B5FA3D077CA3E9/Dry_immersion/Kallisto/transcripts.idx \
-o ${i%.trim.fastq*} \
--single  \
-l 200  \
-s 20  \
-t 8 \
$i;
} ; done

R
txt=readLines('/media/maxpauel/01B5FA3D077CA3E9/Diabetus/Fastq_trim/files.txt')
datalist = lapply(txt, function(x)read.table(x, header=T)) 
datafr = do.call("cbind", datalist)
rownames(datafr)=datafr$target_id
data=datafr[,c(FALSE,FALSE,FALSE,FALSE,TRUE)]
write.table(data,'/media/maxpauel/01B5FA3D077CA3E9/Diabetus/TPM/TPM.txt',quote=F,sep='\t',row.names=T)



