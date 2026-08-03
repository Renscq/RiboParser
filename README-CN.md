<!--
 * @Author: 'rensc' 'rensc0718@163.com'
 * @Date: 2024-10-15 11:44:58
 * @LastEditors: 'rensc' 'rensc0718@163.com'
 * @LastEditTime: 2026-06-10 17:04
 * @FilePath: \RiboParser\README-CN.md 
 * 
-->

# RiboParser

```
Ren, S., Li, Y. & Zhou, Z. 
RiboParser/RiboShiny: An integrated platform for comprehensive analysis and visualization of ribo-seq data. 
Journal of Genetics and Genomics (2025) 
doi:10.1016/j.jgg.2025.04.010.
```

为了简化理解和应用，我们将分析可公开访问的项目数据，分解每个分析步骤以说明完整的工作流程。

这个过程包括一般的分析步骤和专门的分析和可视化技术，由`RiboParser`和`RiboShiny`提供便利。

具体步骤如下：

1. 软件安装
2. 创建参考文件
3. 原始数据下载
4. 原始数据清理
5. 数据比对
6. 测序质量分析
7. 基因水平分析
8. 密码子水平分析

该数据分析的结果可以在`RiboShiny`中进一步分析和可视化。

## 1. 软件配置

首先，需要安装 `miniconda` 或者 `micromamba` 来管理软件环境。

### 1.1 用`conda`创建环境

```bash
conda create -n ribo
conda activate ribo
```

### 1.2 使用 `conda` 安装 `riboparser` 和软件依赖项

```bash
# conda
conda install riboparser -c rensc -c conda-forge -c bioconda

# mamba
micromamba install riboparser -c rensc -c conda-forge -c bioconda

```

### 1.3 按步骤安装 `riboparser` 和软件依赖包

利用单个命令安装依赖包。

```bash
conda install cutadapt bowtie samtools star bedtools subread rsem gffread sra-tools \
 ucsc-genepredtogtf ucsc-gtftogenepred ucsc-gff3togenepred ucsc-bedgraphtobigwig ucsc-bedsort \
 -c bioconda

conda install pigz -c conda-forge
```


利用`pip`安装`RiboParser`

```bash
pip install riboparser

```

或者，我们可以从GitHub下载版本，重新 build，然后安装。

```bash
cd RiboParser

# 安装构建工具
pip install build

# 构建包
python -m build

# 安装本地包
pip install .
```


### 1.4 运行测试
测试软件的依赖性、安装和操作问题。

```bash
# 检查版本信息
riboparser -v
# 检查参考文献信息
riboparser -c
# 检查依赖环境信息
riboparser -d
# 检查包中模块信息
riboparser -m

```

## 2. 准备参考文件

### 2.1 完整的项目目录示例如下所示

完整的数据分析包括参考文件准备、原始数据、RNA-seq数据分析和Ribo-seq数据分析。

```bash
$ cd && cd ./sce/
$ tree -d

.
├── 1.reference
│   ├── genome
│   ├── mrna
│   ├── ncrna
│   ├── norm
│   ├── rrna
│   ├── rsem-index
│   ├── star-index
│   └── trna
├── 2.rawdata
│   ├── ribo-seq
│   └── rna-seq
├── 3.rna-seq
│   ├── 1.cleandata
│   ├── 2.bowtie
│   ├── 3.star
│   ├── 4.quantification
│   └── 5.riboparser
│       ├── 01.qc
│       ├── 02.digestion
│       ├── 03.offset
│       ├── 04.density
│       ├── 05.merge
│       ├── 06.periodicity
│       ├── 07.metaplot
│       ├── 08.coverage
│       ├── 09.correlation
│       ├── 10.shuffle
│       └── 11.retrieve
└── 4.ribo-seq
    ├── 1.cleandata
    ├── 2.bowtie
    ├── 3.star
    ├── 4.quantification
    └── 5.riboparser
        ├── 01.qc
        ├── 02.digestion
        ├── 03.offset
        ├── 04.density
        ├── 05.merge
        ├── 06.periodicity
        ├── 07.metaplot
        ├── 08.coverage
        ├── 09.correlation
        ├── 10.quantification
        ├── 11.pausing_score
        ├── 12.codon_occupancy
        ├── 13.codon_decoding_time
        ├── 14.codon_selection_time
        ├── 15.coefficient_of_variation
        ├── 16.meta_codon
        ├── 17.shuffle
        ├── 18.retrieve
        └── 19.frame_shift
```

### 2.2 准备参考基因组索引
#### 2.2.1 创建目录

创建文件夹来保存不同类型的引用序列文件。

```bash
$ mkdir -p ./sce/1.reference/
$ cd ./sce/1.reference/
$ mkdir cdna genome gtf mrna ncrna rrna trna norm rsem-index
```

#### 2.2.2 从NCBI下载参考文件

使用最常用的数据分析文件格式，基因组序列采用fasta格式，参考文件采用GTF或GFF3格式。

```bash
# genome sequence
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

# GTF or GFF3
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz

# cDNA sequence
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz

# feature table
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_feature_table.txt.gz

# decompression
$ gunzip *.gz

$ gffread -g GCF_000146045.2_R64_genomic.fna GCF_000146045.2_R64_genomic.gff -F -w cdna.fa
```

#### 2.2.3 使用`bowtie`创建`genome`索引

```bash
$ bowtie-build ../GCF_000146045.2_R64_genomic.fna ./genome/genome --threads 12 &>> ./genome/genome_build.log
```

#### 2.2.4 使用`bowtie`创建`mRNA`索引

这里使用自定义脚本提取相应的序列信息,根据序列名称从`fasta`文件中查找。

```bash
$ retrieve_seq -h

usage: retrieve_seq [-h] [-v] -i INPUT -n NAME [-u UNMAPPED] -o OUTPUT

This script is used to retrieve the fasta sequence by name.

options:
  -h, --help     show this help message and exit
  -v, --version  show programs version number and exit
  -i INPUT       input the fasta file
  -n NAME        gene ids in txt format
  -u UNMAPPED    output the unmapped gene ids
  -o OUTPUT      prefix of output file name (default results_peaks.txt)
```


```bash
# filter the mrna sequence
$ grep -i 'gbkey=mRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./mrna/mrna.ids
$ retrieve_seq -i ./cdna.fa -n ./mrna/mrna.ids -o ./mrna/mrna.fa &>> ./mrna/mrna_build.log
# build the mrna index
$ bowtie-build ./mrna/mrna.fa ./mrna/mrna --threads 12 &>> ./mrna/mrna_build.log
```

#### 2.2.5 使用`bowtie`创建`rRNA`索引
```bash
# filter the rrna sequence
$ grep -i 'gbkey=rRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./rrna/rrna.ids
$ retrieve_seq -i ./cdna.fa -n ./rrna/rrna.ids -o ./rrna/rrna.fa &>> ./rrna/rrna_build.log
# build the rrna index
$ bowtie-build ./rrna/rrna.fa ./rrna/rrna --threads 12 &>> ./rrna/rrna_build.log
```

#### 2.2.6 使用`bowtie`创建`tRNA`索引
```bash
# filter the trna sequence
$ grep -i 'gbkey=tRNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./trna/trna.ids
$ retrieve_seq -i ./cdna.fa -n ./trna/trna.ids -o ./trna/trna.fa &>> ./trna/trna_build.log
# build the trna index
$ bowtie-build ./trna/trna.fa ./trna/trna --threads 12 &>> ./trna/trna_build.log
```


#### 2.2.7 使用“`bowtie`创建`ncRNA`索引
```bash
# filter the ncrna sequence
$ grep -iE 'gbkey=ncRNA|gbkey=lnc_RNA|gbkey=miRNA|gbkey=snoRNA|gbkey=snRNA|gbkey=misc_RNA' ./cdna.fa | cut -d ' ' -f 1 | cut -c 2- > ./ncrna/ncrna.ids
$ retrieve_seq -i ./cdna.fa -n ./ncrna/ncrna.ids -o ./ncrna/ncrna.fa &>> ./ncrna/ncrna_build.log
# build the ncrna index
$ bowtie-build ./ncrna/ncrna.fa ./ncrna/ncrna --threads 12 &>> ./ncrna/ncrna_build.log
```

#### 2.2.8 标准化的`gtf`或`gff3`文件

- `rpf_Reference`的解释

```bash
$ rpf_Reference -h

usage: rpf_Reference [-h] -g GENOME -t GTF -o OUTPUT [-u UTR] [-c] [-l] [-w]

This script is used to build the references for the RiboParser.

options:
  -h, --help  show this help message and exit
  -u UTR      add the pseudo UTR to the leaderless transcripts (default: 0 nt).
  -c          only retain the protein coding transcripts (default: False).
  -l          only retain the longest protein coding transcripts, it is recommended to select 
              the longest transcript for subsequent analysis. (default: False).
  -w          output whole message (default: False).

Required arguments:
  -g GENOME   the input file name of genome sequence
  -t GTF      the input file name of gtf file
  -o OUTPUT   the prefix of output file. (prefix + _norm.gtf)
```

- 利用GTF文件和基因组fasta文件创建参考文件

```bash
$ rpf_Reference \
 -g ../GCF_000146045.2_R64_genomic.fna \
 -t ../GCF_000146045.2_R64_genomic.gff \
 -u 30 -o ./norm/gene &>> ./norm/norm_build.log
```

#### 2.2.9 使用`STAR`创建`genome`索引

```bash
$ STAR \
 --genomeSAindexNbases 11 \
 --runThreadN 12 \
 --runMode genomeGenerate \
 --genomeDir ./star-index \
 --genomeFastaFiles GCF_000146045.2_R64_genomic.fna \
 --sjdbGTFfile ./norm/gene.norm.gtf
```

#### 2.2.10 使用`rsem`创建`transcriptome`索引

```bash
$ rsem-prepare-reference \
 -p 12 \
 --gtf ../norm/gene.norm.gtf ../GCF_000146045.2_R64_genomic.fna ./rsem-index/rsem
```


## 3. 数据预处理和比对

为了介绍`RiboParser`的分析过程和使用方法，本文以`GSE67387`数据集的RNA-seq和Ribo-seq数据为例。

```shell
# dataset
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67387

# reference
Nedialkova DD, Leidel SA. Optimization of Codon Translation Rates via tRNA Modifications Maintains Proteome Integrity. Cell 2015 Jun 18;161(7):1606-18. 
PMID: 26052047
```


### 3.1 `GSE67387`数据集RNA-seq数据基础分析

#### 3.1.1 下载RNA-seq原始数据

使用`sra-tools`中的`prefetch`下载原始sra格式数据并将其提取为`fastq`文件。

```bash
$ mkdir -p ./sce/2.rawdata/rna-seq/
$ cd ./sce/2.rawdata/rna-seq/

#################################################
# download rna-seq
$ prefetch -o SRR1944925.sra SRR1944925
$ prefetch -o SRR1944926.sra SRR1944926
$ prefetch -o SRR1944927.sra SRR1944927
$ prefetch -o SRR1944928.sra SRR1944928
$ prefetch -o SRR1944929.sra SRR1944929
$ prefetch -o SRR1944930.sra SRR1944930
$ prefetch -o SRR1944931.sra SRR1944931
$ prefetch -o SRR1944932.sra SRR1944932
$ prefetch -o SRR1944933.sra SRR1944933
$ prefetch -o SRR1944934.sra SRR1944934
$ prefetch -o SRR1944935.sra SRR1944935

# decompression
for sra in *.sra
do

fastq-dump $sra
pigz *fastq

done
```

#### 3.1.2 RNA-seq数据清洗

因为来自gse项目的数据是被清理的，所以它不包括适配器和索引序列。所以下面只是展示一般步骤，不需要运行。

1. RNA-seq数据清洗

```bash
$ cd
$ mkdir -p ./sce/3.rna-seq/1.cleandata/
$ cd ./sce/3.rna-seq/1.cleandata/

#################################################
# run the cutadapt
for fq in ../../2.rawdata/rna-seq/*fastq.gz
do
cutadapt --match-read-wildcards \
 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC \
 -m 25 -O 6 -j 12 \
 -o `\basename $fq fastq.gz`clean.fastq.gz $fq &> $fq".log"

done
```

#### 3.1.3 将干净的数据与不同类型的参考文件进行比对

为了评估文库质量并消除来自不同非编码RNA(ncRNAs)的reads对后续分析的影响，我们采用`bowtie`对测序数据中的reads进行分类。

在正常情况下，特别是使用oligo（dT）方法构建的RNA-seq文库，大多数reads来自mRNA。因此，对于RNA-seq分析，一般不需要这一步。适用于rRNA耗竭法构建的文库。

1. 将RNA-seq数据与参考文献进行比对

```bash
$ cd
$ mkdir -p ./sce/3.rna-seq/2.bowtie/
$ cd ./sce/3.rna-seq/2.bowtie/

#################################################
# set database
rrna='../../sce/1.reference/rrna/rrna'
trna='../../sce/1.reference/trna/trna'
ncrna='../../sce/1.reference/ncrna/ncrna'
mrna='../../sce/1.reference/mrna/mrna'
chrom='../../sce/1.reference/genome/genome'

threads=12
mismatch=1

# alignment reads to reference
for fq in ../1.cleandata/*fastq.gz
do
fqname=`\basename $fq .fastq.gz`

## rrna
bowtie -p $threads -v $mismatch --un="$fqname".norrna.fq --al="$fqname".rrna.fq \
 -x $rrna $fq -S "$fqname".rrna.sam 2>> "$fqname".log

## trna
bowtie -p $threads -v $mismatch --un="$fqname".notrna.fq --al="$fqname".trna.fq \
 -x $trna "$fqname".norrna.fq -S "$fqname".trna.sam 2>> "$fqname".log

## ncrna
bowtie -p $threads -v $mismatch --un="$fqname".noncrna.fq --al="$fqname".ncrna.fq \
 -x $ncrna "$fqname".notrna.fq -S "$fqname".ncrna.sam 2>> "$fqname".log

## mrna
bowtie -p $threads -v $mismatch --un="$fqname".nomrna.fq --al="$fqname".mrna.fq \
 -x $mrna "$fqname".noncrna.fq -S "$fqname".mrna.sam 2>> "$fqname".log

## genome
bowtie -p $threads -v $mismatch --un="$fqname".nogenome.fq --al="$fqname".genome.fq \
 -x $chrom "$fqname".nomrna.fq -S "$fqname".genome.sam 2>> "$fqname".log

## compress fastq
pigz *fq

## compress sam
for sam in *.sam
do

samtools view -h -F 4 $sam | samtools sort -@ $threads -o `\basename $sam sam`bam
rm $sam

done

done
```


2. 所有参考文献的统计比对结果。

-  `merge_bwt_log`的解释

```bash
$ merge_bwt_log -h

Step1: Checking the input Arguments.
usage: merge_bwt_log [-h] -l LIST [LIST ...] -o OUTPUT [-n NAME]

This script is used to statstic the mapped reads from log files

options:
  -h, --help            show this help message and exit
  -n NAME, --name NAME  set the name of each database (default: rRNA,tRNA,ncRNA,mRNA,Genome).

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for bowtie mapping log files (e.g., '*log').
  -o OUTPUT             prefix of output file name.
```

- 统计比对结果

```bash
#################################################
# merge all log files
merge_bwt_log -n rRNA,tRNA,ncRNA,mRNA,Genome -l *log -o RNA_seq &>> merge_bowtie.log
```

#### 3.1.4 使用`STAR`对mRNA reads进行比对

在去除ncRNA读段后，使用`STAR`将剩余的干净reads重新与酵母基因组比对。

1. 使用`STAR`对mRNA reads （RNA-seq）进行比对
```bash
$ cd
$ mkdir -p ./sce/3.rna-seq/3.star/
$ cd ./sce/3.rna-seq/3.star/

#################################################
# set the option and database
genome='../../1.reference/star-index/'
threads=12

#################################################
# map the all rna-seq reads to genome and transcriptome region
for fastq in ../2.bowtie/*.noncrna.fq.gz
do

## get file name
output=$(basename $fastq .noncrna.fq.gz)

#################################################
## run the alignment
STAR --runThreadN $threads \
 --readFilesCommand zcat \
 --genomeDir $genome \
 --readFilesIn $fastq \
 --outFileNamePrefix $output \
 --outSAMtype BAM Unsorted \
 --outFilterType BySJout \
 --quantMode TranscriptomeSAM GeneCounts \
 --outReadsUnmapped Fastx \
 --outSAMattributes All \
 --alignEndsType Local \
 --outFilterMultimapNmax 3 \
 --outFilterMismatchNmax 1 \
 --alignIntronMax 10000 \
 --outFilterMatchNmin 20
# --outWigType wiggle --outWigNorm RPM

pigz *mate1

#################################################
## sort the bam file
samtools sort -@ $threads $output"Aligned.out.bam" -o $output"Aligned.sortedByCoord.out.bam"
samtools index -@ $threads $output"Aligned.sortedByCoord.out.bam"
rm $output"Aligned.out.bam"

done
```

#### 3.1.5 用`RSEM`或`featureCounts`量化基因表达水平

`RSEM`和`featureCounts`都可以用来量化基因表达水平。为了分析的目的，我们将使用`RSEM`作为代表性工具。

1. 利用RNA-seq数据量化转录物丰度

```bash
$ cd
$ mkdir -p ./sce/3.rna-seq/4.quantification/
$ cd ./sce/3.rna-seq/4.quantification/

#################################################
# quantify the gene expression
for bam in ../3.star/*Aligned.toTranscriptome.out.bam
do
rsem-calculate-expression -p 10 \
 --no-bam-output --alignments \
 -q $bam ../../1.reference/rsem-index/rsem `\basename $bam Aligned.toTranscriptome.out.bam`
# rsem-calculate-expression -p 10 \
# --paired-end --no-bam-output --alignments \
# -q $bam ../../1.reference/rsem-index/rsem `\basename $bam Aligned.toTranscriptome.out.bam`

done
```

2. 整合所有样品的RNA-seq定量值

- `merge_rsem`的解释

```bash
$ merge_rsem -h

usage: merge_rsem [-h] -l LIST [LIST ...] -o OUTPUT [-c {expected_count,TPM,FPKM}]

This script is used to merge specified columns from result files

options:
  -h, --help            show this help message and exit
  -c {expected_count,TPM,FPKM}, --column {expected_count,TPM,FPKM}
                        Column name to merge (e.g., 'expected_count').

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for result files (e.g., '*results').
  -o OUTPUT             output file name.
```

- 整合RNA-seq定量值

```bash
#################################################
# merge the gene expression
merge_rsem -c expected_count -l *.genes.results -o gene.expected_count.txt &>> merge_rsem.log
merge_rsem -c TPM -l *.genes.results -o gene.TPM.txt &>> merge_rsem.log
merge_rsem -c FPKM -l *.genes.results -o gene.FPKM.txt &>> merge_rsem.log

#################################################
# merge the isoforms expression
merge_rsem -c expected_count -l *.isoforms.results -o isoforms.expected_count.txt &>> merge_rsem.log
merge_rsem -c TPM -l *.isoforms.results -o isoforms.TPM.txt &>> merge_rsem.log
merge_rsem -c FPKM -l *.isoforms.results -o isoforms.FPKM.txt &>> merge_rsem.log
```


### 3.2 `GSE67387`数据集Ribo-seq数据的基本分析

#### 3.2.1 下载Ribo-seq原始数据

使用`sra-tools`中的`prefetch`下载原始sra格式数据并将其提取为`fastq`格式文件。

```bash
$ cd
$ mkdir -p ./sce/2.rawdata/ribo-seq/
$ cd ./sce/2.rawdata/ribo-seq/

#################################################
# download ribo-seq
prefetch -o SRR1944912.sra SRR1944912
prefetch -o SRR1944913.sra SRR1944913
prefetch -o SRR1944914.sra SRR1944914
prefetch -o SRR1944915.sra SRR1944915
prefetch -o SRR1944916.sra SRR1944916
prefetch -o SRR1944917.sra SRR1944917
prefetch -o SRR1944918.sra SRR1944918
prefetch -o SRR1944919.sra SRR1944919
prefetch -o SRR1944920.sra SRR1944920
prefetch -o SRR1944921.sra SRR1944921
prefetch -o SRR1944922.sra SRR1944922
prefetch -o SRR1944923.sra SRR1944923

# decompression
for sra in *.sra
do

fastq-dump $sra
pigz *fastq

done
```

#### 3.2.2 Ribo-seq数据清洗

```bash
$ cd
$ mkdir -p ./sce/4.ribo-seq/1.cleandata/
$ cd ./sce/4.ribo-seq/1.cleandata/

#################################################
# run the cutadapt
for fq in ../../2.rawdata/ribo-seq/*fastq.gz
do
cutadapt --match-read-wildcards \
 -a AAAAAAAA \
 -m 25 -O 6 -j 10 \
 -o `\basename $fq fastq.gz`clean.fastq.gz $fq &> $fq".log"
done
```


#### 3.2.3 将干净的数据与不同类型的参考文件进行比对

为了评估文库质量并消除来自不同非编码RNA(ncRNAs)的reads对后续分析的影响，我们采用`bowtie`对测序数据中的reads进行分类。

1. Ribo-seq数据比对
```bash
$ cd
$ mkdir -p ./sce/4.ribo-seq/2.bowtie/
$ cd ./sce/4.ribo-seq/2.bowtie/

#################################################
# set database
rrna='../../sce/1.reference/rrna/rrna'
trna='../../sce/1.reference/trna/trna'
ncrna='../../sce/1.reference/ncrna/ncrna'
mrna='../../sce/1.reference/mrna/mrna'
chrom='../../sce/1.reference/genome/genome'

threads=12
mismatch=1

# alignment reads to reference
for fq in ../1.cleandata/*fastq.gz
do
fqname=`\basename $fq .fastq.gz`

## rrna
bowtie -p $threads -v $mismatch --un="$fqname".norrna.fq --al="$fqname".rrna.fq \
 -x $rrna $fq -S "$fqname".rrna.sam 2>> "$fqname".log

## trna
bowtie -p $threads -v $mismatch --un="$fqname".notrna.fq --al="$fqname".trna.fq \
 -x $trna "$fqname".norrna.fq -S "$fqname".trna.sam 2>> "$fqname".log

## ncrna
bowtie -p $threads -v $mismatch --un="$fqname".noncrna.fq --al="$fqname".ncrna.fq \
 -x $ncrna "$fqname".notrna.fq -S "$fqname".ncrna.sam 2>> "$fqname".log

## mrna
bowtie -p $threads -v $mismatch --un="$fqname".nomrna.fq --al="$fqname".mrna.fq \
 -x $mrna "$fqname".noncrna.fq -S "$fqname".mrna.sam 2>> "$fqname".log

## genome
bowtie -p $threads -v $mismatch --un="$fqname".nogenome.fq --al="$fqname".genome.fq \
 -x $chrom "$fqname".nomrna.fq -S "$fqname".genome.sam 2>> "$fqname".log

## compress fastq
pigz *fq

## compress sam
for sam in *.sam
do

samtools view -h -F 4 $sam | samtools sort -@ $threads -o `\basename $sam sam`bam
rm $sam

done

done
```

2. 统计所有数据库的比对结果。
```bash
#################################################
# merge all log files
merge_bwt_log -n rRNA,tRNA,ncRNA,mRNA,Genome -l *log -o sce &>> merge_bowtie.log
```


#### 3.2.4 使用`STAR`对mRNA reads (Ribo-seq)进行比对

在去除ncRNA读段后，使用`STAR`将剩余的干净reads重新与酵母基因组进行比对。

```bash
$ cd
$ mkdir -p ./sce/4.ribo-seq/3.star/
$ cd ./sce/4.ribo-seq/3.star/

#################################################
# set the option and database
genome='../../1.reference/star-index/'
threads=12

#################################################
# map the all rna-seq reads to genome and transcriptome region
for fastq in ../2.bowtie/*.noncrna.fq.gz
do

## get file name
output=$(basename $fastq .noncrna.fq.gz)

#################################################
## run the alignment
STAR --runThreadN $threads \
 --readFilesCommand zcat \
 --genomeDir $genome \
 --readFilesIn $fastq \
 --outFileNamePrefix $output \
 --outSAMtype BAM Unsorted \
 --outFilterType BySJout \
 --quantMode TranscriptomeSAM GeneCounts \
 --outReadsUnmapped Fastx \
 --outSAMattributes All \
 --alignEndsType Local \
 --outFilterMultimapNmax 3 \
 --outFilterMismatchNmax 1 \
 --alignIntronMax 10000 \
 --outFilterMatchNmin 20
# --outWigType wiggle --outWigNorm RPM

pigz *mate1

#################################################
## sort the bam file
samtools sort -@ $threads $output"Aligned.out.bam" -o $output"Aligned.sortedByCoord.out.bam"
samtools index -@ $threads $output"Aligned.sortedByCoord.out.bam"
rm $output"Aligned.out.bam"

done
```


#### 3.2.5 用`RSEM`或`featureCounts`量化基因表达水平

`RSEM`和`featureCounts`都可以用来量化基因表达水平。为了分析的目的，我们将使用`RSEM`作为代表性工具。

1. 利用Ribo-seq数据量化转录物丰度

```bash
$ cd
$ mkdir -p ./sce/4.ribo-seq/4.quantification/
$ cd ./sce/4.ribo-seq/4.quantification/

#################################################
# quantify the isoforms expression
for bam in ../3.star/*Aligned.toTranscriptome.out.bam
do

rsem-calculate-expression -p 12 \
 --no-bam-output --alignments \
 -q $bam ../../1.reference/rsem-index/rsem `\basename $bam Aligned.toTranscriptome.out.bam`

done
```

2. 整合所有样品的Ribo-seq定量值

```bash
#################################################
# merge the gene expression
merge_rsem -c expected_count -l *.genes.results -o gene.expected_count.txt &>> merge_rsem.log
merge_rsem -c TPM -l *.genes.results -o gene.TPM.txt &>> merge_rsem.log
merge_rsem -c FPKM -l *.genes.results -o gene.FPKM.txt &>> merge_rsem.log

#################################################
# merge the isoforms expression
merge_rsem -c expected_count -l *.isoforms.results -o isoforms.expected_count.txt &>> merge_rsem.log
merge_rsem -c TPM -l *.isoforms.results -o isoforms.TPM.txt &>> merge_rsem.log
merge_rsem -c FPKM -l *.isoforms.results -o isoforms.FPKM.txt &>> merge_rsem.log
```


## 4. 使用`RiboParser`对`GSE67387`进行RNA-seq数据分析

### 4.0 准备用于存储结果的目录

```bash
$ cd
$ mkdir -p ./3.rna-seq/5.riboparser/
$ cd ./3.rna-seq/5.riboparser/
$ mkdir 01.qc 02.digestion 03.offset 04.density 05.merge \
 06.periodicity 07.metaplot 08.coverage 09.correlation 10.shuffle
```

### 4.1 测序数据的质量检查

1. RNA-seq数据的质量检查

```bash
$ cd ./01.qc

#################################################
# check the ribo-seq quality
for bam in ../../3.star/*Aligned.toTranscriptome.out.bam
do
prefix_name=$(basename $bam Aligned.toTranscriptome.out.bam)

rpf_Check -b $bam -s --thread 10 \
 -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done
```

2. 整合所有样品的RNA-seq质量检查结果
```bash

#################################################
# merge the rna-seq quality results
merge_length -l *length_distribution.txt -o sce
merge_saturation -l *gene_saturation.txt -o sce

cd ..
```

### 4.2 NGS文库制备中的酶切偏好
1. 测序数据中限制性内切酶酶切和连接的偏好

```bash
$ cd ./02.digestion/

#################################################
# check the reads digestion
for bam in ../01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Digest -b $bam -m 25 -M 50 --scale \
 -s ../../../1.reference/norm/gene.norm.rna.fa \
 -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done
```

2. 整合读取所有样品的消化结果

```bash
#################################################
# merge the rpf digestion
merge_digestion -l ./02.digestion/*pwm.txt -o sce

cd ..
```

### 4.3 使用`RiboParser`创建`offset`表

1. 创建RNA-seq的`offset`表

`Offset`预测在RNA-seq分析中是不必要的。一个12的常量`offset`可以分配给表中的所有条目。

```bash
$ cd ./03.offset/

#################################################
# set the offset table
for bam in ../01.qc/*.bam
do

prefix_name=$(basename $bam .bam)
rna_Offset -m 27 -M 50 -e 12 -o $prefix_name &> $prefix_name".log"

done
```

### 4.4 将`BAM`文件转换为reads密度

将`BAM`文件中的读取计数转换为密度值，并将其保存在`TXT`格式文件中。

1. 转换RNA-seq数据输出

```bash
$ cd ./04.density/

#################################################
# convert the reads to density
for bam in ../01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rna_Density -b $bam -m 27 -M 33 -l --thread 10 \
 -p ../03.offset/$prefix_name"_offset.txt" \
 -s ../../../1.reference/norm/gene.norm.rna.fa \
 -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

cd ..
```


### 4.5 为所有样本整合密度文件

不同批次和样品的数据可以集成，统一分析。如果不集成，则需要对每个样品进行单独分析，从而增加操作步骤的数量。

1. 整合所有样品的RNA-seq密度结果
```bash
$ cd ./05.merge/

#################################################
# create the samples file: RNA.file.list
merge_dst_list -l ../04.density/*_rna.txt -o RNA.file.list

cat RNA.file.list

Name File  Type
wt_rna_YPD1 /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944924_rna.txt RNA
wt_rna_YPD2 /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944925_rna.txt RNA
wt_rna_YPD3 /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944926_rna.txt RNA
ncs2d_rna_YPD1  /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944927rna.txt RNA
ncs2d_rna_YPD2  /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944928_rna.txt RNA
ncs2d_rna_YPD3  /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944929_rna.txt RNA
elp6d_rna_YPD1  /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944930_rna.txt RNA
elp6d_rna_YPD2  /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944931_rna.txt RNA
elp6d_rna_YPD3  /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944932_rna.txt RNA
ncs2d_elp6d_rna_YPD1  /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944933_rna.txt RNA
ncs2d_elp6d_rna_YPD2  /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944934_rna.txt RNA
ncs2d_elp6d_rna_YPD3  /home/sce/3.rna-seq/5.riboparser/04.density/SRR1944935_rna.txt RNA

#################################################
# merge all the RNA-seq files
rpf_Merge -l RNA.file.list -o RNA &> RNA.log

cd ..
```

### 4.6 计算三核苷酸周期

1. 检查RNA-seq的三核苷酸周期性
```bash
$ cd ./06.periodicity/

#################################################
# check the periodicity
rpf_Periodicity \
 -r ../05.merge/RNA_merged.txt \
 -m 30 --tis 0 --tts 0 -o RNA &> RNA.log
```

### 4.7 `Meta-gene`分析

利用`meta-gene`分析研究起始和终止密码子附近的reads密度。

1. RNA-seq的`meta-gene`分析

```bash
$ cd ./07.metaplot/

#################################################
# metagene analysis
rpf_Metaplot \
 -t ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RNA_merged.txt \
 -m 50 --mode bar -o RNA &> RNA.log
```

### 4.8 基因覆盖度

检查reads密度沿基因体的分布。

1. 检查RNA-seq基因密度

```bash
$ cd ./08.coverage/

#################################################
# check the reads density along with the gene body
rpf_Coverage \
 -t ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RNA_merged.txt \
 -m 50 --outlier \
 -b 10,100,10 \
 -n --heat \
 -o RNA &> RNA.log
```

### 4.9 检查样品的重复性

1. 检查RNA-seq的重复性

```bash
$ cd ./09.correlation/

#################################################
# calculate the samples replication of RNA-seq
rpf_Corr \
 -r ../05.merge/RNA_merged.txt \
 -o RNA &> RNA.log
```


## 5. 使用 `RiboParser` 对 `GSE67387` 数据集的 Ribo-seq 进行分析

### 5.0 准备用于存储结果的目录

```bash
$ cd
$ mkdir -p ./4.ribo-seq/5.riboparser/
$ cd ./4.ribo-seq/5.riboparser/
$ mkdir 01.qc 02.digestion 03.offset 04.density 05.merge \
 06.periodicity 07.metaplot 08.coverage 09.correlation 10.quantification \
 11.pausing_score 12.codon_occupancy 13.codon_decoding_time 14.codon_selection_time \
 15.coefficient_of_variation 16.meta_codon 17.shuffle 18.retrieve
```

### 5.1 测序数据的质量检查

为了确保下游分析的可靠性，核糖体分析数据的严格质量控制（QC）必须包括对峰值检测和基因检测率的系统评估。

`峰值检测`:\
计算优势RPF长度（期望范围：26-34 nt）及其占总reads的比例（QC通过阈值：期望范围≥70%）。

`基因检出率`：\
从5%到95%的数据中取样，计算被覆盖基因的数量及其表达水平。当基因计数曲线斜率趋近于0时，通常表明测序数据接近饱和。

1. `rpf_Check`的解释

```bash
$ rpf_Check -h

Check the RPFs mapping condition.

Step1: Checking the input Arguments.

usage: rpf_Check [-h] -t TRANSCRIPT -b BAM -o OUTPUT [--thread THREAD] [-g {0,1}] [-a {star,hisat2,bowtie2}] [-r] [-l] [-s]

This script is used to summary the BAM condition.

options:
  -h, --help            show this help message and exit.
  --thread THREAD       the number of threads (default: 1). Suitable for large bam files > 1G.
                        It will take a lot of memory.
  -g {0,1}              filter the number of reads mapped loci (default: 0). [0]: all reads will be
                        used; [1]: reads with unique mapped loci will be used
  -a {star,hisat2,bowtie2}
                        screen the reads mapped loci from BAM file generate with different reads alignment methods (default: star).
  -r                    reads aligned to negative strand will also be counted. (default: False).
  -l                    only keep the longest transcripts (default: False).
  -s                    whether to calculate RPF saturation. (default: False). This step will take 
                        a lot of time and memory.

Required arguments:
  -t TRANSCRIPT         the input file name of gene annotation.
  -b BAM                the input file name of bam.
  -o OUTPUT             the prefix of output file.
```

2. ribo-seq数据的质量检查

```bash
$ cd ./01.qc/

#################################################
# check the ribo-seq quality
for bam in ../../3.star/*Aligned.toTranscriptome.out.bam
do
prefix_name=$(basename $bam Aligned.toTranscriptome.out.bam)

rpf_Check -b $bam -s --thread 10 \
 -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done
```

3. `rpf_Check`的结果

```bash
# results of sample SRR1944912
SRR1944912.bam # Filtered and sorted BAM file
SRR1944912.bam.bai # Index file for the BAM file
SRR1944912_gene_saturation.pdf # Barplot for gene saturation analysis
SRR1944912_gene_saturation.png # Barplot for gene saturation analysis
SRR1944912_gene_saturation.txt # Statistical results of gene saturation analysis
SRR1944912_length_distribution.pdf # Line graph of reads length distribution
SRR1944912_length_distribution.png # Line graph of reads length distribution
SRR1944912_length_distribution.txt # Statistical results of reads length distribution
SRR1944912.log # Log file of program execution
SRR1944912_reads_saturation.pdf # Boxplot for reads saturation analysis
SRR1944912_reads_saturation.png # Boxplot for reads saturation analysis
SRR1944912_reads_saturation.txt # Statistical results of reads saturation analysis
```

4. `merge_length`和`merge_saturation`的解释

```bash
$ merge_length -h

Step1: Checking the input Arguments.
usage: merge_length [-h] -l LIST [LIST ...] -o OUTPUT

This script is used to merge the length distribution files

options:
  -h, --help            show this help message and exit

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for length distribution files (e.g., '*length_distribution.txt').
  -o OUTPUT             prefix of output file name.
```

```bash
$ merge_saturation -h

Step1: Checking the input Arguments.
usage: merge_saturation [-h] -l LIST [LIST ...] -o OUTPUT

This script is used to merge the saturation files

options:
  -h, --help            show this help message and exit

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for saturation files (e.g., '*_gene_saturation.txt').
  -o OUTPUT             prefix of output file name.
```

5. 整合所有样品的Ribo-seq质量检查结果

```bash
#################################################
# merge the ribo-seq quality results
merge_length -l *length_distribution.txt -o RIBO
merge_saturation -l *gene_saturation.txt -o RIBO

cd ..
```


### 5.2 NGS文库制备中的酶切偏好

核糖核酸酶酶切和连接步骤在Ribo-seq建库步骤中会引入偏好，影响RNA片段的代表性。

`核糖核酸酶消化偏好`：\
核糖核酸酶可能会优先切割某些RNA序列或结构，导致数据中RNA种类的不均匀断裂和扭曲。

`连接偏好`:\
连接效率可能因序列背景、RNA片段长度或二级结构而异。更有利的连接位点可能导致某些片段的过度表达。

为了减少这些偏差，重要的是分析消化和连接效率，优化协议，并包括控制以确保数据的准确性和可靠性。

1. `rpf_Digest`的解释

```bash
$ rpf_Digest -h

Step1: Checking the input Arguments.

usage: rpf_Digest [-h] -t TRANSCRIPT -s SEQUENCE -b BAM -o OUTPUT [-l] [--scale] [-m MIN] [-M MAX]

This script is used to Detect the digestion sites.

options:
  -h, --help     show this help message and exit
  -l             only retain the transcript with longest CDS of each gene (default: False).Recommended : True
  --scale        scale the motif matrix (default: False).
  -m MIN         the minimum reads length to keep (default: 20 nt).
  -M MAX         the maximum reads length to keep (default: 100 nt).

Required arguments:
  -t TRANSCRIPT  the name of input transcript file in TXT format.
  -s SEQUENCE    the name of input transcript sequence file in FA format.
  -b BAM         the name of mapping file in BAM format.
  -o OUTPUT      the name of output file. (prefix + _digestion_sites.txt)
```

2. 测序数据中限制性内切酶酶切和连接的偏好

```bash
$ cd ./02.digestion/

#################################################
# check the reads digestion
for bam in ../01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Digest -b $bam -m 27 -M 33 --scale \
 -s ../../../1.reference/norm/gene.norm.rna.fa \
 -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done
```

3. `rpf_Digest`的结果

```bash
# results of SRR1944912
SRR1944912_3end_counts.txt # Statistical values of different bases at the 3'-end of reads
SRR1944912_3end_pwm.txt # PWM matrix of bases at the 3'-end of reads
SRR1944912_3end_seqlogo2.pdf # Seqlogo of bases at the 3'-end of reads
SRR1944912_5end_counts.txt # Statistical values of different bases at the 5'-end of reads
SRR1944912_5end_pwm.txt # PWM matrix of bases at the 5'-end of reads
SRR1944912_5end_seqlogo2.pdf # Seqlogo of bases at the 5'-end of reads
SRR1944912_digestion_sites.txt # Open reading frame aligned at the ends of reads
SRR1944912.log # Log file of program execution
SRR1944912_scaled_digestion_sites_plot.pdf # Heatmap of open reading frame aligned at the ends of reads
```


4. `merge_digestion`的解释

```bash
$ merge_digestion -h

Step1: Checking the input Arguments.
usage: merge_digestion [-h] -l LIST [LIST ...] -o OUTPUT

This script is used to merge the reads digestion files

options:
  -h, --help            show this help message and exit

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for digestion files (e.g., '*_5end_pwm.txt').
  -o OUTPUT             prefix of output file name.
```

5. 整合读取所有样品的消化结果
```bash
#################################################
# merge the rpf digestion
merge_digestion -l *pwm.txt -o RIBO

cd ..
```


### 5.3 使用`RiboParser`来预测最佳偏移量

Ribo-seq分析密码子解析翻译的能力依赖于准确识别每个RPF位于核糖体`A`,`P`和`E`位点的密码子。

偏移量表示RPF的5 '端到`P`位点密码子的第一个核苷酸的距离。

常用的两种方法是基于核糖体结构的模型（`RSBM`）和基于启动/停止的模型（`SSCBM`）。

1. `rpf_Offset`的解释

```bash
$ rpf_Offset -h

Step1: Checking the input Arguments.

usage: rpf_Offset [-h] -t TRANSCRIPT -b BAM -o OUTPUT [--mode {SSCBM,RSBM}] [-a {both,tis,tts}] [-l] [-m MIN] [-M MAX] [-p EXP_PEAK] [-s SHIFT] [--silence] [-d]

This script is used to detect the P-site offset.

options:
  -h, --help           show this help message and exit
  --mode {SSCBM,RSBM}  specify the mode of offset detect [SSCBM, RSBM]. (default: SSCBM).
  -a {both,tis,tts}    specify the alignment of reads for offset detect [both, tis, tts]. (default: both).
  -l                   only retain the transcript with longest CDS of each gene (default: False).
  -m MIN               the minimum reads length to keep (default: 27 nt).
  -M MAX               the maximum reads length to keep (default: 33 nt).
  -p EXP_PEAK          expected RPFs length fitted to ribosome structure [~30 nt] (default: 30 nt).
  -s SHIFT             psite shift for different RPFs length. (default: 2 nt).
  --silence            discard the warning information. (default: True).
  -d                   output the details of offset (default: False).

Required arguments:
  -t TRANSCRIPT        the name of input transcript filein TXT format.
  -b BAM               the name of mapping file in BAM format.
  -o OUTPUT            the prefix of output file. (prefix + _offset.txt)
```

2. 预测Ribo-seq的最佳偏移量

```bash
$ cd ../03.offset/

#################################################
# predict the offset table
for bam in ../01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Offset -b $bam -m 27 -M 33 -p 30 -d \
 --mode SSCBM \
 -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done
```

3. `rpf_Offset`的结果

```bash
# results of SRR1944912
SRR1944912.log # Log file of program execution
SRR1944912_RSBM_offset.pdf # Heatmap of RSBM model calculation results
SRR1944912_RSBM_offset.png # Heatmap of RSBM model calculation results
SRR1944912_RSBM_offset.txt # Calculation results of the RSBM model
SRR1944912_SSCBM_offset.pdf # Heatmap of SSCBM model calculation results
SRR1944912_SSCBM_offset.png # Heatmap of SSCBM model calculation results
SRR1944912_SSCBM_offset_scale.pdf # Row-normalized heatmap of SSCBM model calculation results
SRR1944912_SSCBM_offset_scale.png # Row-normalized heatmap of SSCBM model calculation results
SRR1944912_SSCBM_offset.txt # Calculation results of the SSCBM model
SRR1944912_tis_3end.txt # Distribution statistics of reads 3'-end at the start codon
SRR1944912_tis_5end.txt # Distribution statistics of reads 5'-end at the start codon
SRR1944912_tts_3end.txt # Distribution statistics of reads 3'-end at the stop codon
SRR1944912_tts_5end.txt # Distribution statistics of reads 5'-end at the stop codon
```

4. `merge_offset`的解释

```bash
$ merge_offset -h

Step1: Checking the input Arguments.
usage: merge_offset [-h] -l LIST [LIST ...] -o OUTPUT

This script is used to merge the reads offset files

options:
  -h, --help            show this help message and exit

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for RSBM/SSCBM offset files (e.g., '*RSBM_offset.txt').
  -o OUTPUT             prefix of output file name (default: prefix + _offset.txt).
```

5. 对所有样本的偏移表结果进行积分

```bash
#################################################
# merge the ribo-seq offset results
merge_offset_detail -l *end.txt -o RIBO
merge_offset -l *SSCBM_offset.txt -o RIBO_SSCBM
merge_offset -l *RSBM_offset.txt -o RIBO_RSBM

cd ..
```


### 5.4 将`BAM`文件转换为reads密度

将`BAM`文件中的reads计数转换为密度值，并将其保存在`TXT`格式文件中。

1. `rpf_Density`的解释

```bash
$ rpf_Density -h

Convert reads to RPFs density.

Step1: Checking the input Arguments.

usage: rpf_Density [-h] -t TRANSCRIPT -s SEQUENCE -b BAM -p PSITE -o OUTPUT [-l] [-m MIN] [-M MAX] [--period PERIODICITY] [--silence] [--thread THREAD]

This script is used to convert Ribo-seq bam to p-site density.

options:
  -h, --help            show this help message and exit
  -l                    only retain the transcript with longest CDS of each gene (default: False).Recommended : True
  -m MIN                the minimum reads length to keep (default: 27 nt).
  -M MAX                the maximum reads length to keep (default: 33 nt).
  --silence             discard the warning information. (default: True).
  --thread THREAD       the number of threads (default: 1). It will take a lot of memory.

Required arguments:
  -t TRANSCRIPT         the name of input transcript file in TXT format.
  -s SEQUENCE           the name of input transcript sequence file in FA format.
  -b BAM                the name of mapping file in BAM format.
  -p PSITE              the name of p-site offset file in TXT format.
  -o OUTPUT             the prefix of output file. (output = prefix + _rpf.txt)
  --period PERIODICITY  the minimum 3nt periodicity to keep. (default: 40).
```

2. 转换Ribo-seq数据输出

```bash
#################################################
# convert the rpf to density
for bam in ../01.qc/*.bam
do
prefix_name=$(basename $bam .bam)

rpf_Density -b $bam -m 27 -M 33 --period 40 -l --thread 12 \
 -p ../03.offset/$prefix_name"_SSCBM_offset.txt" \
 -s ../../../1.reference/norm/gene.norm.rna.fa \
 -t ../../../1.reference/norm/gene.norm.txt \
 -o $prefix_name &> $prefix_name".log"

done

cd ..
```

3. `rpf_Density`的结果

```bash
# results of SRR1944912
SRR1944912.log # Log file of program execution
SRR1944912_rpf.txt # File contains RPFs density on each gene
```


### 5.5 为所有样本整合密度文件

不同批次和样品的数据可以集成，统一分析。如果不集成，则需要对每个样品进行单独分析，从而增加操作步骤的数量。

1. `merge_dst_list`的解释

```bash
$ merge_dst_list -h

Step1: Checking the input Arguments.
usage: merge_dst_list [-h] -l LIST [LIST ...] [-o OUTPUT]

This script is used to create the list of density files

options:
  -h, --help            show this help message and exit

Required arguments:
  -l LIST [LIST ...], --list LIST [LIST ...]
                        List for density files (e.g., '*_rpf.txt').
  -o OUTPUT             output file name (default: RIBO.file.list).
```

2. 创建示例列表

```bash
$ cd ./05.merge/

#################################################
# create the samples file: Ribo.file.list
merge_dst_list -l ../04.density/*_rpf.txt -o RIBO.file.list

cat RIBO.file.list

Name File  Type
wt_ribo_YPD1	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944912_rpf.txt Ribo
wt_ribo_YPD2	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944913_rpf.txt Ribo
wt_ribo_YPD3	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944914_rpf.txt Ribo
ncs2d_ribo_YPD1	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944915_rpf.txt Ribo
ncs2d_ribo_YPD2	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944916_rpf.txt Ribo
ncs2d_ribo_YPD3	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944917_rpf.txt Ribo
elp6d_ribo_YPD1	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944918_rpf.txt Ribo
elp6d_ribo_YPD2	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944919_rpf.txt Ribo
elp6d_ribo_YPD3	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944920_rpf.txt Ribo
ncs2d_elp6d_ribo_YPD1	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944921_rpf.txt Ribo
ncs2d_elp6d_ribo_YPD2	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944922_rpf.txt Ribo
ncs2d_elp6d_ribo_YPD3	/home/sce/4.ribo-seq/5.riboparser/04.density/SRR1944923_rpf.txt Ribo
```

3. `rpf_Merge`的解释

```bash
$ rpf_Merge -h

Merge RPFs files from different samples.

Step1: Checking the input Arguments.
usage: rpf_Merge [-h] -l LIST -o OUTPUT

This script is used to merge the density file.

options:
  -h, --help  show this help message and exit

Required arguments:
  -l LIST     the sample list in TXT format.
  -o OUTPUT   the prefix of output file. (prefix + _merged.txt)
```


4. 整合所有样品的Ribo-seq密度结果

```bash
#################################################
# merge all the Ribo-seq files
rpf_Merge -l RIBO.file.list -o RIBO &> RIBO.log

cd ..
```

5. `rpf_Merge`的结果

```bash
RIBO.log # Log file of program execution
RIBO.file.list # RPFs density file list of samples
RIBO_merged.txt # Merged RPFs density file
```

### 5.6 计算三核苷酸周期

三核苷酸周期是Ribo-seq数据分析的关键质量指标，从根本上决定了密码子解析结果的生物学可解释性。

高质量的周期性（通常为>0.6相位相干分数）是稳健密码子级分析的必要先决条件，而周期性不足（<0.45）会系统性地引入帧模糊伪像，从而影响平移测量。

读帧同步的丢失导致核糖体`A/P/E`位点的错误分配，产生帧外读聚合的假阳性暂停位点。

1. `rpf_Periodicity`的解释

```bash
$ rpf_Periodicity -h

Draw the periodicity plot.

Step1: Checking the input Arguments.

usage: rpf_Periodicity [-h] -r RPF -o OUTPUT [-t TRANSCRIPT] [-m MIN] [--tis TIS] [--tts TTS]

This script is used to draw the periodicity plot.

options:
  -h, --help     show this help message and exit
  -m MIN         retain transcript with more than minimum RPFs. (default: 50).
  --tis TIS      The number of codons after TIS will be discarded.. (default: 0 AA).
  --tts TTS      The number of codons before TTS will be discarded.. (default: 0 AA).

Required arguments:
  -r RPF         the name of input RPFs file in TXT format.
  -o OUTPUT      the prefix of output file.
  -t TRANSCRIPT  the name of input transcript filein TXT format.
```

2. 检查Ribo-seq的三核苷酸周期性

```bash
$ cd ./06.periodicity/

#################################################
# check the periodicity
rpf_Periodicity \
 -r ../05.merge/RIBO_merged.txt \
 -m 30 --tis 0 --tts 0 -o RIBO &> RIBO.log

cd ..
```

3. `rpf_Periodicity`的结果

```bash
RIBO_count_periodicity_plot.pdf # Barplot of read counts showing 3-nucleotide periodicity
RIBO_count_periodicity_plot.png # Barplot of read counts showing 3-nucleotide periodicity
RIBO.log # Log file of program execution
RIBO_periodicity.txt # Statistical values of 3-nucleotide periodicity of all samples
RIBO_ratio_periodicity_plot.pdf # Barplot of read ratios showing 3-nucleotide periodicity
RIBO_ratio_periodicity_plot.png # Barplot of read ratios showing 3-nucleotide periodicity
```


### 5.7 `Meta-gene`分析

我们的`Meta-gene`分析框架能够系统地研究近端起始和终止密码子的翻译动态。
它需要聚合核糖体保护片段（RPF）密度相对于起始/终止位点跨越-15和+60个密码子窗口，通过转录物丰度标准化。

这揭示了:\
翻译起始效率的峰/谷模式，\
Reads积累梯度反映终止动力学，\
三核苷酸周期性定量。

1. `rpf_Metaplot`的解释

```bash
$ rpf_Metaplot -h

Draw the metaplot.

Step1: Checking the input Arguments.

usage: rpf_Metaplot [-h] -t TRANSCRIPT -r RPF -o OUTPUT [-m MIN] [--utr5 UTR5] [--cds CDS] [--utr3 UTR3] [-n] [--mode {line,bar}]

This script is used to draw the meta plot.

options:
  -h, --help         show this help message and exit
  -m MIN             delete transcript with less than minimum RPFs. (default: 50).
  --utr5 UTR5        the codon number in 5-utr region (default: 20 AA).
  --cds CDS          the codon number in cds region (default: 50 AA).
  --utr3 UTR3        the codon number in 3-utr region (default: 20 AA).
  -n                 normalize the RPFs count to RPM. (default: False).
  --mode {line,bar}  specify the mode of metaplot. (default: bar).

Required arguments:
  -t TRANSCRIPT      the name of input transcript filein TXT format.
  -r RPF             the name of input RPFs file in TXT format.
  -o OUTPUT          the prefix name of output file.
```

2. Ribo-seq的`Meta-gene`分析

```bash
$ cd ./07.metaplot/

#################################################
# metagene analysis
rpf_Metaplot \
 -t ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -m 50 --mode bar -o RIBO &> RIBO.log

cd ..
```

3. `rpf_Metaplot`的结果

```bash
RIBO.log # Log file of program execution
RIBO_tis_tts_metaplot.txt # Metagene statistical values at start and stop codons across all samples
RIBO_SRR1944912_meta_bar_plot.pdf # Metaplot of all samples at start and stop codons
RIBO_SRR1944912_meta_bar_plot.png # Metaplot of all samples at start and stop codons
RIBO_SRR1944912_tis_tts_metaplot.txt # Metagene statistical values at start and stop codons for sample SRR1944912
```


### 5.8 基因覆盖度

为了系统地评估Ribo-seq数据中潜在的技术或翻译偏差，我们开发了一个分析管道，用于量化基因体间全基因组RPF覆盖均匀性。

基因被划分为固定长度的间隔（默认：CDS长度的10%），记录每个`bin`的原始RPF计数。
这些值使用每百万读数（RPM）进行库大小标准化，以解释转录本长度和测序深度的变化。

处理后的数据通过三个互补的分析视图呈现：\
`总体密度曲线`：平滑的线形图显示所有基因的平均RPF密度，突出显示全局覆盖趋势。\
`基因特异性热图`：按基因表达水平排序的规范化RPF计数（log2转换）的矩阵可视化，揭示个体基因覆盖模式。\
`分形覆盖率分布`：堆叠条形图显示在每个区间内达到阈值覆盖率（≥预期读数的50%）的基因百分比。

1. `rpf_Coverage`的解释

```bash
$ rpf_Coverage -h

Draw the metagene coverage.

Step1: Checking the input Arguments.
usage: rpf_Coverage [-h] -t TRANSCRIPT -r RPF [-o OUTPUT] [-f {0,1,2,all}] [-m MIN] [-b BIN] [-n] [--thread THREAD] [--outlier] [--set {intersect,union}] [--heat] [--bar]

This script is used to draw the coverage meta plot.

options:
  -h, --help            show this help message and exit
  -f {0,1,2,all}        set the reading frame for occupancy calculation. (default: all).
  -m MIN                retain transcript with more than minimum RPFs. (default: 50).
  -b BIN                adjust the transcript to specified bins. 30 for 5'-UTRand 3'-UTR, 100 for CDS. (default: 30,100,30).
  -n                    normalize the RPFs count to RPM. (default: False).
  --thread THREAD       the number of threads. (default: 1).
  --outlier             filter the outliers (default: False).
  --set {intersect,union}
                        filter the gene list with 5-UTR / CDS / 3-UTR. (default: union).
  --heat                draw the coverage heatmap of whole gene. (default: False).
  --bar                 draw the coverage barplot of whole gene. (default: False).

Required arguments:
  -t TRANSCRIPT         the name of input transcript filein TXT format.
  -r RPF                the name of input RPFs file in TXT format.
  -o OUTPUT             the prefix of output file.
```

2. 检查Ribo-seq的基因密度

```bash
$ cd ./08.coverage/

#################################################
# check the rpf density along with the gene body
rpf_Coverage \
 -t ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -m 50 --outlier \
 -b 10,100,10 \
 -n --heat \
 -o RIBO &> RIBO.log

rpf_Percent \
 -t ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -n -m 50 \
 -f 0 \
 -o RIBO &>> RIBO.log

cd ..
```

3. `rpf_Coverage` 的结果

```bash
RIBO.log # Log file of program execution
RIBO_SRR1944912_10_150_10_coverage.txt # Statistical results of RPF density distribution across genes
RIBO_SRR1944912_10_150_10_heat_plot.png # Heatmap of RPF density distribution across genes
RIBO_SRR1944912_coverage_bar_plot.pdf # Percentage barplot of RPF density distribution across genes
RIBO_SRR1944912_coverage_bar_plot.png # Percentage barplot of RPF density distribution across genes
RIBO_SRR1944912_coverage_line_plot.pdf # Percentage lineplot of RPF density distribution across genes
RIBO_SRR1944912_coverage_line_plot.png # Percentage lineplot of RPF density distribution across genes
```


### 5.9 检查样品的重复性

我们的分析管道提供了一个层次框架，用于评估核糖体分析研究中的样品可重复性，实现双水平相关分析。

`gene水平的再现性`:\
定量整个基因体的核糖体保护片段（RPFs），然后计算样本间Pearson相关系数。
这种方法评估全局翻译一致性，特别适用于具有强大核糖体覆盖的高表达基因。

`ORF水平的再现性`:\
在单个开放阅读框（ORF）中执行RPFs的核苷酸分辨率定量，然后计算重复样本之间的Pearson相关性。
这种细粒度分析检测局部平移变化，同时通过帧内reads滤波保持相位感知。

1. `rpf_Corr`的解释

```bash
$ rpf_Corr -h

Draw the correlation of samples.

Step1: Checking the input Arguments.

usage: rpf_Corr [-h] -r RPF -o OUTPUT

This script is used to draw the correlation of rpf density.

options:
  -h, --help  show this help message and exit

Required arguments:
  -r RPF      the name of input RPFs file in TXT format.
  -o OUTPUT   the prefix of output file. (prefix + _rpf_merged.txt)
```

2. 检查Ribo-seq的重复性

```bash
$ cd ./09.correlation/

#################################################
# calculate the samples replication of Ribo-seq
rpf_Corr \
 -r ../05.merge/RIBO_merged.txt \
 -o RIBO &> RIBO.log
```

3. `rpf_Corr`的结果

```bash
RIBO_gene_corr_f0.txt # Pearson correlation of total RPFs on gene frame 0
RIBO_gene_corr_f1.txt # Pearson correlation of total RPFs on gene frame 1
RIBO_gene_corr_f2.txt # Pearson correlation of total RPFs on gene frame 2
RIBO_gene_corr_frame.txt # Pearson correlation of RPFs across all gene frames
RIBO_gene_correlation_plot.pdf # Heatmap of Pearson correlation of RPFs across all gene frames
RIBO_gene_correlation_plot.png # Heatmap of Pearson correlation of RPFs across all gene frames
RIBO.log # Log file of program execution
RIBO_rpf_correlation_plot.pdf # Heatmap of Pearson correlation of RPFs on gene frame 0
RIBO_rpf_correlation_plot.png # Heatmap of Pearson correlation of RPFs on gene frame 0
RIBO_rpf_corr_f0.txt # Pearson correlation between RPFs on gene frame 0
RIBO_rpf_corr_f1.txt # Pearson correlation between RPFs on gene frame 1
RIBO_rpf_corr_f2.txt # Pearson correlation between RPFs on gene frame 2
RIBO_rpf_corr_frame.txt # Pearson correlation between RPFs across all gene frames
```


### 5.10 基因表达和翻译水平的量化

Ribo-seq定量与RNA-seq方法有本质区别。

RNA-seq通过转录本的读取覆盖率来量化表达，而Ribo-seq专门测量编码序列（CDS）内的核糖体占用密度。

为了减少特定区域中动态翻译造成的误差，标准分析管道排除了起始密码子后15个密码子和终止密码子前5个密码子。

为了提高精度，可以实施严格的过滤协议，仅保留帧内RPFs。

系统地去除可能代表随机噪声而不是真实伸长事件的帧外读数。

1. `rpf_Quant`的解释

```bash
$ rpf_Quant -h

Quantify the RPFs in the different region.

Step1: Checking the input Arguments.
usage: rpf_Quant [-h] -r RPF -o OUTPUT [-f {0,1,2,all}] [--tis TIS] [--tts TTS] [--utr5] [--utr3]

This script is used to quantify RPFs in the CDS region.

options:
  -h, --help      show this help message and exit
  -f {0,1,2,all}  set the reading frame for occupancy calculation. (default: all).
  --tis TIS       The number of codons after TIS will be discarded. (default: 0).
  --tts TTS       The number of codons before TES will be discarded. (default: 0).
  --utr5          quantification of 5'-utr. (default: False).
  --utr3          quantification of 3'-utr. (default: False).

Required arguments:
  -r RPF          the name of input RPFs file in TXT format.
  -o OUTPUT       the prefix of output file. (default: prefix + _rpf_quant.txt)
```

2. Ribo-seq的定量

```bash
$ cd ./10.quantification/

#################################################
# quantify the gene expression
rpf_Quant \
 -r ../05.merge/RIBO_merged.txt \
 --tis 15 \
 --tts 5 \
 -o RIBO &> RIBO.log 

cd ..
```

3. `rpf_Quant`的结果

```bash
RIBO_cds_rpm_bar_plot.pdf # Barplot of RPFs distribution in the gene CDS region
RIBO_cds_rpm_cdf_plot.pdf # Cumulative distribution plot of RPFs in the gene CDS region
RIBO_cds_rpm_heatmap.pdf # Expression heatmap of RPFs in the gene CDS region
RIBO_cds_rpm_pca_plot.pdf # Principal component analysis (PCA) plot of RPFs in the gene CDS region
RIBO_cds_rpm_pca.txt # Principal component analysis (PCA) results of RPFs in the gene CDS region
RIBO_cds_rpf_quant.txt # Statistical results of RPFs in the gene CDS region
RIBO_cds_rpkm_quant.txt # Statistical results of RPKM in the gene CDS region
RIBO_cds_rpm_quant.txt # Statistical results of RPM in the gene CDS region
RIBO_cds_tpm_quant.txt # Statistical results of TPM in the gene CDS region
RIBO.log # Log file of program execution
RIBO_total.txt # Total RPFs across all samples
```


### 5.11 计算`Codon Pausing score`

核糖体分析研究能够捕捉到核糖体停顿，因为当一个特定的密码子平均翻译速度较慢时，
核糖体占用率的增加表明在全基因组范围内有更多与该密码子相关的RPFs。

`Pausing score`是通过对基因的平均读取密度进行标准化的每核苷酸读取计数来计算的。对于密码子水平分析，
可以使用每个密码子的所有三个核苷酸位置的读取计数进行计算。最终暂停分数表示目标密码子所有出现的平均值。

1. `rpf_Pausing`的解释

```bash
$ rpf_Pausing -h

Calculate the relative codon pausing score.

Step1: Checking the input Arguments.
usage: rpf_Pausing [-h] -r RPF [-l LIST] -o OUTPUT [-s {E,P,A}] [-f {0,1,2,all}] [-b BACKGROUND] [-m MIN] [--tis TIS] [--tts TTS] [-n] [--scale {zscore,minmax}] [--stop]
                   [--fig {none,png,pdf}] [--all]

This script is used to calculate the relative pausing score.

options:
  -h, --help            show this help message and exit
  -s {E,P,A}            set the E/P/A-site for pausing calculation. (default: P).
  -f {0,1,2,all}        set the reading frame for pausing calculation. (default: all).
  -b BACKGROUND         set the codon number before and after p-site as the background. (default: 2).
  -m MIN                retain transcript with more than minimum RPFs. (default: 50).
  --tis TIS             The number of codons after TIS will be discarded. (default: 10 AA).
  --tts TTS             The number of codons before TTS will be discarded. (default: 5 AA).
  -n                    normalize the RPFs count to RPM. (default: False).
  --scale {zscore,minmax}
                        normalize the pausing score. (default: minmax).
  --stop                rmove the stop codon. (default: False).
  --fig {none,png,pdf}  draw the rpf pausing score of each gene (it will takes a lot of time). (default: none).
  --all                 output all pausing score of each gene. (default: False).

Required arguments:
  -r RPF                the name of input RPFs file in TXT format.
  -l LIST               the gene name list in TXT format. (default: whole).
  -o OUTPUT             the prefix of output file.
```

2. 计算Ribo-seq数据中的密码子水平`Pausing score`

```bash
$ cd ./11.pausing_score/

#################################################
# calculate the codon pausing score of E/P/A site
for sites in E P A
do
rpf_Pausing \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -b 0 --stop \
 -m 30 \
 -s $sites \
 -f 0 \
 --scale minmax \
 -o "$sites"_site &> "$sites"_site.log
done

cd ..
```

3. `rpf_Pausing`的结果

```bash
A_site_cds_codon_pausing_score.txt # Pausing score for each codon at the A-site in the gene CDS region
A_site_cds_pausing_score.txt # Sum of pausing score at the A-site in the gene CDS region
A_site.log # Log file of program execution
A_site_sum_codon_pausing_score.txt # Sum of pausing score for all codons at the A-site in the gene CDS region
A_site_total_pausing_heatplot.pdf # Heatmap of pausing score for all codons at the A-site in the gene CDS region
A_site_total_pausing_heatplot.png # Heatmap of pausing score for all codons at the A-site in the gene CDS region
A_site_valid_pausing_heatplot.pdf # Heatmap of pausing score for valid codons at the A-site in the gene CDS region
A_site_valid_pausing_heatplot.png # Heatmap of pausing score for valid codons at the A-site in the gene CDS region
```


### 5.12 计算`Codon occupancy`

对于全局`codon occupancy`分析，`A（P/E）`位点密码子按既定标准分配,
每个密码子的reads计数根据其各自开放阅读帧（ORF）内的平均每个密码子reads密度进行归一化。

1. `rpf_Occupancy`的解释

```bash
$ rpf_Occupancy -h

Calculate the codon occupancy.

Step1: Checking the input Arguments.
usage: rpf_Occupancy [-h] -r RPF [-l LIST] -o OUTPUT [-s {E,P,A}] [-f {0,1,2,all}] [-m MIN] [-n] [--tis TIS] [--tts TTS] [--scale {zscore,minmax}] [--stop] [--all]

This script is used to draw the codon occupancy plot.

options:
  -h, --help            show this help message and exit
  -s {E,P,A}            set the E/P/A-site for occupancy calculation. (default: P).
  -f {0,1,2,all}        set the reading frame for occupancy calculation. (default: all).
  -m MIN                retain transcript with more than minimum RPFs. (default: 30).
  -n                    normalize the RPFs count to RPM. (default: False).
  --tis TIS             The number of codons after TIS will be discarded. (default: 15 AA).
  --tts TTS             The number of codons before TTS will be discarded.. (default: 5 AA).
  --scale {zscore,minmax}
                        normalize the occupancy. (default: minmax).
  --stop                rmove the stop codon. (default: False).
  --all                 output all RPFs density. (default: False).

Required arguments:
  -r RPF                the name of input RPFs file in TXT format.
  -l LIST               the gene name list in TXT format. (default: whole).
  -o OUTPUT             the prefix of output file.
```

2. 计算Ribo-seq数据中的密码子水平`Codon occupancy`

```bash
$ cd ./12.codon_occupancy/

#################################################
# calculate the codon occupancy of E/P/A site
for sites in E P A
do
rpf_Occupancy \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -m 30 \
 -s "$sites" \
 -f 0 --stop \
 --scale minmax \
 -o "$sites"_site &> "$sites"_site.log
 
done

cd ..
```

3. `rpf_Occupancy`的结果

```bash
A_site_codon_density.txt # Codon occupancy for each codon at the A-site in the gene CDS region
A_site_codon_occupancy.txt # Codon occupancy for all codon at the A-site in the gene CDS region
A_site.log # Log file of program execution
A_site_occupancy_corrplot.pdf # Correlation heatmap of codon occupancy at the A-site
A_site_occupancy_corrplot.png # Correlation heatmap of codon occupancy at the A-site
A_site_occupancy_corr.txt # Correlation of codon occupancy at the A-site
A_site_occupancy_heatplot.pdf # Heatmap of codon occupancy at the A-site
A_site_occupancy_heatplot.png # Heatmap of codon occupancy at the A-site
A_site_occupancy_relative_heatplot.pdf # Heatmap of relative codon occupancy at the A-site
A_site_occupancy_relative_heatplot.png # Heatmap of relative codon occupancy at the A-site
A_site_occupancy_relative_lineplot.pdf # Lineplot of relative codon occupancy at the A-site
A_site_occupancy_relative_lineplot.png # Lineplot of relative codon occupancy at the A-site
```


### 5.13 计算`Codon decoding time`

除了在Ribo-seq水平上进行翻译延伸分析外，RNA-seq还可以用于校正，以解释mRNA固有状态对翻译的影响。

1. `rpf_CDT`的解释

```bash
$ rpf_CDT -h

Calculate the codon decoding time.

Step1: Checking the input Arguments.
usage: rpf_CDT [-h] --rpf RPF --rna RNA -l LIST -o OUTPUT [-s {E,P,A}] [-f {0,1,2,all}] [-m MIN] [--tis TIS] [--tts TTS] [--scale {zscore,minmax}] [--stop]

This script is used to draw the codon decoding time plot.

options:
  -h, --help            show this help message and exit
  -s {E,P,A}            set the E/P/A-site for codon decoding time calculation. (default: P).
  -f {0,1,2,all}        set the reading frame for codon decoding time calculation. (default: all).
  -m MIN                retain transcript with more than minimum RPFs. (default: 30).
  --tis TIS             The number of codons after TIS will be discarded.. (default: 15 AA).
  --tts TTS             The number of codons before TTS will be discarded.. (default: 5 AA).
  --scale {zscore,minmax}
                        normalize the codon decoding time. (default: minmax).
  --stop                rmove the stop codon. (default: False).

Required arguments:
  --rpf RPF             the name of input RPFs file in TXT format.
  --rna RNA             the name of input reads file in TXT format.
  -l LIST               the gene name list in TXT format. (default: whole).
  -o OUTPUT             the prefix of output file.
```

2. 计算Ribo-seq数据中的密码子水平`Codon decoding time`

```bash
$ cd ./13.codon_decoding_time/

#################################################
# calculate the codon decoding time of E/P/A site
for sites in E P A
do
rpf_CDT \
 -l ../../../1.reference/norm/gene.norm.txt \
 --rna ../../../3.rna-seq/5.riboparser/05.merge/RNA_merged.txt \
 --rpf ../05.merge/RIBO_merged.txt \
 --stop \
 -m 50 \
 -f 0 \
 -s $sites \
 --tis 10 \
 --tts 5 \
 -o "$sites"_site &> "$sites"_site.log

done

cd ..
```

3. `rpf_CDT`的结果

```bash
A_site_cdt_corrplot.pdf # Correlation heatmap of codon decoding time at the A-site
A_site_cdt_corrplot.png # Correlation heatmap of codon decoding time at the A-site
A_site_cdt_corr.txt # Correlation of codon decoding time at the A-site
A_site_cdt_heatplot.pdf # Heatmap of codon decoding time at the A-site
A_site_cdt_heatplot.png # Heatmap of codon decoding time at the A-site
A_site_cdt.txt # Codon decoding time for all codon at the A-site in the gene CDS region
A_site.log # Log file of program execution
```


### 5.14 计算`Codon selection time`

密码子的使用在整个基因组中是不平等的，一些密码子比其他密码子使用得更频繁，这被认为是为了提高翻译效率。

然而，先前的研究表明，翻译效率是通过基于tRNA浓度的比例密码子使用机制来优化的。

这些结果为蛋白质翻译提供了新的见解，解释了密码子的不平等使用，并强调了翻译效率的自然选择。`Codon selection time`可以用来量化这一过程。

1. `rpf_CST`的解释

```bash
$ rpf_CST -h

Calculate the codon decoding time.

Step1: Checking the input Arguments.
usage: rpf_CST [-h] --rpf RPF --rna RNA [-l LIST] -o OUTPUT [-s {E,P,A}] [-f {0,1,2,all}] [-m MIN] [-t TIMES] [--tis TIS] [--tts TTS] [--scale {zscore,minmax}] [--stop]

This script is used to draw the codon decoding time plot.

options:
  -h, --help            show this help message and exit
  -s {E,P,A}            set the E/P/A-site for codon decoding time calculation. (default: P).
  -f {0,1,2,all}        set the reading frame for codon decoding time calculation. (default: all).
  -m MIN                retain transcript with more than minimum RPFs. (default: 30).
  -t TIMES              Specify the number of iteration times required for computation. (default: 10).
  --tis TIS             The number of codons after TIS will be discarded. (default: 0 AA).
  --tts TTS             The number of codons before TTS will be discarded. (default: 0 AA).
  --scale {zscore,minmax}
                        normalize the codon selection time. (default: minmax).
  --stop                rmove the stop codon. (default: False).

Required arguments:
  --rpf RPF             the name of input RPFs file in TXT format.
  --rna RNA             the name of input reads file in TXT format.
  -l LIST               the gene name list in TXT format. (default: whole).
  -o OUTPUT             the prefix of output file.
```

2. 计算Ribo-seq数据中的密码子水平选择时间`selection time`

```bash
$ cd ./14.codon_selection_time/

#################################################
# calculate the codon selection time of E/P/A site
for sites in E P A
do
rpf_CST \
 -l ../../../1.reference/norm/gene.norm.txt \
 --rna ../../../3.rna-seq/5.riboparser/05.merge/RNA_merged.txt \
 --rpf ../05.merge/RIBO_merged.txt \
 --stop \
 -m 50 \
 -f 0 \
 -s $sites \
 --tis 10 \
 --tts 5 \
 -o "$sites"_site &> "$sites"_site.log

done

cd ..
```

3. `rpf_CST`的结果

```bash
A_site_codon_selection_time.txt # Codon selection time for all codon at the A-site in the gene CDS region
A_site_cst_corrplot.pdf # Correlation heatmap of codon selection time at the A-site
A_site_cst_corrplot.png # Correlation heatmap of codon selection time at the A-site
A_site_cst_corr.txt # Correlation of codon selection time at the A-site
A_site_cst_heatplot.pdf # Heatmap of codon selection time at the A-site
A_site_cst_heatplot.png # Heatmap of codon selection time at the A-site
A_site_iterative_codon_selection_time.txt # Codon selection time for all codon at the A-site in the gene CDS region
A_site.log # Log file of program execution
```


### 5.15 计算`Coefficient of Variation`

虽然密码子水平的分析可以表明翻译延伸的改变，但由于基因表达动力学的变化，核糖体分析覆盖范围使基因水平的验证变得复杂。
由于噪声与覆盖率成反比，低覆盖率基因可能人为地表现出更多的翻译停顿。 

为了解决这个问题，Peter等人开发了一种基因水平的分析方法，该方法明确地模拟了噪声对覆盖率的依赖。他们将以下双参数模型拟合到数据中，以适应计数噪声的各种统计行为：

$$log_2(CV) = \frac{1}{2} log_2 (\frac{β}{μ} + α)$$

其中`CV`为给定基因核糖体谱的变异系数，`μ`为平均覆盖率（RPF读取过密码子），`α`和`β`为拟合参数。
重要的是，当`α = 0`和`β = 1`时，该方程为泊松分布；`α > 0`和`β = 1`时，该方程为负二项分布。

1. `rpf_CoV`的解释

```bash
$ rpf_CoV -h

Calculate CoV in the CDS region.

Step1: Checking the input Arguments.
usage: rpf_CoV [-h] -r RPF [-g GROUP] [-l LIST] -o OUTPUT [-f {0,1,2,all}] [-m MIN] [-n] [--tis TIS] [--tts TTS] [--fig]

This script is used to calculate CoV in the CDS region.

options:
  -h, --help      show this help message and exit
  -f {0,1,2,all}  set the reading frame for occupancy calculation. (default: all).
  -m MIN          retain transcript with more than minimum RPFs. (default: 5).
  -n              normalize the RPFs count to RPM. (default: False).
  --tis TIS       The number of codons after TIS will be discarded. (default: 15).
  --tts TTS       The number of codons before TES will be discarded. (default: 5).
  --fig           show the figure. (default: False).

Required arguments:
  -r RPF          the name of input RPFs file in TXT format.
  -g GROUP        specify the list of sample group. (default: None)
  -l LIST         the gene name list in TXT format. (default: whole).
  -o OUTPUT       the prefix of output file. (prefix + _Cov.txt)

The group file needs to contain at least two column:
+----------+---------+
| name     | group  |
+==========+=========+
| wt1      | wt      |
| wt2      | wt      |
| treat1   | treat   |
| treat2   | treat   |
| ko1      | ko      |
| ko2      | ko      |
+----------+---------+
```

2. 计算Ribo-seq数据中的基因变异系数

```bash
$ cd ./15.coefficient_of_variation/

#################################################
# Here we can configure the design file to calculate differences between different groups.
$ cat design.txt
Name	Group
WT_ribo_YPD1	WT_ribo_YPD
WT_ribo_YPD2	WT_ribo_YPD
WT_ribo_YPD3	WT_ribo_YPD
ncs2d_ribo_YPD1	ncs2d_ribo_YPD
ncs2d_ribo_YPD2	ncs2d_ribo_YPD
ncs2d_ribo_YPD3	ncs2d_ribo_YPD
elp6d_ribo_YPD1	elp6d_ribo_YPD
elp6d_ribo_YPD2	elp6d_ribo_YPD
elp6d_ribo_YPD3	elp6d_ribo_YPD
ncs2d_elp6d_ribo_YPD1	ncs2d_elp6d_ribo_YPD
ncs2d_elp6d_ribo_YPD2	ncs2d_elp6d_ribo_YPD
ncs2d_elp6d_ribo_YPD3	ncs2d_elp6d_ribo_YPD

#################################################
# calculate the coefficient of variation
rpf_CoV \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -f 0 \
 -m 30 \
 --tis 10 \
 --tts 5 \
 --fig \
 -g design.txt \
 -o RIBO &> RIBO.log

cd ..
```
3. `rpf_Cumulative_CoV` 的解释

```bash
$ rpf_Cumulative_CoV -h

Retrieve the RPFs with gene list.

Step1: Checking the input Arguments.

usage: rpf_Cumulative_CoV [-h] -r RPF [-o OUTPUT] [-l LIST] [-m MIN] [-n] [-t TRIM] [-s] [-z]

This script is used to calculate the cumulative CoV.

options:
  -h, --help  show this help message and exit
  -l LIST     the list of input genes for transcript id.
  -m MIN      retain transcript with more than minimum RPFs (default: 0).
  -n          normalize the RPFs count to RPM (default: False).
  -t TRIM     trim transcript with specific length (default: 50 nt).
  -s          split gene rpf to each TXT file (default: False).
  -z          set the start site to zero (default: False).

Required arguments:
  -r RPF      the name of input RPFs density file in TXT format.
  -o OUTPUT   prefix of output file name (default: filename + '_cumulative_CoV.txt'.
```


4. 计算Ribo-seq数据中的基因变异系数的累积分布

```bash
#################################################
# calculate the cumulative coefficient of variation
rpf_Cumulative_CoV \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -m 50 \
 -n \
 -z \
 -t 300 \
 -o RIBO &>> RIBO.log
```

5. `rpf_CoV` 的结果

```bash
gene_compared_CoV.txt # Compared gene coefficient of variation
gene_CoV.txt # Gene coefficient of variation
gene.log # Log file of program execution
gene_WT_ribo_YPD_vs_ncs2d_ribo_YPD_CoV_fitplot.pdf # Fitted lineplot of gene coefficient of variation
gene_WT_ribo_YPD_vs_ncs2d_ribo_YPD_CoV_fitplot.png # Fitted lineplot of gene coefficient of variation
```


### 5.16 `Meta-codon`分析

为了演示不同密码子的暂停模式，这里使用`Meta-plot`来显示指定密码子周围20个核苷酸窗口内的RPF密度。

该程序还支持不同帧的可视化，并包含平滑功能，以减轻由数据不稳定引起的尖峰信号。

1. `rpf_Meta_Codon`的解释

```bash
$ rpf_Meta_Codon -h

Draw the meta-codon plot.

Step1: Checking the input Arguments.
usage: rpf_Meta_Codon [-h] [-l LIST] -r RPF [-c CODON] -o OUTPUT [-f {0,1,2}] [-a AROUND] [-m MIN] [--tis TIS] [--tts TTS] [-n] [-u] [-s] [--smooth SMOOTH] [--thread THREAD] [--fig]

This script is used to draw the meta codon plot.

options:
  -h, --help       show this help message and exit
  -f {0,1,2}       set the reading frame for occupancy calculation. (default: all).
  -a AROUND        retrieve length of codon upstream and downstream. (default: 20).
  -m MIN           retain transcript with more than minimum RPFs. (default: 50).
  --tis TIS        The number of codons after TIS will be discarded.. (default: 0 AA).
  --tts TTS        The number of codons before TTS will be discarded.. (default: 0 AA).
  -n               normalize the RPFs count to RPM. (default: False).
  -u               delete the cross repetition codon in different window. (default: False).
  -s               scale the window density with gene density. (default: False).
  --smooth SMOOTH  smooth the window density [eg, 3,1]. (default: None).
  --thread THREAD  the number of threads (default: 1).
  --fig            output the figure. (default: False).

Required arguments:
  -l LIST          the gene name list in TXT format. (default: whole).
  -r RPF           the name of input RPFs file in TXT format.
  -c CODON         the codon list in TXT format.
  -o OUTPUT        the prefix of output file.
```

2. 计算Ribo-seq数据的`meta-codon`密度

```bash
$ cd ./16.meta_codon/

#################################################
# Here we can configure the codon list.
$ cat codon_list.txt
AAA
AAC
AAG
AAT
AAGAAG
ATGATG
CCCGGG
...

#################################################
# codon meta analysis
rpf_Meta_Codon \
 -r ../05.merge/RIBO_merged.txt \
 -m 50 -f 0 \
 -c codon_list.txt \
 -a 15 -u -n \
 -o RIBO &> RIBO.log

cd ..
```

3. `rpf_Meta_Codon`的结果

```bash
RIBO.log # Log file of program execution
RIBO_AAA_97591_8146_meta_density.txt # RPFs density of AAA codon
RIBO_AAA_97591_8146_meta_sequence.txt # Upstream and downstream sequence around AAA codon
RIBO_AAA.pdf # Metaplot of AAA codon
RIBO_AAA.png # Metaplot of AAA codon
```

## 6. smORF 鉴定
### 6.1 全转录组扫描 ORF

1. `smorf_scanner` 的解释

使用全部的转录组序列，从头扫描所有潜在的 ORF，结果会输出所有已知和未知的 ORF。

```bash
$ smorf_scanner -h

usage: smorf_scanner 
[-h] -g GENOME -a ANNOTATION [-o OUT_PREFIX] [--orf-prefix ORF_PREFIX] [--start-codons START_CODONS]
[--min-aa MIN_AA] [--max-aa MAX_AA] [--scan-strand {sense,antisense,both}] [--kozak-up KOZAK_UP]
[--kozak-down KOZAK_DOWN] [-t THREADS] [--mark-overlap] [--remove-discarded] [--include-stop]

Scan transcript-centric smORFs from genome FASTA and genePred annotation.

options:
  -h, --help            show this help message and exit
  -g, --genome GENOME   Input genome FASTA file.
  -a, --annotation ANNOTATION
                        Input genePred annotation file.
  -o, --out-prefix OUT_PREFIX
                        Output prefix.
  --orf-prefix ORF_PREFIX
                        Prefix for ORF IDs.
  --start-codons START_CODONS
                        Comma-separated start codons, such as ATG,CTG,GTG,TTG.
  --min-aa MIN_AA       Minimum ORF length in amino acids.
  --max-aa MAX_AA       Maximum ORF length in amino acids.
  --scan-strand {sense,antisense,both}
                        Scan sense, antisense, or both strands.
  --kozak-up KOZAK_UP   Number of upstream nucleotides for Kozak sequence.
  --kozak-down KOZAK_DOWN
                        Number of downstream nucleotides after start codon for Kozak sequence.
  -t, --threads THREADS
                        Number of worker processes for parallel ORF scanning.
  --mark-overlap        Mark nested or overlapping ORFs.
  --remove-discarded    Remove same-frame internal ORFs.
  --include-stop        Keep stop codon symbol in peptide sequence.

# example
smorf_scanner \
  --genome ../genome/GCF_mine_genomic.fna \
  --annotation ../norm/mine.genepred \
  --out-prefix mine \
  --start-codons ATG \
  --min-aa 8 \
  --max-aa 10000 \
  --scan-strand both \
  --kozak-up 6 \
  --kozak-down 6 \
  --mark-overlap \
  --threads 20

```


2. 使用 `smorf_filter` 对扫描的 ORF 过滤

全转录组扫描的结果仅基于开放阅读框的位置，因此存在大量的假阳性结果。所以需要进行假阳性的过滤。
同时可以筛选一些 smORF 用于研究。

```bash
$ smorf_filter -h
usage: smorf_filter 
[-h] -i INPUT [-o OUT_PREFIX] [--keep-start-codons KEEP_START_CODONS] [--min-aa MIN_AA] [--max-aa MAX_AA]
[--keep-categories KEEP_CATEGORIES] [--remove-categories REMOVE_CATEGORIES] [--keep-antisense] [--keep-secondary]
[--keep-partial] [--kozak-mode {none,annotated,builtin,pwm,sequence}]
[--builtin-kozak {arabidopsis,drosophila,maize,plant,rice,terrestrial_plant,vertebrate,yeast}]
[--kozak-pwm KOZAK_PWM] [--kozak-seq KOZAK_SEQ] [--annotated-categories ANNOTATED_CATEGORIES]
[--min-annotated-kozak MIN_ANNOTATED_KOZAK]
[--fallback-builtin-kozak {arabidopsis,drosophila,maize,plant,rice,terrestrial_plant,vertebrate,yeast}]
[--no-kozak-fallback] [--min-kozak-pwm-score MIN_KOZAK_PWM_SCORE] [--export-kozak-pwm EXPORT_KOZAK_PWM]
[--list-builtin-kozak]

Filter ORF.message.txt using rule-based criteria and Kozak PWM scoring.

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     Input ORF.message.txt file.
  -o, --out-prefix OUT_PREFIX
                        Output prefix. (default: ORF.filtered).
  --keep-start-codons KEEP_START_CODONS
                        Comma-separated start codons to keep. (default: ATG,CTG,GTG,TTG).
  --min-aa MIN_AA       Minimum ORF peptide length. (default: 8).
  --max-aa MAX_AA       Maximum ORF peptide length. (default: 10000).
  --keep-categories KEEP_CATEGORIES
                        Comma-separated ORF categories to keep. (default:
                        uORF,dORF,lncORF,iORF,emORF,overlap_uORF,overlap_dORF,other_ORF,annotated_ORF).
  --remove-categories REMOVE_CATEGORIES
                        Comma-separated ORF categories to remove. (default: same_frame_iORF,antisense_ORF).
  --keep-antisense      Keep antisense ORFs. (default: False).
  --keep-secondary      Keep secondary ORFs. (default: False).
  --keep-partial        Keep incomplete ORFs. (default: False).
  --kozak-mode {none,annotated,builtin,pwm,sequence}
                        Kozak PWM mode. (default: annotated).
  --builtin-kozak {arabidopsis,drosophila,maize,plant,rice,terrestrial_plant,vertebrate,yeast}
                        Built-in Kozak PWM name. (default: plant).
  --kozak-pwm KOZAK_PWM
                        Custom Kozak PWM matrix file. (default: None).
  --kozak-seq KOZAK_SEQ
                        Aligned Kozak sequence file used to build PWM. (default: None).
  --annotated-categories ANNOTATED_CATEGORIES
                        Categories used to build annotated ORF Kozak PWM. (default: annotated_ORF).
  --min-annotated-kozak MIN_ANNOTATED_KOZAK
                        Minimum annotated ORF Kozak sequences required to build PWM. (default: 100).
  --fallback-builtin-kozak {arabidopsis,drosophila,maize,plant,rice,terrestrial_plant,vertebrate,yeast}
                        Fallback built-in Kozak PWM if annotated mode fails. (default: plant).
  --no-kozak-fallback   Do not fallback to built-in PWM if annotated PWM construction fails. (default: False).
  --min-kozak-pwm-score MIN_KOZAK_PWM_SCORE
                        Minimum normalized Kozak PWM score. Range: 0-1. (default: 0.0).
  --export-kozak-pwm EXPORT_KOZAK_PWM
                        Export loaded or constructed Kozak PWM. (default: None).
  --list-builtin-kozak  List built-in Kozak PWM models and exit.

# example
smorf_filter \
 -i mine.message.txt \
 -o mine.reliable \
 --min-aa 8 \
 --max-aa 10000 \
 --kozak-mode annotated \
 --keep-categories uORF,dORF,lncORF,overlap_uORF,overlap_dORF

```

3. 使用 `smorf_evidence` 检查 smORF 的翻译情况

通常情况下，稳定存在的 smORF 可以被灵敏的 Ribo-seq 捕获到，因此可以使用 Ribo-seq 的数据对 smORF 进行检查，搜索可靠的 smORF。这里兼容多个数据的输入，用以检测不同时间稳定存在的 smORF。

```bash
$ smorf_evidence -h
usage: smorf_evidence 
[-h] -i ORF_TABLE -o OUTPUT [--genepred GENEPRED] [--chrom-sizes CHROM_SIZES] [--density-list DENSITY_LIST]
[--density-plus DENSITY_PLUS] [--density-minus DENSITY_MINUS] [--density DENSITY] [--sample SAMPLE]
[--density-format {auto,wig,bedgraph}] [--coord-mode {0based-half-open,1based-closed}]
[--post-stop-codons POST_STOP_CODONS] [--pseudocount PSEUDOCOUNT] [--min-rpf-sum MIN_RPF_SUM]
[--min-covered-codon MIN_COVERED_CODON] [--min-coverage-ratio MIN_COVERAGE_RATIO]
[--strong-periodicity STRONG_PERIODICITY] [--moderate-periodicity MODERATE_PERIODICITY]
[--strong-start-pause STRONG_START_PAUSE] [--moderate-start-pause MODERATE_START_PAUSE]
[--strong-stop-pause STRONG_STOP_PAUSE] [--moderate-stop-pause MODERATE_STOP_PAUSE]
[--strong-release STRONG_RELEASE] [--moderate-release MODERATE_RELEASE]
[--uniform-coverage-ratio UNIFORM_COVERAGE_RATIO] [--uniform-gini UNIFORM_GINI]
[--uniform-max-to-mean UNIFORM_MAX_TO_MEAN] [--skewed-max-to-mean SKEWED_MAX_TO_MEAN]
[--skewed-top-fraction SKEWED_TOP_FRACTION] [--disperse-coverage-ratio DISPERSE_COVERAGE_RATIO]
[--progress-every PROGRESS_EVERY] [--keep-no-evidence]

Evaluate smORF translation evidence using Ribo-seq P-site density.

options:
  -h, --help            show this help message and exit
  -i, --orf-table ORF_TABLE
                        Filtered smORF table from smorf_filter.
  -o, --output OUTPUT   Output evidence table in TSV format.
  --genepred GENEPRED   Optional genePred file for ORF exon blocks. (default: None)
  --chrom-sizes CHROM_SIZES
                        Optional two-column chromosome size file. If not provided, chromosome sizes will be inferred from density files.
                        (default: None)
  --density-list DENSITY_LIST
                        TSV with columns: sample, strand, path, optional format. (default: None)
  --density-plus DENSITY_PLUS
                        Plus-strand P-site density file. (default: None)
  --density-minus DENSITY_MINUS
                        Minus-strand P-site density file. (default: None)
  --density DENSITY     Unstranded P-site density file. (default: None)
  --sample SAMPLE       Sample name for direct density input. (default: sample1)
  --density-format {auto,wig,bedgraph}
                        Density file format. (default: auto)
  --coord-mode {0based-half-open,1based-closed}
                        Coordinate mode for ORF table and genePred-like blocks. (default: 0based-half-open)
  --post-stop-codons POST_STOP_CODONS
                        Number of codons after stop codon used for release signal. (default: 10)
  --pseudocount PSEUDOCOUNT
                        Pseudocount for ratio calculation. (default: 0.1)
  --min-rpf-sum MIN_RPF_SUM
                        Minimum ORF-level RPF sum for evidence scoring. (default: 3.0)
  --min-covered-codon MIN_COVERED_CODON
                        Minimum covered codon count. (default: 2)
  --min-coverage-ratio MIN_COVERAGE_RATIO
                        Minimum nucleotide-level coverage ratio. (default: 0.1)
  --strong-periodicity STRONG_PERIODICITY
                        Frame-0 ratio threshold for strong periodicity. (default: 0.7)
  --moderate-periodicity MODERATE_PERIODICITY
                        Frame-0 ratio threshold for moderate periodicity. (default: 0.55)
  --strong-start-pause STRONG_START_PAUSE
                        Start pausing ratio threshold for strong signal. (default: 1.5)
  --moderate-start-pause MODERATE_START_PAUSE
                        Start pausing ratio threshold for moderate signal. (default: 1.2)
  --strong-stop-pause STRONG_STOP_PAUSE
                        Pre-stop pausing ratio threshold for strong signal. (default: 1.5)
  --moderate-stop-pause MODERATE_STOP_PAUSE
                        Pre-stop pausing ratio threshold for moderate signal. (default: 1.2)
  --strong-release STRONG_RELEASE
                        Release ratio threshold for strong signal. (default: 3.0)
  --moderate-release MODERATE_RELEASE
                        Release ratio threshold for moderate signal. (default: 1.5)
  --uniform-coverage-ratio UNIFORM_COVERAGE_RATIO
                        Coverage ratio threshold for Uniform shape. (default: 0.4)
  --uniform-gini UNIFORM_GINI
                        Gini threshold for Uniform shape. (default: 0.5)
  --uniform-max-to-mean UNIFORM_MAX_TO_MEAN
                        Max/mean threshold for Uniform shape. (default: 5.0)
  --skewed-max-to-mean SKEWED_MAX_TO_MEAN
                        Max/mean threshold for Skewed shape. (default: 10.0)
  --skewed-top-fraction SKEWED_TOP_FRACTION
                        Top 10 percent density fraction threshold for Skewed shape. (default: 0.7)
  --disperse-coverage-ratio DISPERSE_COVERAGE_RATIO
                        Coverage ratio threshold below which coverage is Disperse. (default: 0.2)
  --progress-every PROGRESS_EVERY
                        Print progress every N ORFs per chromosome. (default: 10000)
  --keep-no-evidence    Keep ORFs with no RPF evidence in output. (default: False)


# example
smorf_evidence \
  -i mine.reliable.passed.message.txt \
  --genepred mine.genePred \
  --density-list ribo.bedgraph.list \
  -o mine.smorf.riboseq_evidence.txt

# cat ribo.bedgraph.list
# sample	strand	path	format
# ribo1	+	/project/mine/ribo/bedgraph/ribo1_plus.rpf.bedgraph	bedgraph
# ribo1	-	/project/mine/ribo/bedgraph/ribo1_minus.rpf.bedgraph	bedgraph
# ribo2	+	/project/mine/ribo/bedgraph/ribo2_plus.rpf.bedgraph	bedgraph
# ribo2	-	/project/mine/ribo/bedgraph/ribo2_minus.rpf.bedgraph	bedgraph

# 上面的 bedgraph 文件通过命令 rpf_Bam2bw 生成。

```

4. 使用 `smorf_integrate` 整合所有高置信度的 smORF

把通过了不同样本 Ribo-seq 数据过滤的结果进行整合，合并为一张大表。
后续可以通过对输出表格和第一步扫描的注释文件的筛选，然后通过 RiboParser 的 rpf 函数模块进行完整的质控和定量分析。

```bash
$ smorf_integrate -h
usage: smorf_integrate 
[-h] -i INPUT [--output-matrix OUTPUT_MATRIX] [--output-integrated OUTPUT_INTEGRATED] 
[--capture-labels CAPTURE_LABELS] [--pass-labels PASS_LABELS] [--excellent-min-samples EXCELLENT_MIN_SAMPLES]

Integrate long-format smORF Ribo-seq evidence into matrices and ORF-level summary tables.

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     Long-format smORF Ribo-seq evidence table generated by smorf_evidence.py.
  --output-matrix OUTPUT_MATRIX
                        Output ORF-by-sample matrix table for frame density and selected sample metrics. (default: None)
  --output-integrated OUTPUT_INTEGRATED
                        Output integrated ORF-level evidence table. (default: None)
  --capture-labels CAPTURE_LABELS
                        Translation evidence labels used to define captured samples. (default: LowConfidence,MediumConfidence,HighConfidence)
  --pass-labels PASS_LABELS
                        Translation evidence labels used to define reliable/pass samples. (default: MediumConfidence,HighConfidence)
  --excellent-min-samples EXCELLENT_MIN_SAMPLES
                        Minimum number of pass samples required to mark an ORF as Excellent. (default: 2)

# example
smorf_integrate \
 -i mine.smorf.riboseq_evidence.txt \
 --output-matrix mine.smorf.frame_density_matrix.txt \
 --output-integrated mine.smorf.integrated_evidence.txt

```




## 7. 其他工具包
### 7.1 Data shuffling

一些分析过程需要随机分配的数据进行控制，所以这里添加了一个步骤来重新排列rfp密度文件。

1. `rpf_Shuffle`的解释

```bash
$ rpf_Shuffle -h

Shuffle the RPFs data.

Step1: Checking the input Arguments.

usage: rpf_Shuffle [-h] -r RPF -o OUTPUT [-l LIST] [-s SEED] [-i]

This script is used to shuffle the RPFs data.

options:
  -h, --help  show this help message and exit
  -l LIST     the list of input genes for transcript id.
  -s SEED     the random seed for shuffle. (default: 0).
  -i          shuffle the RPFs data for each samples. (default: False).

Required arguments:
  -r RPF      the name of input RPFs density file in TXT format.
  -o OUTPUT   the prefix of output file. (prefix + _shuffle.txt)
```

2. Shuffle Ribo-seq数据的基因密度值

```bash
$ cd
$ cd ./sce/4.ribo-seq/5.riboparser/17.shuffle/

#################################################
# codon meta analysis
rpf_Shuffle \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -s 0 \
 -i \
 -o RIBO &> RIBO.log
```

3. Shuffle RNA-seq数据的基因密度值

```bash
$ cd
$ cd ./sce/3.rna-seq/5.riboparser/10.shuffle/

#################################################
# retrieve and format the gene density
rpf_Shuffle \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RNA_merged.txt \
 -s 0 \
 -i \
 -o RNA &> RNA.log
```

4. `rpf_Shuffle`的结果

```bash
# results of ribo-seq
RIBO.log # Log file of program execution
RIBO_shuffle.txt # Shuffled RPFs density file
```


### 7.2 检索并格式化基因密度

在许多情况下，需要对rfp密度文件中的基因集进行一些额外的操作。

如过滤、RPM标准化、长宽数据格式转换等。

这里为这些操作提供了专门的工具。

1. `rpf_Retrieve`的解释

```bash
$ rpf_Retrieve -h

Retrieve the RPFs with gene list.

Step1: Checking the input Arguments.

usage: rpf_Retrieve [-h] -r RPF [-o OUTPUT] [-l LIST] [-m MIN] [-n] [-f] [-s]

This script is used to retrieve density files.

options:
  -h, --help  show this help message and exit
  -l LIST     the list of input genes for transcript id.
  -m MIN      retain transcript with more than minimum RPFs (default: 0).
  -n          normalize the RPFs count to RPM (default: False).
  -f          melt three column data of each sample to one column (default: False).
  -s          split gene rpf to each TXT file (default: False).

Required arguments:
  -r RPF      the name of input RPFs density file in TXT format.
  -o OUTPUT   prefix of output file name (default: filename + '_retrieve.txt'.
```

2. 从Ribo-seq数据中提取并格式化基因密度

```bash
$ cd
$ cd ./sce/4.ribo-seq/5.riboparser/18.retrieve/

#################################################
# retrieve and format the gene density with gene covered more than 50 reads
rpf_Retrieve \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 -m 50 \
 -f \
 -n \
 -o RIBO &> RIBO.log

cd ..
```

3. 从RNA-seq数据中提取并格式化基因密度

```bash
$ cd
$ cd ./sce/3.rna-seq/5.riboparser/11.retrieve/

#################################################
# retrieve and format the gene density with gene covered more than 50 reads
rpf_Retrieve \
 -l ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RNA_merged.txt \
 -m 50 \
 -f \
 -n \
 -o RNA &> RNA.log

cd ..
```

4. `rpf_Retrieve`的结果

```bash
# results of ribo-seq
RIBO.log
RIBO_retrieve.txt
```

### 7.3 过滤移码基因

当核糖体在mRNA序列中移动一个或多个核苷酸时，翻译中的移码发生，导致密码子的误读。
这导致下游的氨基酸序列完全改变，通常导致过早终止或无功能的蛋白质。

由于突变或翻译错误，帧移可以自然发生，它们可以显著影响蛋白质的功能。

当一个稳定的移码发生时，我们可以从Ribo-seq数据中检测到不同读帧的核糖体占用。

1. `rpf_Shift`的解释

```bash
$ rpf_Shift -h

Draw the frame shifting plot.

Step1: Checking the input Arguments.

usage: rpf_Shift.py [-h] -r RPF -o OUTPUT [-t TRANSCRIPT] [-p PERIOD] [-m MIN] [--tis TIS] [--tts TTS]

This script is used to draw the frame shift plot.

options:
  -h, --help     show this help message and exit
  -p PERIOD      the minimum in-frame value for frame shifting screen, range [0 - 100]. (default: 45).
  -m MIN         retain transcript with more than minimum RPFs. (default: 50).
  --tis TIS      the number of codons after TIS will be discarded.. (default: 0 AA).
  --tts TTS      the number of codons before TTS will be discarded.. (default: 0 AA).

Required arguments:
  -r RPF         the name of input RPFs file in TXT format.
  -o OUTPUT      the prefix of output file.
  -t TRANSCRIPT  the name of input transcript filein TXT format.

```

2. 过滤移码基因

```bash
$ cd
$ cd ./sce/4.ribo-seq/5.riboparser/19.frame_shift/

#################################################
# filter the frame shifting genes
rpf_Shift \
 -t ../../../1.reference/norm/gene.norm.txt \
 -r ../05.merge/RIBO_merged.txt \
 --tis 5 --tts 5 \
 -m 50 \
 -p 45 \
 -o RIBO &> RIBO.log

cd ..
```

4. `rpf_Shift`的结果

```bash
# results of ribo-seq
RIBO.log
RIBO_gene_frame_shift_count_plot.pdf
RIBO_gene_frame_shift_count_plot.png
RIBO_gene_frame_shift_count.txt
RIBO_gene_periodicity.txt
RIBO_SRR1944912_gene_frame_shift.txt
```

## 8. shell 包装流程

### 8.0 为项目准备目录结构和实验设计文件

1. 创建用于存储原始数据和分析结果的目录结构

```bash
$ cd
$ mkdir sce
$ cd ./sce/
$ mkdir -p ./sce/1.reference
$ mkdir -p ./sce/2.rawdata/ribo-seq ./sce/2.rawdata/rna-seq
$ mkdir -p ./sce/3.rna-seq
$ mkdir -p ./sce/4.ribo-seq
```

2. 为您的 RNA-seq 和 Ribo-seq 数据准备设计文件

设计文件必须至少包含两列：`Name` 和 `Group`，如果有其他信息，可以在后面添加列。
该文件内容可保存 excel 格式用于 `RiboShiny` 的分析使用。

```bash
$ cat design.txt

Name    Group
SRR1944912      WT_ribo_YPD
SRR1944913      WT_ribo_YPD
SRR1944914      WT_ribo_YPD
SRR1944915      ncs2d_ribo_YPD
SRR1944916      ncs2d_ribo_YPD
SRR1944917      ncs2d_ribo_YPD
```


### 8.1 run_step1.sh

此步骤用于构建数据库，这是进行 reads 比对及后续使用 `RiboParser` 进行分析的关键步骤。  

该步骤适用于大多数来自 `NCBI` 的基因组和基因注释文件。  

**注意**：`Step 0.0` 中的下载文件需要根据您的项目所使用的 `物种` 进行修改.

```bash
$ nohup sh run_step1.sh &
```

### 8.2 run_step2.sh

此步骤用于分析 `RNA-seq` 数据，包括数据清洗、比对和表达定量。  

**注意**：`Step 0.0` 中的 `接头` 信息需要根据您的项目所使用的测序方法进行修改.

```bash
$ nohup sh run_step2.sh &
```

### 8.3 run_step3.sh

此步骤用于分析 `Ribo-seq` 数据，包括数据清洗、比对和表达定量。  

**注意**：`Step 0.0` 中的 `接头` 信息需要根据您的项目所使用的测序方法进行修改.

```bash
$ nohup sh run_step3.sh &
```

### 8.4 run_step4.sh

此步骤用于分析 `RNA-seq` 数据，利用 `RiboParser` 检查 `RNA-seq` 数据的测序质量，并准备格式化文件，以便与 `Ribo-seq` 进行后续联合分析。  

**注意**：`Step 1.0` 中的 `BAM` 文件和 `参考基因组信息` 文件可能需要根据您的项目文件进行修改.

```bash
$ nohup sh run_step4.sh &
```

### 8.5 run_step5.sh

此步骤用于分析 `Ribo-seq` 数据，利用 `RiboParser` 检查 `Ribo-seq` 数据的测序质量。

**注意**：`Step 0.0` 中的 `BAM` 文件、`parameters` 和 `reference genome information` 文件可以
需要根据为您的项目定义的文件进行修改.

```bash
$ nohup sh run_step5.sh &
```


## 9. Computational performance of the RiboParser
我们在CentOS 7系统上使用12个线程评估工作流，使用来自三个不同物种（S. cerevisiae， M. musus和H. sapiens）的RNA-seq和Ribo-seq数据。
软件中的多线程使用 python 开发，因为众所周知的原因，不建议使用太高的线程，收益较低。

| | | | | | | | | | | |
|-|-|-|-|-|-|-|-|-|-|-|
||||||Index building| |Preprocessing & Alignment| |Riboparser| |
|species|Dataset|library|sample number|sample size|Elapsed time|Disk usage|Elapsed time|Disk usage|Elapsed time|Disk usage|
|S. cerevisiae|GSE67387|Ribo-seq|6|32 G|38 s|357 M|43 m 21 s|30 G|1 h 26 m 23s|3.6 G|
| | |RNA-seq|6|17 G|38 s|357 M|32 m 52 s|26 G|37 m 42 s|2.8 G|
|M. musculus|GSE114064|Ribo-seq|6|43 G|59 m 8 s|36G|50 m 27 s|31 G|32 m 3 s|7.9 G|
| | |RNA-seq|6|60 G|59 m 8 s|36G|4 h 14 m 45 s|62 G|29 m 30 s|7.6 G|
|H. sapiens|GSE131650|Ribo-seq|6|42G|1 h 55 m 56 s|44 G|2 h 11 m 42 s|29 G|2 h 18 m 57 s|14 G|
| | |RNA-seq|6|54G|1 h 55 m 56 s|44 G|1 h 15 m 15 s|30 G|40 m 35 s|11 G|

RiboParser系统推荐：
为了获得最佳性能，我们建议在基于linux的系统上部署RiboParser（已经在Ubuntu 20.04 LTS/CentOS 7上进行了测试）。

较低配置
- Memory: ≥ 16 GB RAM
- Processor: ≥ 4-core CPU (Intel Xeon E5-2600+ or equivalent)
- Storage: ≥ 512 GB HDD (SATA III)

较高配置
- Memory: ≥ 32 GB RAM
- Processor: ≥ 8-core CPU (AMD EPYC 7B12/Intel i9-10900X)
- Storage: ≥ 512 GB NVMe SSD for rapid I/O and 2 TB HDD (SATA III)


## 10. 贡献

感谢在这个过程中使用的所有开源工具。

感谢Nedialkova DD和Leidel SA提供了高质量的数据集。

欢迎提交问题和代码对我们的项目做出贡献。

更多信息请联系 `rensc0718@163.com`。

## 11. License

GPL License.
