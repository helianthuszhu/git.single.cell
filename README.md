# git.single.cell
# umitools
#! /bin/env bash
# Step 1: get data

# Step 2: Identify correct cell barcodes
umi_tools whitelist --stdin SYL1_S1_L008_R1_001.fastq.gz --bc-pattern='(?P<cell_1>.{16})(?P<umi_1>.{10})T{3}.*' --plot-prefix=expect_whitelist --extract-method=regex --set-cell-number=994 --method=umis --log2stderr > whitelist.txt
                    
# Step 3: Extract barcdoes and UMIs and add to read names
--bc-pattern="(?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT){s<=2}(?P<cell_2>.{8})(?P<umi_1>.{6})T{3}.*"

umi_tools extract --bc-pattern= '(?P<cell_1>.{16})(?P<umi_1>.{10})T{3}.*' \
                  --extract-method=regex \
                  --stdin SYL1_S1_L008_R1_001.fastq.gz \
                  --stdout SYL1_S1_L008_R1_001_extracted.fastq.gz \
                  --read2-in SYL1_S1_L008_R2_001.fastq.gz \
                  --read2-out=SYL1_S1_L008_R2_001_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist=whitelist.txt;
                  			  
zcat hgmm_100_R1_extracted.fastq.gz | head -n4
# Step 4: Map reads
STAR --runThreadN 4 \
       --genomeDir ./STAR.index.GRCh37 \
       --readFilesIn SYL1_S1_L008_R2_001_extracted.fastq.gz \
       --readFilesCommand zcat \
       --outFilterMultimapNmax 1 \
       --outSAMtype BAM SortedByCoordinate
nohup STAR --runThreadN 12 --genomeDir /home/zhuxq/nas/Xiaoqiang/SYL_scRNA/FASTQ/CB8E5ANXX/SYL1.umi.tools/STAR.index.GRCh37 --readFilesIn SYL1_S1_L008_R2_001_extracted.fastq.gz --readFilesCommand zcat --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate

# Step 5: Assign reads to genes
featureCounts -a /home/zhuxq/nas/Xiaoqiang/SYL_scRNA/FASTQ/CB8E5ANXX/SYL1.eoulsan/annotation.GRCh37/Homo_sapiens.GRCh37.85.gtf \
              -o SYL1_gene_assigned \
              -R BAM Aligned.sortedByCoord.out.bam \
              -T 12;            
samtools sort Aligned.sortedByCoord.out.bam.featureCounts.bam -o SYL1.assigned_sorted.bam;
samtools index SYL1.assigned_sorted.bam;
              
# Step 6: Count UMIs per gene per cell
umi_tools count --per-gene --wide-format-cell-counts --gene-tag=XT --assigned-tag=XS --per-cell -I SYL1.assigned_sorted.bam -S SYL1.counts.tsv.gz
            