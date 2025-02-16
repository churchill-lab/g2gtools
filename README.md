# g2gtools

g2gtools is a new suite of tools that creates personal diploid genomes. It allows us to easily extract regions from personal genomes so we can create individualized alignment indexes for next-generation sequencing reads. g2gtools can liftover alignments on personal genome coordinates back to that of reference so we can compare alignments from among samples in a population. Unlike other liftover tools, g2gtools does not throw away alignments that land on indel regions.

* Free software: MIT License
* Documentation for Release Ver. 1.0.XX: BELOW
* Documentation for Release Ver. 0.2.XX: http://churchill-lab.github.io/g2gtools/.
* Documentation for Release Ver. 0.1.XX: https://g2gtools.readthedocs.org.

Development Lead
----------------

* Algorithm design and software engineering: **Kwangbom "KB" Choi, Ph. D.**, The Jackson Laboratory <kb.choi@jax.org>
* Software engineering: **Matthew J. Vincent**, The Jackson Laboratory <matt.vincent@jax.org>

Contributors
------------

* Narayanan Raghupathy, The Jackson Laboratory <Narayanan.Raghupathy@jax.org>
* Anuj Srivastava, The Jackson Laboratory <Anuj.Srivastava@jax.org>
* Mike Lloyd, The Jackson Laboratory <Mike.Lloyd@jax.org>


## Installation

*Note: To avoid conflicts among dependencies, we highly recommend using a [Python virtual environment](https://realpython.com/python-virtual-environments-a-primer/).*

g2gtools requires Python 3+ to run.  Install g2gtools and all its dependencies from the command line:

```
pip install git+https://github.com/churchill-lab/g2gtools
```

The **g2gtools** script should now be installed and you should be able to run g2gtools from the command line. 

## Usage

To demonstrate g2gtools, here is a typical workflow for creating diploid genome, exome, and transcriptome using a human example: NA19670 from 1000 Genomes (phase 3).


##### Step 1: vcf2vci -  Create VCI file from VCF file(s)


The first step is to create a VCI (Variant Call Information) file from a VCF (Variant Call Format) File. 

```
g2gtools vcf2vci --vci NA19670.vci --strain NA19670 --diploid \
--vcf ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz \
--vcf ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz \
--vcf ALL.chrMT.phase3_callmom.20130502.genotypes.vcf.gz
```

##### Step 2: patch - Patch SNPs onto the reference sequence
#

The next step is to replace the bases in the input (reference) sequence using SNPs in the 
sample-specific VCI file. This command supports patching of local subsequence on the fly.

```
g2gtools patch --fasta hs37d5.fa \
               --vci NA19670.vci.gz \
               --output NA19670.patched.fa
```

##### Step 3: transform - Incorporate indels onto the input sequence
#

Next, indels specified in a sample's VCI file are incorporated into the input (reference) sequence. This
command also supports the transformation of local subsequence on the fly.

```
g2gtools transform --fasta NA19670.patched.fa \
                   --vci NA19670.vci.gz \
                   --out NA19670.fa
```

##### Step 3: convert - Convert coordinates of BAM|SAM|GTF|GFF|BED files
#

Next, convert (lift over) coordinates from reference to custom personal genome. This function supports BAM/SAM, GTF, and BED formats. We implemented g2gtools since every available liftover tools discards alignments that land on indel locations in BAM/SAM files with no justifiable reason. We wanted to avoid removing such alignments since they are possibly enriched with allele-specific information.

```
g2gtools convert --in Homo_sapiens.GRCh37.75.gtf \
                 --vci NA19670.vci.gz \
                 --out NA19670.gtf
```

##### Step 4: gtf2db - Convert a GTF file to a G2G DB file
#

In this step we create a G2G Database File using gene annotation (GTF format). The G2G Database File is used for 'extract'ing annotated elements -- generally exons, transcripts, and gene regions -- from the input genomes.

```
g2gtools gtf2db --gtf NA19670.gtf \
                --db NA19670.db
```

##### Step 5: extract - Extract subsequence from a fasta file given a region
#

Next we extract genomic regions from the input genome. This is useful for creating FASTA file of exons, transcripts, and genes when used with G2G Database file. We can also specify putative peak regions for ChIP-Seq or ATAC-Seq in a BED file to create a "peakome".
                    
```
g2gtools extract --fasta NA19670.fa --db NA19670.db --exons --out NA19670.exons.fa
g2gtools extract --fasta NA19670.fa --db NA19670.db --transcripts --out NA19670.transcripts.fa
g2gtools extract --fasta NA19670.fa --db NA19670.db --genes --out NA19670.genes.fa
```







