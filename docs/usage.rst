=====
Usage
=====

To use g2gtools in command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Note:** We will assume you installed g2gtools in its own conda virtual environment. First of all, you have to "activate" the virtual environment by doing the following::

    source activate g2gtools

To create a custom genome, we need the following information::

    ${REF}         reference genome in fasta format
    ${STRAIN}      strain name (usually a column name in the vcf file), e.g., CAST_EiJ
    ${VCF_INDELS}  vcf file for indels
    ${VCF_SNPS}    vcf file for snps
    ${GTF}         gene annotation file in gtf format

First of all, we need to create a chain file for mapping bases between two genomes. In this case, between reference and some other strain, like CAST_EiJ::

    $ g2gtools vcf2chain -f ${REF} -i ${VCF_INDELS} -s ${STRAIN} -o ${STRAIN}/REF-to-${STRAIN}.chain

Then we patch snps onto reference genome::

    $ g2gtools patch -i ${REF} -s ${STRAIN} -v ${VCF_SNPS} -o ${STRAIN}/${STRAIN}.patched.fa

We incorporate indels onto the snp-patched genome::

    $ g2gtools transform -i ${STRAIN}/${STRAIN}.patched.fa -c ${STRAIN}/REF-to-${STRAIN}.chain -o ${STRAIN}/${STRAIN}.fa

Create custom gene annotation with respect to the new custom genome. We also create custom annotation database (so we can extract from custom genome) in the following steps::

    $ g2gtools convert -c ${STRAIN}/REF-to-${STRAIN}.chain -i ${GTF} -f gtf -o ${STRAIN}/${STRAIN}.gtf
    $ g2gtools gtf2db -i ${STRAIN}/${STRAIN}.gtf -o ${STRAIN}/${STRAIN}.gtf.db

We can also extract regions of interest from the custom genome. For example,::

    $ g2gtools extract --transcripts -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.transcripts.fa
    $ g2gtools extract --genes -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.genes.fa
    $ g2gtools extract --exons -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.exons.fa

If snp-patched genome is not of interest, we can remove it::

    $ rm ${STRAIN}/${STRAIN}.patched.fa
    $ rm ${STRAIN}/${STRAIN}.patched.fa.fai


To use g2gtools in a project
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All the features are available as a python module. You can simply::

    import g2gtools


