Usage
#######

The main aim of Simtools is to provide an easy to use platform to simulate and manipulate
genetic and phenotype data.

The package was mainly written as a collection of tools for myself to explore various scenarios.
It provides facilities to simulate simple categorical and continuous phenotypes from genetic data, 
as well as more complicated multiple dependent phenotypes.

In addition, the package has a connection to plink to perform various different tasks, such as allelic score computation,
clumping and pruning.
I also provide a vcf reader for rare variant analysis.

Here I will briefly outline the usage of these classes.

genotypes
------------

The genotype module provides two classes to read and sample genotypes.

For plink file you simply can::

    import simtools.genotypes as gp
    n = 100; p=1000
    seedfile = 'plink_stem'
    plink = gp.ReadPlink(seedfile)
    genotypematrix = plink.sample(n, p) 

This  will randomly chose n subjects and p SNPs from a given plink file.

Similar you can also you a vcf file::

    import simtools.genotypes as gp
    n = 100; p=1000
    vcf = gp.ReadVCF('example.vcf.gz')
    vcf.read_vcf()
    genotypematrix = vcf.sample(n, p)

simtools
--------

Simtools is the work horse of the package.
It provides two main functions to simulate phenoytpes.

You can simulate a single phenoytpe on the basis of an existing genotypematrix.::

    import simtools.simtools as st
    heratibiltiy = 0.4
    num_causal_snps = 10
    plink_file = 'plink_stem'
    sims = st.Simtools(plink_file)
    # continuous phenotype
    pheno = sims.simple_phenotype(num_causal_snps, heratibiltiy)
    # binary phenotype
    n_cases = 100
    n_controls = 100
    liability_threshold = 0.1
    pheno = sims.simple_phenotype(num_causal_snps,
        heratibiltiy, (n_cases, n_controls, liability_threshold))

In addition to a simple phenotype the package is also able to simulate multiple inter-related phenotypes.
For example, the module allows you to simulate phenotype A, which is caused by genetic factors as well as phenotype B and C.
To simulate 3 different interrelated phenotypes one can simply::

    import simtools.simtools as st
    import numpy as np
    heratibiltiy = 0.4
    num_causal_snps = [10, 10, 10] #for each phenotype
    # adjacency matrix
    lamb = np.eye(3)
    lamb[1,0] = 0.1
    lamb[2,0] = 0.1
    # genetc effect matrix
    B = np.eye(3); B = B * heratibiltiy5
    
    # simulate phenoytpes
    phenos = sims.multi_phenotype(lamb, B, num_causal_snps)


tools
-------

The tools module allows you to perform common tasks.
This includes:

- GWAS (single or multicore)
- Computation of the genomic inflation factor
- QQ-Plots
- Randomly chose adjacency matrices
- compute GWAS, do clumping and pruning, as well as calculate allelic scores with Plink
