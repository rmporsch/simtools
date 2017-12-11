Usage
#######

genotypes
------------

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

The package provides a number of mechanisms to simulate different kinds of data.
You can simulate a single phenoytpe on the basis of an existing genotypematrix.::

    import simtools.simtools as st
    heratibiltiy = 0.4
    num_causal_snps = 10
    sims = st.Simtools(genotypematrix)
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
    

