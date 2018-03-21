####################################
Genetic Simulation Tools (GST)
####################################

When developing new methods it is often easier to start with some simulated data.
Here I provide some tools to simulate various phenotypes from imported or artificial genetic data.

Overall the library is separated into 3 components

    - genotypes: facilities to import genetic data
    - simtools: simulation of phenotypes
    - tools: handy tools


***************************
Loading Genotypes
***************************

The genotype module provides two classes to read and sample genotypes for both plink and vcf files.

For plink file you simply can

.. code-block:: python

    import simtools.genotypes as gp
    n = 100; p=1000
    seedfile = 'plink_stem'
    plink = gp.ReadPlink(seedfile)
    genotypematrix = plink.sample(n, p)

This  will randomly chose n subjects and p SNPs from a given plink file.
Furthermore, `sample` also allows to write the samples subjects and variants to disk in the current working directory.

You can also read specific variants and subjects into a genotype-matrix.

.. code-block:: python

    import simtools.genotypes as gp
    import numpy as np
    seedfile = 'plink_stem'
    plink = gp.ReadPlink(seedfile)
    marker = np.random.choice(plink.bim.index.values, 5)
    subjects = np.random.choice(plink.bim.index.values, 5)
    genotypematrix = plink.read_bed(marker, subjects)

If `read_plink` is not provided with any marker or subjects it will attempted to read the entire bed file into memory.

Similar you can also you a vcf file

.. code-block:: python

    import simtools.genotypes as gp
    n = 100; p=1000
    vcf = gp.ReadVCF('example.vcf.gz')
    genotypematrix = vcf.sample(n, p)

Similar to the plink class, vcf also has a `read_vcf` method to obtain the genotype-matrix of specific individuals and
variants.

Furthermore, the vcf class also has a simple method to obtain the MAF of all variants in a vcf file

.. code-block:: python

    import simtools.genotypes as gp
    output = 'path.to.output'
    vcf = gp.ReadVCF('example.vcf.gz')
    vcf.get_allele_frequency(output)

***************************
Phenotype Simulations
***************************

One of the main purposes of this library is to provide a simulation framework to simulate phenotypes.
To start off one first needs a seed population in form of a plink file (support for vcf will be added).

.. code-block:: python

    import simtools.simtools as st
    import numpy as np
    plink_file = 'plink_stem'
    sims = st.Simtools(plink_file)
    sims.chunk_size = 100 # default

The chunk size sets how many variants `simtools` should process at a time.
This might be useful if your seed population is large and cannot be loaded into memory.
In addition, you can also select specific subjects within the seed population when initializing the class.

Single Phenotype Simulation
==============================

Multiple Phenotype Simulation
================================

Use
-----

It is relatively easy to simulate multiple related phenotypes from existing genetic data.

.. code-block:: python

    import simtools.simtools as st
    import numpy as np
    plink_file = 'plink_stem'
    sims = st.Simtools(plink_file)
    heratibiltiy = 0.4
    num_causal_snps = [10, 10, 10] #for each phenotype
    # adjacency matrix
    lamb = np.eye(3)
    lamb[1,0] = 0.1
    lamb[2,0] = 0.1
    # genetc effect matrix
    b = np.eye(3); b = b * heratibiltiy

    # simulate phenoytpes
    phenos = sims.multi_phenotype(lamb, b, num_causal_snps)

The resulting matrix `phenos` is a matrix containing the phenotypes of all subjects.
Please be aware that currently the `multi_phenotype` method does not produce binary traits.

Method
-------

It is often desirable to simulate multiple related phenotypes.
Lets say we have three phenotypes :math:`A, B, C` and we would like to simulate inter-causal relationship among those.
Then we can express the relationship between :math:`A, B, C` in matrix form (also called a adjacency matrix) as

.. math::

    \Lambda =
    \begin{bmatrix}
    0 & \lambda_{A,B} & \lambda_{A,C} \\
    \lambda_{B,A} & 0 & \lambda_{B,C} \\
    \lambda_{C,A} & \lambda_{C,B} & 0
    \end{bmatrix}

in which the effect of A on B is indicated by :math:`\lambda_{B,A}`, the effect of B on A as :math:`\lambda_{A,B}`, and so on.
The phenotype for subject i for t phenotypes is then

.. math::

    y_i &= \Lambda y_i + \beta X_i  + e_i \\
     &= {(I - \Lambda)}^{-1} [\beta X_i+ \varepsilon_i]

where

    - :math:`X_i` represents the vector of variants for a subject i
    - :math:`\beta` the corresponding effect sizes of size :math:`t\times p`
    - :math:`e_i` is the vector of random effects
    - :math:`y_i` is a vector of the phenotypes for subjects i
    - :math:`\Lambda` is a :math:`t\times t` matrix of structural coefficients (see above);

Hence it is easy to see that this approach flexible and allows the simulation of :math:`t` number of arbitrary phenotypes.
However, in the current form it is unable to specify specific genetic and environmental effects.
In order to do so one can generalize the previous equation to all subjects

.. math::

    y &= {(I - \Lambda)}^{-1} [BG + DE] \\
    &= {(I - \Lambda)}^{-1}BG + {(I - \Lambda)}^{-1}DE

in which

    - :math:`G = (X\beta')'` is the genetic for the genotype matrix :math:`X_{n\times p}` and the effect sizes :math:`\beta_{t\times p}`
    - :math:`y` the :math:`t \times n` matrix of phenotypes
    - :math:`E` be the :math:`t \times n` matrix of errors (:math:`E = (\varepsilon_1, \ldots, \varepsilon_t) \quad with \quad \varepsilon \sim N(0,1)`)
    - :math:`B` and :math:`D` the diagonal scalar matrices of size :math:`t \times t` for the genetic risk and error term respectively.

In order to constrain the variance of each phenotype to :math:`1`, let :math:`C = {(I - \Lambda)}^{-1} B` and :math:`F = {(I - \Lambda)}^{-1} D` then we can constrain

.. math::

    \Sigma_y = C\times \Sigma_G \times C^{T} + F\times \Sigma_E \times F^{T} =
    \begin{bmatrix}
    1 & \lambda_{1,2} & \lambda_{1,3} \\
    \lambda_{2,1} & 1 & \lambda_{2,3} \\
    \lambda_{3,1} & \lambda_{3,2} & 1
    \end{bmatrix}

in which :math:`\Sigma_y`, :math:`\Sigma_G` and :math:`\Sigma_E` are the variance-covariance matrices of y, G and E respectively.

Since we constrain the variance of all traits to be :math:`1` the sum of the causal and genetic effect cannot be greater than :math:`1`.
Hence the following constrain has to be respected when constructing both :math:`\Lambda` and :math:`B`

.. math::

    (\Lambda + B) \boldsymbol{1} \leq \boldsymbol{1}



Contents:

.. toctree::
   :maxdepth: 3

   usage
   references
