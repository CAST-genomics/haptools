.. _api-examples:


examples
========

Converting a ``.blocks.det`` file into a ``.hap`` file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can use the :ref:`data API <api-data>` to easily convert `a PLINK 1.9 .blocks.det file <https://www.cog-genomics.org/plink/1.9/formats#blocks>`_ into a ``.hap`` file.

As an example, let's say we would like to convert the following ``.blocks.det`` file.

.. code-block::

  CHR          BP1          BP2           KB  NSNPS SNPS
   1      2313888      2331789       17.902      3 rs7527871|rs2840528|rs7545940
   1      2462779      2482556       19.778      2 rs2296442|rs2246732
   1      2867411      2869431        2.021      2 rs10752728|rs897635
   1      2974991      2979823        4.833      3 rs10489588|rs9661525|rs2993510

.. code-block:: python

    from haptools import data

    # load the genotypes file
    # you can use either a VCF or PGEN file
    gt = data.GenotypesRefAlt.load("input.vcf.gz")
    gt = data.GenotypesPGEN.load("input.pgen")

    # load the haplotypes
    hp = data.Haplotypes("output.hap")
    hp.data = {}

    # iterate through lines of the .blocks.det file
    with open("input.blocks.det") as blocks_file:
        for idx, line in enumerate(blocks_file.readlines()):
            # initialize variables and parse line from the blocks file
            hap_id = f"H{idx}"
            chrom, bp1, bp2, kb, nsnps, snps = line.split("\t")

            # create a haplotype line in the .hap file
            hp.data[hap_id] = data.Haplotype(chrom=chrom, start=bp1, end=bp2, id=hap_id)

            # fetch alleles from the genotypes file
            snp_gts = gt.subset(variants=tuple(snps.split("|")))

            # create variant lines for each haplotype
            # Note that the .blocks.det file doesn't specify an allele, so
            # we simply choose the REF allele for this example
            hp.data[hap_id].variants = tuple(
                data.Variant(start=v["pos"], end=v["pos"]+len(v["ref"]), id=v["id"], allele=v["ref"])
                for v in snp_gts.variants
            )

    hp.write()

.. _api-examples-snps2hap:

Creating a ``.hap`` file of SNPs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :ref:`simphenotype <commands-simphenotype>` command requires a ``.hap`` file containing haplotypes, but what if you want to give it SNPs, instead?

Well, you can encode each SNP as a haplotype containing only a single allele. For example, let's say you have two SNPs, rs429358 and rs7412, with ALT alleles C and T respectively. Then your ``.hap`` file might look something like this.

.. code-block::

    #	orderH	beta
    #	version	0.1.0
    #H	beta	.2f	Effect size in linear model
    H	19	45411941	45411942	rs429358	0.73
    H	19	45412079	45412080	rs7412	0.30
    V	rs429358	45411941	45411942	rs429358	C
    V	rs7412	45412079	45412080	rs7412	T

You can easily use the :ref:`data API <api-data>` and the :ref:`simphenotype API <api-haptools-sim_phenotype>` to create such a file.

.. code-block:: python

    from haptools import data
    from haptools.sim_phenotype import Haplotype

    # which variants do we want to write to the haplotype file?
    variants = {"rs429358", "rs7412"}

    # load the genotypes file
    # you can use either a VCF or PGEN file
    gt = data.GenotypesRefAlt("tests/data/apoe.vcf.gz")
    gt.read(variants=variants)
    # the advantage of using a PGEN file is that you can use read_variants() to load
    # the variants quickly w/o having to load the genotypes, too
    gt = data.GenotypesPGEN("tests/data/apoe.pgen")
    gt.read_variants(variants=variants)

    # initialize an empty haplotype file
    hp = data.Haplotypes("output.hap", haplotype=Haplotype)
    hp.data = {}

    for variant in gt.variants:
        ID, chrom, pos, alt = variant[["id", "chrom", "pos", "alt"]]
        end = pos + len(alt)

        # create a haplotype line in the .hap file
        # you should fill out "beta" with your own value
        hp.data[ID] = Haplotype(chrom=chrom, start=pos, end=end, id=ID, beta=0.5)

        # create variant lines for each haplotype
        hp.data[ID].variants = (data.Variant(start=pos, end=end, id=ID, allele=alt),)

    hp.write()
