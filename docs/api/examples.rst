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
