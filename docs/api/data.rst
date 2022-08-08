.. _api-data:


data
====

Overview
~~~~~~~~

The ``data`` module is a submodule of ``haptools`` that handles IO for common file types including genotypes, haplotypes, phenotypes, and model covariates.

.. note::
	This page documents common use-cases. To implement more advanced patterns, take a look at our detailed :ref:`API docs <api-haptools-data>`.

Motivation
----------
Using the ``data`` module makes your code agnostic to the type of file being used. For example, the :class:`Genotypes` class provides the same interface for reading VCFs as it does for PLINK2 PGEN files.

This module also helps reduce common boilerplate, since users can easily *extend* the classes in the ``data`` module for their own specific goals.

Data
~~~~
Classes in the ``data`` module inherit from an abstract class called Data, providing some level of standardization across classes within the module. All classes are initialized with the path to the file containing the data and, optionally, a `python Logger <https://docs.python.org/3/howto/logging.html>`_ instance.

The abstract class requires that all classes contain methods for...

1. reading the contents of a file into a ``data`` property of each class
2. iterating over lines of a file without loading all of the data into memory at once

.. code-block:: python

	from haptools import data
	data.Data

genotypes.py
~~~~~~~~~~~~
Overview
--------
This module supports reading and writing files that follow `the VCF <https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format>`_ and `PLINK2 PGEN <https://www.cog-genomics.org/plink/2.0/formats#pgen>`_ file format specifications. We may also offer support for `BGEN <https://www.well.ox.ac.uk/~gav/bgen_format>`_ and `Hail <https://hail.is/docs/0.2/methods/impex.html#native-file-formats>`_ files in the future.

Documentation
-------------
The :ref:`genotypes.py API docs <api-haptools-data-genotypes>` contain example usage of the :class:`Genotypes` class.
See the documentation for the :class:`GenotypesRefAlt` class for an example of extending the :class:`Genotypes` class so that it loads REF and ALT alleles as well.

Classes
-------
Genotypes
+++++++++
Properties
**********
The ``data`` property of a :class:`Genotypes` object is a numpy array representing the genotype matrix. Rows of the array are samples and columns are variants. Each entry in the matrix is a tuple of values -- one for each chromosome. Each value is an integer denoting the index of the allele (0 for REF, 1 for the first ALT allele, 2 for the next ALT allele, etc).

There are two additional properties that contain variant and sample metadata. The ``variants`` property is a numpy structured array and the ``samples`` property is a simple tuple of sample IDs. The ``variants`` structured array has four named columns: "id", "chrom", "pos", and "aaf" (alternate allele frequency).

Reading a file
**************
Extracting genotypes from a VCF file is quite simple:

.. code-block:: python

	genotypes = data.Genotypes.load('tests/data/simple.vcf')
	genotypes.data     # a numpy array of shape n x p x 2
	genotypes.variants # a numpy structured array of shape p x 4
	genotypes.samples  # a tuple of strings of length n

The ``load()`` method initializes an instance of the :class:`Genotypes` class, calls the ``read()`` method, and then performs some standard :ref:`quality-control checks <api-data-genotypes-quality-control>`. You can also call the ``read()`` method manually if you'd like to forego these checks.

.. code-block:: python

	genotypes = data.Genotypes('tests/data/simple.vcf')
	genotypes.read()
	genotypes.data     # a numpy array of shape n x p x 3
	genotypes.variants # a numpy structured array of shape p x 4
	genotypes.samples  # a tuple of strings of length n

	# check that all genotypes are phased and remove the phasing info (in the third dimension)
	genotypes.check_phase()
	genotypes.data     # a numpy array of shape n x p x 2

Both the ``load()`` and ``read()`` methods support ``region``, ``samples``, and ``variants`` parameters that allow you to request a specific region, list of samples, or set of variant IDs to read from the file.

.. code-block:: python

	genotypes = data.Genotypes('tests/data/simple.vcf.gz')
	genotypes.read(
	    region="1:10115-10117",
	    samples=["HG00097", "HG00100"],
	    variants={"1:10117:C:A"},
	)

The ``region`` parameter only works if the file is indexed, since in that case, the ``read()`` method can take advantage of the indexing to parse the file a bit faster.

Iterating over a file
*********************
If you're worried that the contents of the VCF file might be large, you may opt to parse the file line-by-line instead of loading it all into memory at once.

In cases like these, you can use the ``__iter__()`` method in a for-loop:

.. code-block:: python

	genotypes = data.Genotypes('tests/data/simple.vcf')
	for line in genotypes:
	    print(line)

You'll have to call ``__iter()__`` manually if you want to specify any function parameters:

.. code-block:: python

	genotypes = data.Genotypes('tests/data/simple.vcf.gz')
	for line in genotypes.__iter__(region="1:10115-10117", samples=["HG00097", "HG00100"]):
	    print(line)

.. _api-data-genotypes-quality-control:

Quality control
***************
There are several quality-control checks performed by default (in the ``load()`` method). You can call these methods yourself, if you'd like:

1. ``check_missing()`` - raises an error if any samples are missing genotypes
2. ``check_biallelic()`` - raises an error if any variants have more than one ALT allele
3. ``check_phase()`` - raises an error if any genotypes are unphased

Subsetting
**********
You can index into a loaded :class:`Genotypes` instance using the ``subset()`` function. This works similiar to numpy indexing with the added benefit that you can specify a subset of variants and/or samples by their IDs instead of just their indices.

.. code-block:: python

	genotypes = data.Genotypes.load('tests/data/simple.vcf')
	gts_subset = genotypes.subset(samples=("HG00100", "HG00101"), variants=("1:10114:T:C", '1:10116:A:G'))
	gts_subset # a new Genotypes instance containing only the specified samples and variants

By default, the ``subset()`` method returns a new :class:`Genotypes` instance. The samples and variants in the new instance will be in the order specified.

GenotypesRefAlt
+++++++++++++++
The :class:`Genotypes` class can be easily *extended* (sub-classed) to load extra fields into the ``variants`` structured array. The :class:`GenotypesRefAlt` class is an example of this where I extended the :class:`Genotypes` class to add REF and ALT fields from the VCF to the columns of the structured array. So the ``variants`` array will have named columns: "id", "chrom", "pos", "aaf", "ref", and "alt".

All of the other methods in the :class:`Genotypes` class are inherited, but the :class:`GenotypesRefAlt` class implements an additional method ``write()`` for dumping the contents of the class to the provided file.

.. code-block:: python

	genotypes = data.GenotypesRefAlt.load('tests/data/simple.vcf')
	# make the first sample homozygous for the alt allele of the fourth variant
	genotypes.data[0, 3] = (1, 1)
	genotypes.write()

.. _api-data-genotypesplink:

GenotypesPLINK
++++++++++++++
The :class:`GenotypesPLINK` class offers experimental support for reading and writing PLINK2 PGEN, PVAR, and PSAM files. We are able to read genotypes from a PLINK2 PGEN files in a fraction of the time of VCFs. Reading from VCFs is :math:`O(n*p)`, while reading from PGEN files is approximately :math:`O(1)`.

.. figure:: https://drive.google.com/uc?export=view&id=1_JARKJQ0LX-DzL0XsHW1aiQgLCOJ1ZvC

	The time required to load various genotype file formats.

.. warning::
	Use of this class is not officially supported yet because it relies upon the as-yet-unpublished ``pgenlib`` python library. See `issue #16 <https://github.com/gymrek-lab/haptools/pull/16>`_ for current progress on this challenge. In the meantime, you must install the library manually from Github via ``pip``.

	.. code-block:: bash

		pip install git+https://github.com/chrchang/plink-ng.git#subdirectory=2.0/Python

The :class:`GenotypesPLINK` class inherits from the :class:`GenotypesRefAlt` class, so it has all the same methods and properties. Loading genotypes is the exact same, for example.

.. code-block:: python

	genotypes = data.GenotypesPLINK.load('tests/data/simple.pgen')
	genotypes.data     # a numpy array of shape n x p x 2
	genotypes.variants # a numpy structured array of shape p x 6
	genotypes.samples  # a tuple of strings of length n

In addition to the ``read()`` and ``load()`` methods, the :class:`GenotypesPLINK` class also has methods for reading (or writing) PVAR or PSAM files separately, without having to read (or write) the PGEN file as well.

.. code-block:: python

	genotypes = data.GenotypesPLINK('tests/data/simple.pgen')

	genotypes.read_variants()
	genotypes.variants # a numpy structured array of shape p x 6

	genotypes.read_samples()
	genotypes.samples  # a tuple of strings of length n

	genotypes.data     # simply None

Limiting memory usage
*********************
Unfortunately, reading from PGEN files can require a lot of memory, at least initially. (Once the genotypes have been loaded, they are converted down to a lower-memory form.) To determine whether you may be having memory issues, you can place the module in "verbose mode" by providing a `python Logger <https://docs.python.org/3/howto/logging.html>`_ object at the "DEBUG" level when initializing the :class:`GenotypesPLINK` class.

.. code-block:: python

	import logging
	log = logging.getLogger("debug_plink_mem")
	logging.basicConfig(format="[%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)", level="DEBUG")

	genotypes = data.GenotypesPLINK('tests/data/simple.pgen', log=log)
	genotypes.read()

If you find yourself running out of memory when trying to load a PGEN file, you may want to try loading the genotypes in chunks. You can specify the number of variants to read (and write) together at once via the ``chunk_size`` parameter. This parameter is only available for the :class:`GenotypesPLINK` class.

A large ``chunk_size`` is more likely to result in memory over-use while a small ``chunk_size`` will increase the time it takes to read the file. If the ``chunk_size`` is not specified, all of the genotypes will be loaded together in a single chunk.

.. code-block:: python

	genotypes = data.GenotypesPLINK('tests/data/simple.pgen', chunk_size=500)
	genotypes.read()

haplotypes.py
~~~~~~~~~~~~~
Overview
--------
This module supports reading and writing files that follow the **.hap** file format specification.

Lines from the file are parsed into instances of the :class:`Haplotype` and :class:`Variant` classes. These classes can be *extended* (sub-classed) to support "extra" fields appended to the ends of each line.

Documentation
-------------

1. The **.hap** :ref:`format specification <formats-haplotypes>`
2. The :ref:`haplotypes.py API docs <api-haptools-data-haplotypes>` contain example usage of the :class:`Haplotypes` class and examples of sub-classing the :class:`Haplotype` and :class:`Variant` classes

Classes
-------
Haplotypes
++++++++++
Reading a file
**************
Parsing a basic **.hap** file without any extra fields is as simple as it gets:

.. code-block:: python

	haplotypes = data.Haplotypes.load('tests/data/basic.hap')
	haplotypes.data # returns a dictionary of Haplotype objects

The ``load()`` method initializes an instance of the :class:`Haplotypes` class and calls the ``read()`` method, but if the **.hap** file contains extra fields, you'll need to call the ``read()`` method manually. You'll also need to create :class:`Haplotype` and :class:`Variant` subclasses that support the extra fields and then specify the names of the classes when you initialize the :class:`Haplotypes` object:

.. code-block:: python

	haplotypes = data.Haplotypes('tests/data/basic.hap', Haplotype, Variant)
	haplotypes.read()
	haplotypes.data # returns a dictionary of Haplotype objects

Both the ``load()`` and ``read()`` methods support `region` and `haplotypes` parameters that allow you to request a specific region or set of haplotype IDs to read from the file.

.. code-block:: python

	haplotypes = data.Haplotypes('tests/data/basic.hap.gz', Haplotype, Variant)
	haplotypes.read(region='chr21:26928472-26941960', haplotypes=["chr21.q.3365*10"])

The file must be indexed if you wish to use these parameters, since in that case, the ``read()`` method can take advantage of the indexing to parse the file a bit faster. Otherwise, if the file isn't indexed, the ``read()`` method will assume the file could be unsorted and simply reads each line one-by-one. Although I haven't tested it yet, streams like stdin should be supported by this case.

Iterating over a file
*********************
If you're worried that the contents of the **.hap** file will be large, you may opt to parse the file line-by-line instead of loading it all into memory at once.

In cases like these, you can use the ``__iter__()`` method in a for-loop:

.. code-block:: python

	haplotypes = data.Haplotypes('tests/data/basic.hap')
	for line in haplotypes:
	    print(line)

You'll have to call ``__iter()__`` manually if you want to specify any function parameters:

.. code-block:: python

	haplotypes = data.Haplotypes('tests/data/basic.hap')
	for line in haplotypes.__iter__(region='21:26928472-26941960', haplotypes={"chr21.q.3365*1"}):
	    print(line)

Writing a file
**************
To write to a **.hap** file, you must first initialize a :class:`Haplotypes` object and then fill out the data property:

.. code-block:: python

	haplotypes = data.Haplotypes('tests/data/example-write.hap')
	haplotypes.data = {}
	haplotypes.data['H1'] = Haplotype(chrom='chr1', start=0, end=10, id='H1')
	haplotypes.data['H1'].variants = [Variant(start=0, end=1, id='rs123', allele='A')]
	haplotypes.write()

Obtaining haplotype "genotypes"
*******************************
Using the ``transform()`` function, you can obtain a full instance of the :class:`GenotypesRefAlt` class where haplotypes from a :class:`Haplotypes` object are encoded as the variants in the genotype matrix.

.. code-block:: python

	haplotypes = data.Haplotypes.load('tests/data/example.hap.gz')
	genotypes = data.GenotypesRefAlt.load('tests/data/example.vcf.gz')
	hap_gts = haplotypes.transform(genotypes)
	hap_gts   # a GenotypesRefAlt instance where haplotypes are variants

Haplotype
+++++++++
The :class:`Haplotype` class stores haplotype lines from the **.hap** file. Each property in the object is a field in the line. A separate ``variants`` property stores a tuple of :class:`Variant` objects belonging to this haplotype.

The :class:`Haplotypes` class will initialize :class:`Haplotype` objects in its ``read()`` and ``__iter__()`` methods. It uses a few methods within the :class:`Haplotype` class for this:

1. ``from_hap_spec()`` - this static method initializes a Haplotype object from a line in the **.hap** file.
2. ``to_hap_spec()`` - this method converts a Haplotype object into a line in the **.hap** file

To read "extra" fields from a **.hap** file, one need only *extend* (sub-class) the base :class:`Haplotype` class and add the extra properties that you want to load. For example, let's add an extra field called "ancestry" that is encoded as a string.

.. code-block:: python

    from dataclasses import dataclass, field
    from haptools.data import Haplotype, Extra

    @dataclass
    class CustomHaplotype(Haplotype):
        score: float
        _extras: tuple = field(
            repr=False,
            init=False,
            default=(
                Extra("ancestry", "s", "Local ancestry"),
            ),
        )

    haps = Haplotypes("file.hap", haplotype=CustomHaplotype)
    haps.read()
    haps.write()

Variant
+++++++
The :class:`Variant` class stores variant lines from the **.hap** file. Each property in the object is a field in the line.

The :class:`Haplotypes` class will initialize :class:`Variant` objects in its ``read()`` and ``__iter__()`` methods. It uses a few methods within the :class:`Variant` class for this:

1. ``from_hap_spec()`` - this static method initializes a :class:`Variant` object from a line in the **.hap** file.
2. ``to_hap_spec()`` - this method converts a :class:`Variant` object into a line in the **.hap** file

To read "extra" fields from a **.hap** file, one need only *extend* (sub-class) the base :class:`Variant` class and add the extra properties that you want to load. For example, let's add an extra field called "score" that is encoded as a float with a precision of three decimal places.

.. code-block:: python

    from dataclasses import dataclass, field
    from haptools.data import Haplotype, Extra

    @dataclass
    class CustomVariant(Variant):
        score: float
        _extras: tuple = field(
            repr=False,
            init=False,
            default=(
                Extra("score", ".3f", "Importance of inclusion"),
            ),
        )

    haps = Haplotypes("file.hap", variant=CustomVariant)
    haps.read()
    haps.write()

phenotypes.py
~~~~~~~~~~~~~
Overview
--------
This module supports reading and writing PLINK2-style phenotype files.

Documentation
-------------

1. The **.pheno** `phenotype format specification <https://www.cog-genomics.org/plink/2.0/input#pheno>`_
2. The :ref:`phenotypes.py API docs <api-haptools-data-phenotypes>` contain example usage of the :class:`Phenotypes` class

Classes
-------
Phenotypes
++++++++++
Reading a file
**************
Loading a **.pheno** file is easy:

.. code-block:: python

	phenotypes = data.Phenotypes.load('tests/data/simple.pheno')
	phenotypes.data # returns a np array of shape p x k

The ``load()`` method initializes an instance of the :class:`Phenotypes` class and calls the ``read()`` method as well as the ``standardize()`` method. To forego the standardization, you'll need to call the ``read()`` method manually.

.. code-block:: python

	phenotypes = data.Phenotypes('tests/data/simple.pheno')
	phenotypes.read()
	phenotypes.data # returns a np array of shape p x k

Both the ``load()`` and ``read()`` methods support the `samples` parameter that allows you to request a specific set of sample IDs to read from the file.

.. code-block:: python

	phenotypes = data.Phenotypes('tests/data/simple.pheno')
	phenotypes.read(samples=["HG00097", "HG00099"])

Iterating over a file
*********************
If you're worried that the contents of the **.pheno** file will be large, you may opt to parse the file line-by-line instead of loading it all into memory at once.

In cases like these, you can use the ``__iter__()`` method in a for-loop:

.. code-block:: python

	phenotypes = data.Phenotypes('tests/data/simple.pheno')
	for line in phenotypes:
	    print(line)

You'll have to call ``__iter()__`` manually if you want to specify any function parameters:

.. code-block:: python

	phenotypes = data.Phenotypes('tests/data/simple.pheno')
	for line in phenotypes.__iter__(samples=["HG00097", "HG00099"]):
	    print(line)

Writing a file
**************
To write to a **.pheno** file, you must first initialize a :class:`Phenotypes` object and then fill out the necessary properties:

.. code-block:: python

	phenotypes = data.Phenotypes('tests/data/example-write.pheno')
	phenotypes.data = np.array([[1, 0.2], [1, 0.5], [0, 0.9]], dtype='float64')
	phenotypes.samples = ("HG00097", "HG00099", "HG00100")
	phenotypes.names = ("height", "bmi")
	phenotypes.write()

covariates.py
~~~~~~~~~~~~~
Overview
--------
This module supports reading and writing PLINK2-style covariate files.

Documentation
-------------

1. The **.covar** `covariate format specification <https://www.cog-genomics.org/plink/2.0/input#covar>`_
2. The :ref:`covariates.py API docs <api-haptools-data-covariates>` contain example usage of the :class:`Covariates` class

Classes
-------
Covariates
++++++++++
The :class:`Covariates` class is simply a sub-class of the :class:`Phenotypes` class. It has all of the same methods and properties. There are no major differences between the two classes, except between the file extensions that they use.
