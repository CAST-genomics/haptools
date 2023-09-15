.. _commands-validate:


validate
========

Validate the structure of a ``.hap`` file.

When a ``.hap`` file contains any errors, they will be logged accordingly.

If provided, the SNPs and TRs present in the ``.hap`` file will be confirmed to exist in a ``.pvar`` file.

Usage
~~~~~
.. code-block:: bash

  haptools validate \
  --sort \
  --genotypes PATH \
  --verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
  HAPFILE

Examples
~~~~~~~~
.. code-block:: bash

  haptools validate tests/data/hapfiles/basic.hap

Outputs a message specifying the amount of errors and warnings.

.. code-block::

  [    INFO] Completed HapFile validation with 0 errors and 0 warnings.

All warnings and errors will be logged if there are any.

.. code-block:: bash

  haptools validate tests/data/hapfiles/valhap_with_no_version.hap

.. code-block::

  [ WARNING] No version declaration found. Assuming to use the latest version.
  [    INFO] Completed HapFile validation with 0 errors and 1 warnings.
  [ WARNING] Found several warnings and / or errors in the hapfile

One can use ``--no-sort`` to avoid sorting the file.
This will make it so that all unordered files will get removed, such as out-of-header lines with meta information.

.. code-block:: bash

  haptools validate --no-sort tests/data/hapfiles/valhap_with_out_of_header_metas.hap

Will turn:

.. code-block::

  #   orderH	ancestry	beta
  #	version	0.2.0
  #H	ancestry	s	Local ancestry
  #H	beta	.2f	Effect size in linear model
  #R	beta	.2f	Effect size in linear model
  H	21	26928472	26941960	chr21.q.3365*1	ASW	0.73
  R	21	26938353	26938400	21_26938353_STR	0.45
  H	21	26938989	26941960	chr21.q.3365*10	CEU	0.30
  H	21	26938353	26938989	chr21.q.3365*11	MXL	0.49
  # This should cause an error if the file is sorted
  #V	test_field	s	A field to test with
  V	chr21.q.3365*1	26928472	26928472	21_26928472_C_A	C
  V	chr21.q.3365*1	26938353	26938353	21_26938353_T_C	T
  V	chr21.q.3365*1	26940815	26940815	21_26940815_T_C	C
  V	chr21.q.3365*1	26941960	26941960	21_26941960_A_G	G
  V	chr21.q.3365*10	26938989	26938989	21_26938989_G_A	A
  V	chr21.q.3365*10	26940815	26940815	21_26940815_T_C	T
  V	chr21.q.3365*10	26941960	26941960	21_26941960_A_G	A
  V	chr21.q.3365*11	26938353	26938353	21_26938353_T_C	T
  V	chr21.q.3365*11	26938989	26938989	21_26938989_G_A	A

Into

.. code-block::

  #	orderH	ancestry	beta
  #	version	0.2.0
  #H	ancestry	s	Local ancestry
  #H	beta	.2f	Effect size in linear model
  #R	beta	.2f	Effect size in linear model
  H	21	26928472	26941960	chr21.q.3365*1	ASW	0.73
  R	21	26938353	26938400	21_26938353_STR	0.45
  H	21	26938989	26941960	chr21.q.3365*10	CEU	0.30
  H	21	26938353	26938989	chr21.q.3365*11	MXL	0.49
  V	chr21.q.3365*1	26928472	26928472	21_26928472_C_A	C
  V	chr21.q.3365*1	26938353	26938353	21_26938353_T_C	T
  V	chr21.q.3365*1	26940815	26940815	21_26940815_T_C	C
  V	chr21.q.3365*1	26941960	26941960	21_26941960_A_G	G
  V	chr21.q.3365*10	26938989	26938989	21_26938989_G_A	A
  V	chr21.q.3365*10	26940815	26940815	21_26940815_T_C	T
  V	chr21.q.3365*10	26941960	26941960	21_26941960_A_G	A
  V	chr21.q.3365*11	26938353	26938353	21_26938353_T_C	T
  V	chr21.q.3365*11	26938989	26938989	21_26938989_G_A	A


If the previous example were to be sorted then there would be several errors in the ``.hap`` file.
All sorted files parse the meta information lines first, thus the ``V`` lines would be incomplete.

As mentioned before, one can use the ``--genotypes`` flag to provide a ``.pvar`` file with which to compare the existence of variant IDs.
The following will check if all of the variant IDs in the ``.hap`` appear in the ``.pvar`` file.

.. code-block:: bash

  haptools validate --genotypes tests/data/hapfiles/valhap_test_data.pvar tests/data/hapfiles/valhap_test_data.hap

.. note::

  We accept a PVAR file instead of a VCF in order to avoid reading lots of
  information which is not relevant to the validation process. However, any
  VCF wihtout a FORMAT field is a valid PVAR file. So you can easily create a PVAR file
  using the ``cut`` command or ``plink2 --make-just-pvar``.

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
  :prog: haptools
  :show-nested:
  :commands: validate
