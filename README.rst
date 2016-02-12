
====================
EMPANADA Documentation
====================

EMPANADA is a tool for evidence-based assignment of genes to pathways in metagenomic data,
developed and maintained by the Borenstein group at the University of Washington.

============
Availability
============

EMPANADA is available as a Python module from GitHub or PyPI (see installation instructions below)

=======
License
=======

EMPANADA is distributed under a non-commercial license (see LICENSE).

=========================
Installation Instructions
=========================

Prerequisites for installing:

In order for EMPANADA to run successfully, the following Python modules should be pre-installed on your system:

- Numpy >= 1.6.1 (http://www.numpy.org/)
- Pandas >= 0.14 (http://pandas.pydata.org/)

If you have *pip* installed, you can install these packages by running the following command:

``pip install -U numpy pandas``

**Installing EMPANADA:**

To install EMPANADA, download the package from https://github.com/borenstein-lab/empanada/archive/0.0.1.tar.gz

After downloading EMPANADA, you’ll need to unzip the file. If you’ve downloaded the release version, do this with the following command:

``tar -xzf empanada-0.0.1.tar.gz``

You’ll then change into the new EMPANADA directory as follows:

``cd empanada-0.0.1``

and install using the following command:

``python setup.py install``

ALTERNATIVELY, you can install EMPANADA directly from PyPI by running:

``pip install -U empanada``

============================
Testing the software package
============================

After downloading and installing the software, we recommend testing it by running the following command:

``test_empanada.py``

This will invoke a series of tests. A correct output should end with:

``Ran 1 tests in X.XXXXs``

``OK``

=================================
EMPANADA API via the command line
=================================
The EMPANADA module handles all calculations internally.
EMPANADA offers an interface to the EMPANADA functionality via the command line and the run_empanada script.

Usage:
------

``run_empanada.py -ko KO_ABUNDANCE_FILE -ko2path KO_TO_PATHWAY_FILE [options]``

Required arguments:
-------------------

**-ko KO_ABUNDANCE_FILE**
    Input KO abundance file to aggregate to pathway abundance

**-ko2path KO_TO_PATHWAY_FILE**
    Input file of KO-to-pathway mapping

Optional arguments:
-------------------

**-h, --help**
    show help message and exit

**-o, --output**,
    Output file for resulting pathway abundance (default: out.tab)

**-oc, --output_counts**,
    Output file for number of KOs mapped to each pathway (default: counts.tab)

**-om, --output_mapping**,
    Output the mapping table (either given or generated) to file, works only with pooled mappings (default: mapping.tab)

**-map {naive, by_support, by_sum_abundance, by_avg_abundance}, --mapping_method {naive, by_support, by_sum_abundance, by_avg_abundance}**
    Method to map KOs to Pathway (default: naive)

**-compute {sum}, --compute_method {sum}**
    Method to compute pathway abundance from mapped KOs (default: sum)

**-threshold, --abundance_threshold**
    Abundance threshold to include KOs (default: 0.0)

**-fraction, --fractional_ko_contribution**
    Divide KO contributions such that they sum to 1 for each KO (default: False)

**-remove_ko_with_no_pathway**
    Remove KOs with no pathway from analysis (default: False)

**-remove_ko_with_no_abundance_measurement**
    Remove KOs with no measurements in the abundance table from analysis (default: False)

**-transpose_ko, --transpose_ko_abundance**
    Transpose the ko abundance matrix given (default: False)

**-transpose_output, --transpose_output**
    Transpose the output pathway abundance matrix (default: False)

**-permute_ko_mapping**
    Permute the given KO mapping, i.e., which KO map to which pathways for hypothesis testing (default: False)

**-use_only_non_overlapping_genes**
    If the mapping is by_abundance, compute pathway support by only using non-overlapping genes (default: False)

**-pool_samples_use_median**
    If the mapping is by_abundance, pool samples together using the median KO abundance, and learn the mapping only once (default: False)

**-pool_samples_use_average**
    If the mapping is by_abundance, pool samples together using the average KO abundance, and learn the mapping only once (default: False)

**-leave_one_ko_out_pathway_support**
    If the mapping is by_abundance, compute pathway support for each KO separately by removing it from the computation (default: False)

**-compute_support_with_weighted_double_counting**
    If the mapping is by_abundance, double count KO abundance (weighted by mapping) when computing pathway support (default: False)

**-v, --verbose**
    Increase verbosity of module (default: False)


========
Examples
========

In the *empanada/examples* directory, the file *simulated_ko_relative_abundance.tab* contains simulated KO abundance measurements of 20 samples.
Using this file as input for EMPANADA results in the following files:

- pathway_abundance_empanada.tab

The command used are the following (via command line):

``run_empanada.py -ko examples/simulated_ko_relative_abundance.tab -ko2path data/KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15.tab -o examples/pathway_abundance_empanada.tab -threshold 0 -map by_avg_abundance -fraction -leave_one_ko_out_pathway_support -use_only_non_overlapping_genes``

==================
Citing Information
==================

If you use the EMPANADA software, please cite the following paper:

Functional variability in the human microbiome: More than meets the eye
**Ohad Manor and Elhanan Borenstein.** *In preparation*

