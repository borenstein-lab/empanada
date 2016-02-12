====================
EMPANADA Documentation
====================

EMPANADA is a marker genes based framework for metagenomic normalization and accurate profiling of gene abundances in the microbiome,
developed and maintained by the Borenstein group at the University of Washington.

============
Availability
============

EMPANADA is available through the following sources:

- As a Python module from GitHub or PyPI (see installation instructions below)
- As an online tool at: http://elbo.gs.washington.edu/software_empanada.html.

=======
License
=======

EMPANADA is distributed under a BSD license and can be readily incorporated into custom analysis tools.

=========================
Installation Instructions
=========================

Prerequisites for installing:

In order for EMPANADA to run successfully, the following Python modules should be pre-installed on your system:

- Numpy >= 1.6.1 (http://www.numpy.org/)
- Scipy >= 0.9 (http://www.scipy.org/)
- Scikit-learn >= 0.15.2 (http://scikit-learn.org/stable/)
- Pandas >= 0.14 (http://pandas.pydata.org/)

If you have *pip* installed, you can install these packages by running the following command:

``pip install -U numpy scipy scikit-learn pandas``

**Installing EMPANADA:**

To install EMPANADA, download the package from https://github.com/omanor/empanada/archive/1.0.tar.gz

After downloading EMPANADA, you’ll need to unzip the file. If you’ve downloaded the release version, do this with the following command:

``tar -xzf empanada-1.0.tar.gz``

You’ll then change into the new EMPANADA directory as follows:

``cd empanada-1.0``

and install using the following command:

``python setup.py install``

ALTERNATIVELY, you can install EMPANADA directly from PyPI by running:

``pip install -U empanada``

Note for windows users: Under some windows installations, Scipy may fail when importing the Stats module. Workarounds may be found online, such
as `here <https://code.google.com/p/pythonxy/issues/detail?id=745>`_.

============================
Testing the software package
============================

After downloading and installing the software, we recommend testing it by running the following command:

``test_empanada.py``

This will invoke a series of tests. A correct output should end with:

``Ran 3 tests in X.XXXXs``

``OK``

===============================
EMPANADA API via the command line
===============================
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
In the *empanada/examples* directory, the file *simulated_ko_relative_abundance.tab* contains simulated KO abundance measurements of 20 samples described in the
EMPANADA manuscript. Using this file as input for EMPANADA results in the following files:

- simulated_ko_empanada_Normalized.tab (only normalization)
- simulated_ko_empanada_Normalized_Corrected_use_generic.tab (normalize and correct using the generic model learned from HMP)
- simulated_ko_empanada_Normalized_Corrected_learn_model.tab (normalize and correct learning a new model for each sample)

The commands used were the following (via command line):

``run_empanada.py empanada/examples/simulated_ko_relative_abundance.tab -n -perf -v -o empanada/examples/simulated_ko_empanada_Normalized.tab``

``run_empanada.py empanada/examples/simulated_ko_relative_abundance.tab -n -c use_generic -perf -v -o empanada/examples/simulated_ko_empanada_Normalized_Corrected_use_generic.tab``

``run_empanada.py empanada/examples/simulated_ko_relative_abundance.tab -n -c learn_model -perf -v -o empanada/examples/simulated_ko_empanada_Normalized_Corrected_learn_model.tab``

==================
Citing Information
==================

If you use the EMPANADA software, please cite the following paper:

EMPANADA: A marker genes based framework for metagenomic normalization and accurate profiling of gene abundances in the microbiome.
**Ohad Manor and Elhanan Borenstein.** *Genome Biology*

==================
Question forum
==================
For EMPANADA announcements and questions, including notification of new releases, you can visit the `EMPANADA users forum <https://groups.google.com/forum/#!forum/empanada-users>`_.
