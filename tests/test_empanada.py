#!/usr/bin/env python

"""
This is the testing unit for EMPANADA
"""
# to comply with both Py2 and Py3
from __future__ import absolute_import, division, print_function

import unittest
import subprocess
import os
import pandas as pd

# when the test module will be ready:
import empanada
from empanada.compute_pathway_abundance import main

# for testing:
#import sys
#sys.path.append('/net/gs/vol1/home/ohadm/BorensteinLab/PROJECTS/FUNC_VAR_EMPANADA_OM/GitHub/empanada/empanada')
#from compute_pathway_abundance import main


class EmpanadaTestCase(unittest.TestCase):
    """Tests for `compute_pathway_abundance.py`."""

    # Get the path to the location of the package:
    path_to_data = os.path.dirname(empanada.__file__)

    # for testing:
    # path_to_data = '/net/gs/vol1/home/ohadm/BorensteinLab/PROJECTS/FUNC_VAR_EMPANADA_OM/GitHub/empanada/empanada'

    def test_is_output_correct_for_empanada(self):
        """Does EMPANADA produce the correct output for the example case?"""
        print("Path to examples:" + EmpanadaTestCase.path_to_data)

        # run EMPANADA from shell to make sure it works properly
        # python ~/METAFIT/PyCode/FiShTaCo/fishtaco/compute_pathway_abundance.py -ko Data/t2d_ko_vs_sample.tab -ko2path Data/KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15
        # -o Data/t2d_pathway_vs_sample_BY_AVG_ABUNDANCE_NON_OVERLAPPING_LKO_FRACTIONAL.tab -threshold 0 -map by_avg_abundance -fraction
        # -leave_one_ko_out_pathway_support -use_only_non_overlapping_genes; \
        subprocess.call('run_empanada.py -ko ' + EmpanadaTestCase.path_to_data + '/examples/t2d_ko_vs_sample.tab ' +
                        ' -ko2path ' + EmpanadaTestCase.path_to_data + '/data/KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15.tab' +
                        ' -o test_empanada.tab -threshold 0 -map by_avg_abundance -fraction' +
                        ' -leave_one_ko_out_pathway_support -use_only_non_overlapping_genes', shell=True)

        print("Testing output...")

        # assert that the result is equal to the example (up to small difference due to platform and operating system)
        example = pd.read_table(EmpanadaTestCase.path_to_data + '/examples/empanada.tab', index_col=0)
        output = pd.read_table('test_empanada.tab', index_col=0)
        example_vals = example.values
        output_vals = output.values
        self.assertTrue(example_vals.shape[0] == output_vals.shape[0])
        self.assertTrue(example_vals.shape[1] == output_vals.shape[1])
        for i in range(example_vals.shape[0]):
            for j in range(example_vals.shape[1]):
                self.assertTrue(abs(example_vals[i, j] - output_vals[i, j]) < 0.1)

        print("Deleting temporary files...")
        os.remove('test_empanada.tab')
        os.remove('counts.tab')
        os.remove('mapping.tab')

################################################

if __name__ == '__main__':
    unittest.main()


