#!/usr/bin/env python

"""
This is the running script for EMPANADA
"""
# to comply with both Py2 and Py3
from __future__ import absolute_import, division, print_function

import argparse

# when the test module will be ready:
from empanada.compute_pathway_abundance import main

# for testing:
#import sys
#sys.path.append('/net/gs/vol1/home/ohadm/BorensteinLab/PROJECTS/FUNC_VAR_EMPANADA_OM/GitHub/empanada/empanada')
#from compute_pathway_abundance import main

if __name__ == "__main__":
    # get options from user
    parser = argparse.ArgumentParser(description='Compute the abundance of pathways in metagenomic data')

    # Required arguments:

    parser.add_argument('-ko', '--ko_abundance', dest='ko_abun_file', help='Input file of ko abundance per sample', default=None)

    parser.add_argument('-ko2path', '--ko_to_pathway', dest='ko_to_pathway_file', help='Input file of mapping from ko to pathway', default=None)

    # Optional arguments:

    parser.add_argument('-o', '--output', dest='output_file', help='Output file for resulting pathway abundance (default: out.tab)', default='out.tab')

    parser.add_argument('-oc', '--output_counts', dest='output_counts_file', help='Output file for number of KOs mapped to each pathway (default: counts.tab)', default='counts.tab')

    parser.add_argument('-om', '--output_mapping', dest='output_mapping_table',
                        help='Output the mapping table (either given or generated) to file, works only with pooled mappings (default: mapping.tab)', default='mapping.tab')

    parser.add_argument('-map', '--mapping_method', dest='mapping_method', help='Method to map KOs to Pathway (default: naive)',
                        default='naive', choices=['naive', 'by_support', 'by_sum_abundance', 'by_avg_abundance'])

    parser.add_argument('-compute', '--compute_method', dest='compute_method', help='Method to compute pathway abundance from mapped KOs (default: sum)',
                        default='sum', choices=['sum'])

    parser.add_argument('-threshold', '--abundance_threshold', dest='abundance_threshold', help='Abundance threshold to include KOs (default: 0.0)',
                        default='0.0')

    parser.add_argument('-fraction', '--fractional_ko_contribution', dest='fractional_ko_contribution',
                        help='Divide KO contributions such that they sum to 1 for each KO (default: False)', action='store_true')

    parser.add_argument('-remove_ko_with_no_pathway', dest='remove_ko_with_no_pathway',
                        help='Remove KOs with no pathway from analysis (default: False)', action='store_true')

    parser.add_argument('-remove_ko_with_no_abundance_measurement', dest='remove_ko_with_no_abundance_measurement',
                        help='Remove KOs with no measurements in the abundance table from analysis (default: False)', action='store_true')

    parser.add_argument('-transpose_ko', '--transpose_ko_abundance', dest='transpose_ko_abundance',
                        help='Transpose the ko abundance matrix given (default: False)', action='store_true')

    parser.add_argument('-transpose_output', '--transpose_output', dest='transpose_output',
                        help='Transpose the output pathway abundance matrix (default: False)', action='store_true')

    parser.add_argument('-permute_ko_mapping', dest='permute_ko_mapping',
                        help='Permute the given KO mapping, i.e., which KO map to which pathways for hypothesis testing (default: False)', action='store_true')

    parser.add_argument('-use_only_non_overlapping_genes', dest='use_only_non_overlapping_genes',
                        help='If the mapping is by_abundance, compute pathway support by only using non-overlapping genes (default: False)', action='store_true')

    parser.add_argument('-pool_samples_use_median', dest='pool_samples_use_median',
                        help='If the mapping is by_abundance, pool samples together using the median KO abundance, and learn the mapping only once (default: False)', action='store_true')

    parser.add_argument('-pool_samples_use_average', dest='pool_samples_use_average',
                        help='If the mapping is by_abundance, pool samples together using the average KO abundance, and learn the mapping only once (default: False)', action='store_true')

    parser.add_argument('-leave_one_ko_out_pathway_support', dest='leave_one_ko_out_pathway_support',
                        help='If the mapping is by_abundance, compute pathway support for each KO separately by removing it from the computation (default: False)', action='store_true')

    parser.add_argument('-compute_support_with_weighted_double_counting', dest='compute_support_with_weighted_double_counting',
                        help='If the mapping is by_abundance, double count KO abundance (weighted by mapping) when computing pathway support  (default: False)', action='store_true')

    parser.add_argument('-v', '--verbose', dest='verbose', help='Increase verbosity of module (default: false)', action='store_true')

    given_args = parser.parse_args()
    main(vars(given_args))
