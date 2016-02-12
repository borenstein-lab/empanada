"""
This function aggregates the abundance of KOs into pathways or modules abundances
"""

# to comply with both Py2 and Py3
from __future__ import absolute_import, division, print_function

import argparse
import numpy as np
import pandas as pd
import os
import sys
import warnings


###################################################################################################################
# MAIN FUNCTION
###################################################################################################################
def main(args):

    ##################################################
    # some default args for testing
    # args = {'ko_to_pathway_file': '../MUSiCC/Matrices/KOvsMODULE_KEGG_2013_07_15', 'ko_abun_file': 'HMP_DATA/WGS_KO_vs_SAMPLE_MUSiCC.tab',
    #          'output_file': 'out.tab', 'mapping_method': 'naive', 'compute_method': 'sum',
    #          'transpose_ko_abundance': False, 'transpose_output': False, 'output_counts_file': 'counts.tab'}
    # ##################################################

    if 'verbose' in args.keys() and args['verbose']:
        print("Given parameters: ", args)

    if args['leave_one_ko_out_pathway_support'] and args['mapping_method'] != 'by_sum_abundance' and args['mapping_method'] != 'by_avg_abundance':
        print(args['mapping_method'])
        sys.exit('Error: -leave_one_ko_out_pathway_support only works with mapping_method == by_sum_abundance')

    ###################################################################################################################
    # INPUT
    ###################################################################################################################

    abundance_threshold = float(args['abundance_threshold'])

    print("Reading files...")

    if 'ko_abun_file' in args.keys() and args['ko_abun_file'] is not None:
        if not os.path.isfile(args['ko_abun_file']):
            sys.exit('Error: Input file "' + args['ko_abun_file'] + '" does not exist')
        ko_abun_data = pd.read_table(args['ko_abun_file'], index_col=0, dtype={0: str})

    elif 'ko_abun_pd' in args.keys():
        ko_abun_data = args['ko_abun_pd']

    else:
        sys.exit('Error: No input ko abundance file given to script')

    if args['transpose_ko_abundance']:
        ko_abun_data = ko_abun_data.T

    if 'output_file' in args.keys() and args['ko_to_pathway_file'] is not None:
        if not os.path.isfile(args['ko_to_pathway_file']):
            sys.exit('Error: Input file "' + args['ko_to_pathway_file'] + '" does not exist')
        ko_to_pathway_data = pd.read_table(args['ko_to_pathway_file'], index_col=0, dtype={0: str})

    else:
        if args['ko_to_pathway_pd'] is not None:
            ko_to_pathway_data = args['ko_to_pathway_pd']
        else:
            sys.exit('Error: No input ko to pathway file given to script')

    print("Done.")

    ###################################################################################################################
    # IF REQUESTED BY USER, PERMUTE THE KO TO PATHWAY MAPPING
    ###################################################################################################################
    if args['permute_ko_mapping']:
        permuted_ko_to_path_values = np.random.permutation(ko_to_pathway_data.values)
        ko_to_pathway_data = pd.DataFrame(data=permuted_ko_to_path_values, index=ko_to_pathway_data.index, columns=ko_to_pathway_data.columns)

    ###################################################################################################################
    # CREATE AN "UNKNOWN" PATHWAY FOR KOS THAT ARE NOT MAPPED TO ANY PATHWAY (DEFAULT)
    # OR
    # FILTER OUT KOs THAT ARE NOT PART OF ANY PATHWAY
    ###################################################################################################################
    if not args['remove_ko_with_no_pathway']:
        # find the list of KOs that do not have any pathway
        ko_not_in_path = np.setdiff1d(ko_abun_data.index.values, ko_to_pathway_data.index.values)
        # create an "unknown" pathway for these KOs
        print("Adding 'unknown' pathway assignment for " + str(len(ko_not_in_path)) + " KOs with no pathway...")
        original_ko_to_path_values = ko_to_pathway_data.values
        original_ko_to_path_values_index = ko_to_pathway_data.index
        original_ko_to_path_values_columns = ko_to_pathway_data.columns
        #print(len(original_ko_to_path_values_index))
        new_ko_to_path_values_index = np.hstack((original_ko_to_path_values_index, ko_not_in_path))
        #print(len(new_ko_to_path_values_index))
        new_ko_to_path_values_columns = np.hstack((original_ko_to_path_values_columns, ['unknown']))
        # now, for all the KOs we added to the table, set the "unknown" pathway to be 1
        append_ko_to_path_values = np.zeros((len(ko_not_in_path), len(new_ko_to_path_values_columns)))
        append_ko_to_path_values[:, (len(new_ko_to_path_values_columns)-1)] = 1
        #print(original_ko_to_path_values.shape)
        new_ko_to_path_values = np.hstack((original_ko_to_path_values, np.zeros((original_ko_to_path_values.shape[0], 1))))
        new_ko_to_path_values = np.vstack((new_ko_to_path_values, append_ko_to_path_values))
        # now, for all the KOs that are in the table but have 0 for all pathways (i.e., not mapped), set the "unknown" as 1
        new_ko_to_path_values[np.sum(new_ko_to_path_values, axis=1) == 0, (len(new_ko_to_path_values_columns)-1)] = 1
        #print(new_ko_to_path_values.shape)

        # create the new ko_to_pathway table:
        ko_to_pathway_data = pd.DataFrame(data=new_ko_to_path_values, index=new_ko_to_path_values_index, columns=new_ko_to_path_values_columns)

    if not args['remove_ko_with_no_abundance_measurement']:  # add KOs that have pathway information and no abundance as constant zero abundance
        # find the list of KOs that do not have abundance values
        ko_not_in_abun = np.setdiff1d(ko_to_pathway_data.index.values, ko_abun_data.index.values)
        # create an "unknown" pathway for these KOs
        print("Adding 'zero abundance' vectors for " + str(len(ko_not_in_abun)) + " KOs with no abundance values...")
        original_ko_abun_values = ko_abun_data.values
        original_ko_abun_values_index = ko_abun_data.index
        original_ko_abun_values_columns = ko_abun_data.columns
        #print(len(original_ko_abun_values_index))
        new_ko_abun_values_index = np.hstack((original_ko_abun_values_index, ko_not_in_abun))
        #print(len(new_ko_abun_values_index))
        new_ko_abun_values_columns = original_ko_abun_values_columns  # no change in columns/samples
        # now, for all the KOs we added to the abundance table, set all values to be zero abundance
        append_ko_abun_values = np.zeros((len(ko_not_in_abun), len(new_ko_abun_values_columns)))
        #print(original_ko_abun_values.shape)
        new_ko_abun_values = np.vstack((original_ko_abun_values, append_ko_abun_values))
        #print(new_ko_abun_values.shape)

        # create the new ko_abundance table:
        ko_abun_data = pd.DataFrame(data=new_ko_abun_values, index=new_ko_abun_values_index, columns=new_ko_abun_values_columns)

    # sort tables similarly
    #print(ko_abun_data.shape)
    #print(ko_to_pathway_data.shape)
    ko = np.sort(np.intersect1d(ko_abun_data.index.values, ko_to_pathway_data.index.values))
    ko_abun_data = ko_abun_data.loc[ko]
    ko_to_pathway_data = ko_to_pathway_data.loc[ko]
    #print(ko_abun_data.shape)
    #print(ko_to_pathway_data.shape)

    num_of_kos = ko_abun_data.shape[0]
    num_of_samples = ko_abun_data.shape[1]
    num_of_pathways = ko_to_pathway_data.shape[1]
    print("number of samples is " + str(num_of_samples))
    print("number of KOs is " + str(num_of_kos))
    print("number of pathways is " + str(num_of_pathways))

    ###################################################################################################################
    # ASSERTIONS
    ###################################################################################################################
    # make sure the abundance matrix and the ko mapping matrix have the exact same index
    assert(np.array_equal(ko_abun_data.index.values, ko_to_pathway_data.index.values))

    ###################################################################################################################
    # set KOs < abundance_threshold to be 0.0
    ###################################################################################################################
    #print(ko_abun_data)
    ko_abun_data__threshold = ko_abun_data.values
    ko_abun_data__threshold[ko_abun_data__threshold < abundance_threshold] = 0.0
    ko_abun_data = pd.DataFrame(data=ko_abun_data__threshold, index=ko_abun_data.index, columns=ko_abun_data.columns)
    #print(ko_abun_data)

    ###################################################################################################################
    # FOR EACH PATHWAY, COUNT THE NUMBER OF NON-ZERO KOS THAT MAP TO IT IN EACH SAMPLE
    # AND CALCULATE THE SUPPORT FOR EACH PATHWAY AS THE FRACTION OF NON-ZERO KOS FOR THAT PATHWAY
    ###################################################################################################################
    ko_binary_data = (ko_abun_data.values > 0).astype(float)
    ko_to_pathway_binary_data = (ko_to_pathway_data.values > 0).astype(float)
    pathway_counts = np.dot(ko_binary_data.T, ko_to_pathway_binary_data).T

    #print(np.sum(ko_to_pathway_binary_data, axis=0))
    pathway_support = pathway_counts / np.sum(ko_to_pathway_binary_data, axis=0)[:, None]
    pathway_support[np.isnan(pathway_support)] = 0.0

    #np.set_printoptions(threshold=np.nan)
    #print(pathway_counts)
    #print(pathway_support)

    ###################################################################################################################
    # COMPUTE THE PATHWAY ABUNDANCES
    ###################################################################################################################
    if args['mapping_method'] == 'naive':
        # first, if user requested a fractional version, divide each KO row by its sum
        if args['fractional_ko_contribution']:
            ko_to_pathway_data_fractional = ko_to_pathway_data.values / np.sum(ko_to_pathway_data.values, axis=1)[:, None]
            ko_to_pathway_data_fractional[np.isnan(ko_to_pathway_data_fractional)] = 0.0
            ko_to_pathway_data = pd.DataFrame(data=ko_to_pathway_data_fractional, index=ko_to_pathway_data.index, columns=ko_to_pathway_data.columns)

        # now, given the method requested, create the pathway abundance matrix
        if args['compute_method'] == 'sum':
            pathway_abundance = np.dot(ko_abun_data.T, ko_to_pathway_data).T

    # for each sample, we use the pathway support to compute the weight for each KO
    elif args['mapping_method'] == 'by_support':
        pathway_abundance = np.zeros((num_of_pathways, num_of_samples))
        if args['compute_method'] == 'sum':
            for s in range(num_of_samples):  # num_of_samples
                curr_sample_ko_abun = ko_abun_data.values[:, s].reshape((1, num_of_kos))
                curr_sample_support = pathway_support[:, s]
                curr_sample_support_ko_to_pathway_weighted = ko_to_pathway_binary_data * curr_sample_support[:, None].T

                # if user requested a fractional version, divide each KO row by its sum
                if args['fractional_ko_contribution']:
                    curr_sample_support_ko_to_pathway_weighted = curr_sample_support_ko_to_pathway_weighted / np.sum(curr_sample_support_ko_to_pathway_weighted, axis=1)[:, None]
                    curr_sample_support_ko_to_pathway_weighted[np.isnan(curr_sample_support_ko_to_pathway_weighted)] = 0.0

                pathway_abundance[:, s] = np.dot(curr_sample_ko_abun, curr_sample_support_ko_to_pathway_weighted).T.flatten()

    # for each sample, we use the pathway KO abundance to compute the weight for each KO
    # this is an iterative process since the pathways' abundance changes every iteration
    elif args['mapping_method'] == 'by_sum_abundance' or args['mapping_method'] == 'by_avg_abundance':
        # first make sure the user chose to make the mapping fractional
        if not args['fractional_ko_contribution']:
            print("Mapping by abundance is only available with -fraction")
            exit()

        np.set_printoptions(threshold=np.nan)

        num_of_pathway_per_ko = np.sum(ko_to_pathway_binary_data, axis=1).reshape((num_of_kos, 1))

        if args['use_only_non_overlapping_genes']:
            non_overlapping_kos = (num_of_pathway_per_ko == 1)
            non_overlapping_kos_names = ko_abun_data.index[non_overlapping_kos.flatten()]
            assert(np.sum(non_overlapping_kos) == len(non_overlapping_kos_names))
            print("num of overlapping KOs = " + str(np.sum(~non_overlapping_kos)))
            print("num of non-overlapping KOs = " + str(np.sum(non_overlapping_kos)))
            num_of_kos_per_pathway = np.sum(ko_to_pathway_binary_data[non_overlapping_kos.flatten(), :], axis=0)
        else:
            num_of_kos_per_pathway = np.sum(ko_to_pathway_binary_data, axis=0)

        if args['leave_one_ko_out_pathway_support']:
            num_of_kos_per_pathway -= 1
            num_of_kos_per_pathway[num_of_kos_per_pathway < 0] = 0

        #print(num_of_kos_per_pathway)

        if args['pool_samples_use_median'] or args['pool_samples_use_average']:

            if args['pool_samples_use_median']:
                print("Optimizing mapping on median of pooled samples...")
                pooled_ko_abun = np.median(ko_abun_data.values, axis=1).reshape((num_of_kos, 1))
            else:  # pool using average
                print("Optimizing mapping on mean of pooled samples...")
                pooled_ko_abun = np.mean(ko_abun_data.values, axis=1).reshape((num_of_kos, 1))

            pooled_samples_ko_to_pathway_mapping = ko_to_pathway_binary_data / np.sum(ko_to_pathway_binary_data, axis=1)[:, None]
            # verify that each row sums to 1.0 and non-negative
            if not 0.999 < np.sum(pooled_samples_ko_to_pathway_mapping, axis=1).all() < 1.001:
                print("\npooled_samples_ko_to_pathway_mapping does not sum to one in all rows:")
                print(np.sum(pooled_samples_ko_to_pathway_mapping, axis=1))
                exit()
            if not 0.0 <= pooled_samples_ko_to_pathway_mapping.all() <= 1.0:
                print("\npooled_samples_ko_to_pathway_mapping is not between 0 and 1:")
                exit()

            # now use the current mapping to create the pathway abundance vector
            if args['use_only_non_overlapping_genes']:  # use only non-overlapping KOs to compute pathway abundance/support
                pooled_samples_pathway_abun = np.dot(pooled_samples_ko_to_pathway_mapping[non_overlapping_kos.flatten(), :].T, pooled_ko_abun[non_overlapping_kos.flatten(), :])
            else:
                pooled_samples_pathway_abun = np.dot(pooled_samples_ko_to_pathway_mapping.T, pooled_ko_abun)

            #print(non_overlapping_kos.flatten()[0:10])
            #print(pooled_samples_pathway_abun[0:10])
            print("sum(pooled_samples_pathway_abun) = " + str(np.sum(pooled_samples_pathway_abun)))

            # now use the pathway abun to create the target fractions for each KO, given its binary mapping
            pathway_support = pooled_samples_pathway_abun.flatten()

            #print(pathway_support)
            target_ko_fraction_for_each_pathway = ko_to_pathway_binary_data * pathway_support
            # if we do not want to include each KO in its own computation, then subtract the KO abun (using previous mapping) from the support:
            if args['leave_one_ko_out_pathway_support']:
                if args['use_only_non_overlapping_genes']:  # use only non-overlapping KOs to compute pathway abundance/support
                    target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway - (pooled_samples_ko_to_pathway_mapping * (pooled_ko_abun * non_overlapping_kos))
                    # zero out columns of pathways with zero KOs contributing to them
                    target_ko_fraction_for_each_pathway[:, (num_of_kos_per_pathway == 0)] = 0.0
                else:
                    target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway - (pooled_samples_ko_to_pathway_mapping * pooled_ko_abun)

            #print(target_ko_fraction_for_each_pathway[0, :])
            # if we want to use the average pathway abundance as evidence, then divide by the number of KOs in pathway
            if args['mapping_method'] == 'by_avg_abundance':
                target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway / num_of_kos_per_pathway
                target_ko_fraction_for_each_pathway[np.isnan(target_ko_fraction_for_each_pathway)] = 0.0
            #print(target_ko_fraction_for_each_pathway[0, :])
            # for rows in the target ko to pathway that sum to zero (i.e., no information) copy the original row from the binary matrix
            target_ko_fraction_for_each_pathway[np.sum(target_ko_fraction_for_each_pathway, axis=1) == 0, :] = ko_to_pathway_binary_data[np.sum(target_ko_fraction_for_each_pathway, axis=1) == 0, :]
            target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway / np.sum(target_ko_fraction_for_each_pathway, axis=1)[:, None]
            target_ko_fraction_for_each_pathway[np.isnan(target_ko_fraction_for_each_pathway)] = 0.0

            # compute the current error given the target and current fractions:
            mapping_error_per_ko = np.sum(abs(pooled_samples_ko_to_pathway_mapping - target_ko_fraction_for_each_pathway), axis=1)
            print("starting error:" + str(np.sum(mapping_error_per_ko[pooled_ko_abun.flatten() > 0])))

            counter = 0
            while counter < 100 and np.sum(mapping_error_per_ko) > 0.001:
                counter += 1

                # now, set the current mapping to be the target fraction
                pooled_samples_ko_to_pathway_mapping = target_ko_fraction_for_each_pathway
                # verify that each row sums to 1.0
                if not 0.999 < np.sum(pooled_samples_ko_to_pathway_mapping[pooled_ko_abun.flatten() > 0, :], axis=1).all() < 1.001:
                    print("\npooled_samples_ko_to_pathway_mapping does not sum to one in all rows:")
                    print(np.where(np.sum(pooled_samples_ko_to_pathway_mapping[pooled_ko_abun.flatten() > 0, :], axis=1) <= 0.999))
                    print(np.where(np.sum(pooled_samples_ko_to_pathway_mapping[pooled_ko_abun.flatten() > 0, :], axis=1) >= 1.001))
                    print(target_ko_fraction_for_each_pathway[0,:])
                    exit()
                if not 0.0 <= pooled_samples_ko_to_pathway_mapping.all() <= 1.0:
                    print("\npooled_samples_ko_to_pathway_mapping is not between 0 and 1:")
                    exit()

                # recalculate the pathway abundance given the new mapping
                if args['use_only_non_overlapping_genes']:  # use only non-overlapping KOs to compute pathway abundance/support
                    pooled_samples_pathway_abun = np.dot(pooled_samples_ko_to_pathway_mapping[non_overlapping_kos.flatten(), :].T, pooled_ko_abun[non_overlapping_kos.flatten(), :])
                else:
                    pooled_samples_pathway_abun = np.dot(pooled_samples_ko_to_pathway_mapping.T, pooled_ko_abun)
                #print("sum(pooled_samples_pathway_abun) = " + str(np.sum(pooled_samples_pathway_abun)))

                # now use the pathway abun to create the target fractions for each KO, given its binary mapping
                pathway_support = pooled_samples_pathway_abun.flatten()
                #print(pathway_support)
                target_ko_fraction_for_each_pathway = ko_to_pathway_binary_data * pathway_support
                # if we do not want to include each KO in its own computation, then subtract the KO abun (using previous mapping) from the support:
                if args['leave_one_ko_out_pathway_support']:
                    if args['use_only_non_overlapping_genes']:  # use only non-overlapping KOs to compute pathway abundance/support
                        target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway - (pooled_samples_ko_to_pathway_mapping * (pooled_ko_abun * non_overlapping_kos))
                        # zero out columns of pathways with zero KOs contributing to them
                        target_ko_fraction_for_each_pathway[:, (num_of_kos_per_pathway == 0)] = 0.0
                    else:
                        target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway - (pooled_samples_ko_to_pathway_mapping * pooled_ko_abun)
                # if we want to use the average pathway abundance as evidence, then divide by the number of KOs in pathway
                if args['mapping_method'] == 'by_avg_abundance':
                    target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway / num_of_kos_per_pathway
                    target_ko_fraction_for_each_pathway[np.isnan(target_ko_fraction_for_each_pathway)] = 0.0

                # for rows in the target ko to pathway that sum to zero (i.e., no information) copy the original row from the binary matrix
                target_ko_fraction_for_each_pathway[np.sum(target_ko_fraction_for_each_pathway, axis=1) == 0, :] = ko_to_pathway_binary_data[np.sum(target_ko_fraction_for_each_pathway, axis=1) == 0, :]
                target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway / np.sum(target_ko_fraction_for_each_pathway, axis=1)[:, None]
                target_ko_fraction_for_each_pathway[np.isnan(target_ko_fraction_for_each_pathway)] = 0.0

                # compute the current error given the target and current fractions:
                mapping_error_per_ko = np.sum(abs(pooled_samples_ko_to_pathway_mapping - target_ko_fraction_for_each_pathway), axis=1)
                print('.',end="",flush=True)
                #print("error:" + str(np.sum(mapping_error_per_ko[pooled_ko_abun.flatten() > 0])))

            print("\noptimized error:" + str(np.sum(mapping_error_per_ko[pooled_ko_abun.flatten() > 0])))

            if args['compute_method'] == 'sum':
                pathway_abundance = np.dot(ko_abun_data.T, pooled_samples_ko_to_pathway_mapping).T

            # save Ko to pathway mapping file
            ko_to_pathway_data = pd.DataFrame(data=pooled_samples_ko_to_pathway_mapping, index=ko_to_pathway_data.index, columns=ko_to_pathway_data.columns)

        else:  # learn mapping for each sample separately

            pathway_abundance = np.zeros((num_of_pathways, num_of_samples))

            for s in range(num_of_samples):  # num_of_samples
                print("Optimizing sample " + str(s+1) + "/" + str(num_of_samples))
                # get the current KO abundance vector
                curr_sample_ko_abun = ko_abun_data.values[:, s].reshape((num_of_kos, 1))
                # get the initial mapping by using a uniform distribution across pathways
                curr_sample_ko_to_pathway_mapping = ko_to_pathway_binary_data / np.sum(ko_to_pathway_binary_data, axis=1)[:, None]
                # verify that each row sums to 1.0
                if not 0.999 < np.sum(curr_sample_ko_to_pathway_mapping, axis=1).all() < 1.001:
                    print("\ncurr_sample_ko_to_pathway_mapping does not sum to one in all rows:")
                    print(np.sum(curr_sample_ko_to_pathway_mapping, axis=1))
                    exit()
                if not 0.0 <= curr_sample_ko_to_pathway_mapping.all() <= 1.0:
                    print("\npooled_samples_ko_to_pathway_mapping is not between 0 and 1:")
                    exit()

                # now use the current mapping to create the pathway abundance vector
                if args['use_only_non_overlapping_genes']:  # use only non-overlapping KOs to compute pathway abundance/support
                    curr_sample_pathway_abun = np.dot(curr_sample_ko_to_pathway_mapping[non_overlapping_kos.flatten(), :].T, curr_sample_ko_abun[non_overlapping_kos.flatten(), :])
                else:
                    curr_sample_pathway_abun = np.dot(curr_sample_ko_to_pathway_mapping.T, curr_sample_ko_abun)

                #print(curr_sample_pathway_abun)
                print("sum(pathway_abun) = " + str(np.sum(curr_sample_pathway_abun)))

                # now use the pathway abun to create the target fractions for each KO, given its binary mapping
                target_ko_fraction_for_each_pathway = ko_to_pathway_binary_data * curr_sample_pathway_abun.flatten()
                # if we do not want to include each KO in its own computation, then subtract the KO abundance from the support:
                if args['leave_one_ko_out_pathway_support']:
                    if args['use_only_non_overlapping_genes']:  # use only non-overlapping KOs to compute pathway abundance/support
                        target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway - (curr_sample_ko_to_pathway_mapping * (curr_sample_ko_abun * non_overlapping_kos))
                        # zero out columns of pathways with zero KOs contributing to them
                        target_ko_fraction_for_each_pathway[:, (num_of_kos_per_pathway == 0)] = 0.0
                    else:
                        target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway - (curr_sample_ko_to_pathway_mapping * curr_sample_ko_abun)

                # if we want to use the average pathway abundance as evidence, then divide by the number of KOs in pathway
                if args['mapping_method'] == 'by_avg_abundance':
                    target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway / num_of_kos_per_pathway
                    target_ko_fraction_for_each_pathway[np.isnan(target_ko_fraction_for_each_pathway)] = 0.0

                # for rows in the target ko to pathway that sum to zero (i.e., no information) copy the original row from the binary matrix
                target_ko_fraction_for_each_pathway[np.sum(target_ko_fraction_for_each_pathway, axis=1) == 0, :] = ko_to_pathway_binary_data[np.sum(target_ko_fraction_for_each_pathway, axis=1) == 0, :]
                target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway / np.sum(target_ko_fraction_for_each_pathway, axis=1)[:, None]
                target_ko_fraction_for_each_pathway[np.isnan(target_ko_fraction_for_each_pathway)] = 0.0

                # compute the current error given the target and current fractions:
                mapping_error_per_ko = np.sum(abs(curr_sample_ko_to_pathway_mapping - target_ko_fraction_for_each_pathway), axis=1)
                print("starting error:" + str(np.sum(mapping_error_per_ko[curr_sample_ko_abun.flatten() > 0])))

                counter = 0
                while counter < 100 and np.sum(mapping_error_per_ko) > 0.001:
                    counter += 1

                    # now, set the current mapping to be the target fraction
                    curr_sample_ko_to_pathway_mapping = target_ko_fraction_for_each_pathway
                    # verify that each row sums to 1.0
                    if not 0.999 < np.sum(curr_sample_ko_to_pathway_mapping[curr_sample_ko_abun.flatten() > 0, :], axis=1).all() < 1.001:
                        print("\ncurr_sample_ko_to_pathway_mapping does not sum to one in all rows:")
                        print(np.where(np.sum(curr_sample_ko_to_pathway_mapping[curr_sample_ko_abun.flatten() > 0, :], axis=1) <= 0.999))
                        print(np.where(np.sum(curr_sample_ko_to_pathway_mapping[curr_sample_ko_abun.flatten() > 0, :], axis=1) >= 1.001))
                        exit()
                    if not 0.0 <= curr_sample_ko_to_pathway_mapping.all() <= 1.0:
                        print("\npooled_samples_ko_to_pathway_mapping is not between 0 and 1:")
                        exit()

                    # recalculate the pathway abundance given the new mapping
                    if args['use_only_non_overlapping_genes']:  # use only non-overlapping KOs to compute pathway abundance/support
                        curr_sample_pathway_abun = np.dot(curr_sample_ko_to_pathway_mapping[non_overlapping_kos.flatten(), :].T, curr_sample_ko_abun[non_overlapping_kos.flatten(), :])
                    else:
                        curr_sample_pathway_abun = np.dot(curr_sample_ko_to_pathway_mapping.T, curr_sample_ko_abun)
                    #print("sum(pathway_abun) = " + str(np.sum(curr_sample_pathway_abun)))

                    # now use the pathway abun to create the target fractions for each KO, given its binary mapping
                    if args['compute_support_with_weighted_double_counting']:
                        curr_sample_ko_to_pathway_mapping_times_num_pathways = curr_sample_ko_to_pathway_mapping * num_of_pathway_per_ko
                        pathway_support = np.dot(curr_sample_ko_to_pathway_mapping_times_num_pathways.T, curr_sample_ko_abun).flatten()
                    else:
                        pathway_support = curr_sample_pathway_abun.flatten()

                    #print(pathway_support)
                    target_ko_fraction_for_each_pathway = ko_to_pathway_binary_data * pathway_support
                    # if we do not want to include each KO in its own computation, then subtract the KO abundance from the support:
                    if args['leave_one_ko_out_pathway_support']:
                        if args['use_only_non_overlapping_genes']:  # use only non-overlapping KOs to compute pathway abundance/support
                            target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway - (curr_sample_ko_to_pathway_mapping * (curr_sample_ko_abun * non_overlapping_kos))
                            # zero out columns of pathways with zero KOs contributing to them
                            target_ko_fraction_for_each_pathway[:, (num_of_kos_per_pathway == 0)] = 0.0
                        elif args['compute_support_with_weighted_double_counting']:
                            target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway - (curr_sample_ko_to_pathway_mapping_times_num_pathways * curr_sample_ko_abun)
                        else:
                            target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway - (curr_sample_ko_to_pathway_mapping * curr_sample_ko_abun)

                    # if we want to use the average pathway abundance as evidence, then divide by the number of KOs in pathway
                    if args['mapping_method'] == 'by_avg_abundance':
                        target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway / num_of_kos_per_pathway
                        target_ko_fraction_for_each_pathway[np.isnan(target_ko_fraction_for_each_pathway)] = 0.0

                    # for rows in the target ko to pathway that sum to zero (i.e., no information) copy the original row from the binary matrix
                    target_ko_fraction_for_each_pathway[np.sum(target_ko_fraction_for_each_pathway, axis=1) == 0, :] = ko_to_pathway_binary_data[np.sum(target_ko_fraction_for_each_pathway, axis=1) == 0, :]
                    target_ko_fraction_for_each_pathway = target_ko_fraction_for_each_pathway / np.sum(target_ko_fraction_for_each_pathway, axis=1)[:, None]
                    target_ko_fraction_for_each_pathway[np.isnan(target_ko_fraction_for_each_pathway)] = 0.0

                    # compute the current error given the target and current fractions:
                    mapping_error_per_ko = np.sum(abs(curr_sample_ko_to_pathway_mapping - target_ko_fraction_for_each_pathway), axis=1)
                    print('.',end="",flush=True)
                    #print("error:" + str(np.sum(mapping_error_per_ko[curr_sample_ko_abun.flatten() > 0])))

                print("\noptimized error:" + str(np.sum(mapping_error_per_ko[curr_sample_ko_abun.flatten() > 0])))
                if args['compute_method'] == 'sum':
                    pathway_abundance[:, s] = np.dot(curr_sample_ko_to_pathway_mapping.T, curr_sample_ko_abun).flatten()

            # save Ko to pathway mapping file
            ko_to_pathway_data = pd.DataFrame(data=curr_sample_ko_to_pathway_mapping, index=ko_to_pathway_data.index, columns=ko_to_pathway_data.columns)

    ###################################################################################################################
    # WRITE OUTPUT
    ###################################################################################################################

    print("Writing output...")

    if args['transpose_output']:
        path_pd = pd.DataFrame(data=pathway_abundance.T, index=ko_abun_data.columns, columns=ko_to_pathway_data.columns)

    else:
        path_pd = pd.DataFrame(data=pathway_abundance, index=ko_to_pathway_data.columns, columns=ko_abun_data.columns)

    path_pd.index.name = 'Pathway'

    if 'output_file' in args.keys():
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            path_pd.to_csv(args['output_file'], sep='\t')

    elif 'output_pd' in args.keys():
        args['output_pd'] = path_pd

    else:
        sys.exit('Error: No output destination given')

    if 'output_counts_file' in args.keys():
        counts_pd = pd.DataFrame(data=pathway_counts, index=ko_to_pathway_data.columns, columns=ko_abun_data.columns)
        counts_pd.index.name = 'Pathway'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            counts_pd.to_csv(args['output_counts_file'], sep='\t')

    if 'output_mapping_table' in args.keys():
        ko_to_pathway_data.index.name = 'KO'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            ko_to_pathway_data.to_csv(args['output_mapping_table'], sep='\t')

    print("Done.")

################################################################################################################

if __name__ == "__main__":
    # get options from user
    parser = argparse.ArgumentParser(description='Compute the abundance of pathways in metagenomic data')

    parser.add_argument('-ko', '--ko_abundance', dest='ko_abun_file', help='Input file of ko abundance per sample', default=None)

    parser.add_argument('-ko2path', '--ko_to_pathway', dest='ko_to_pathway_file', help='Input file of mapping from ko to pathway', default=None)

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






