#!/usr/bin/env python


import helpers.helpers as helpers
import argparse, sys
import yaml
import os
import pandas as pd
from itertools import compress
import numpy as np
import time

# Example python snv_feast_cli.py --species Bacteroides_uniformis_57318 --min_reads 5 --start_index 1 --end_index 200 --config_file_path /Users/leahbriscoe/Documents/FEASTX/Signature_SNVs/configs/config.yaml


def main():

    with open(config_file_path, "r") as f:
        config = yaml.safe_load(f)

    snp_dir = config["snp_dir"]
    input_dir = config["input_dir"]
    output_dir = config["output_dir"]
    file_string = (
        "Signature_SNVs_"
        + species
        + "_region_"
        + str(start_index)
        + "_"
        + str(end_index)
    )

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(output_dir + "/" + species):
        os.makedirs(output_dir + "/" + species)

    # SPECFICY PATHS

    depth_path = snp_dir + "/" + species + "/" + "snps_depth.txt.bz2"
    freq_path = snp_dir + "/" + species + "/" + "snps_ref_freq.txt.bz2"
    # Get sample data
    sink_source_config = pd.read_csv(input_dir + "/" + "sink_source.csv")
    all_samples = sink_source_config.iloc[:, 1 : (sink_source_config.shape[1])].values
    all_samples = all_samples.flatten()
    all_source_samples = sink_source_config.iloc[
        :, 2 : (sink_source_config.shape[1])
    ].values
    all_sink_samples = sink_source_config.iloc[:, 1:2].values
    all_source_samples = sorted(list(set(all_source_samples.flatten())))
    all_sink_samples = sorted(list(set(all_sink_samples.flatten())))
    sink_source_config.index = sink_source_config["family_id"]
    unique_family_ids = [f for f in sink_source_config["family_id"]]

    # PRODUCE THE SNP DEPTH AND FREQ MATRICES
    (
        snp_depth_matrix,
        depth_snp_names,
        depth_sample_names,
    ) = helpers.make_np_array_from_file(
        depth_path, start_index, end_index, all_samples=all_samples
    )
    (
        snp_freq_matrix,
        freq_snp_names,
        freq_sample_names,
    ) = helpers.make_np_array_from_file(
        freq_path, start_index, end_index, all_samples=all_samples
    )

    # to make it alternative allele frequency

    snp_freq_matrix = 1 - snp_freq_matrix

    (
        snp_freq_matrix,
        freq_snp_names,
        snp_depth_matrix,
        depth_snp_names,
    ) = helpers.prevalence_filter(
        snp_freq_matrix,
        freq_snp_names,
        freq_sample_names,
        snp_depth_matrix,
        depth_snp_names,
        depth_sample_names,
        species,
        all_sink_samples,
    )

    all_fam_LR = dict()
    all_fam_present = dict()
    fam_count = 0
    for fam in unique_family_ids:
        print("Experiment number:" + str(fam))
        start_time = time.time()
        fam_count += 1
        fam_file_string = file_string + "_fam_" + str(fam)
        all_fam = list(sink_source_config.loc[fam][1 : sink_source_config.shape[1]])
        all_sources = list(sink_source_config.loc[fam][2 : sink_source_config.shape[1]])

        source_bit_mask = [s in all_sources for s in freq_sample_names]
        family_bit_mask = [s in all_fam for s in freq_sample_names]
        familywise_bit_mask = [s in freq_sample_names for s in all_fam]
        all_fam_present[fam] = pd.Series(familywise_bit_mask)
        sink_bit_mask = [s in all_fam[0:1] for s in freq_sample_names]

        # skip this iteration if none of infants are here
        if not any(
            familywise_bit_mask[0:1]
        ):  # seed1 is this, seed2 is any(familywise_bit_mask[0])
            print("sink sample is not here, skip analysis for this family")
            continue
        else:
            print("present")
        family_samples_present = [(s) for s in freq_sample_names if s in all_fam]
        # print("family present")
        # print(family_samples_present)

        if len(family_samples_present) < 2:
            print("too little")
            continue

        # IMPOSE FILTER BASED ON READ DEPTH, snp_freq_matrix, snp_depth_matrix
        snp_depth_matrix_family = snp_depth_matrix[:, family_bit_mask]
        snp_freq_matrix_family = snp_freq_matrix[:, family_bit_mask]

        snp_freq_matrix_sink = snp_freq_matrix[:, sink_bit_mask]
        snp_depth_matrix_sink = snp_depth_matrix[:, sink_bit_mask]

        min_depth_bit_mask = snp_depth_matrix_family >= min_reads
        num_cols = snp_depth_matrix_family.shape[1]

        feature_filter = np.where(
            (min_depth_bit_mask.sum(axis=1) >= num_cols)
            & (snp_freq_matrix_sink.sum(axis=1) > 0),
            True,
            False,
        )
        # print("LENGTH row names")
        # print(len(freq_snp_names))
        freq_snp_names_family = np.array(freq_snp_names)[feature_filter]
        # print("shape before filter")
        # print(snp_freq_matrix_family.shape)
        # print(snp_freq_matrix_family[0:10])
        # print(snp_depth_matrix_family[0:10])
        snp_freq_matrix_family_filtered = snp_freq_matrix_family[feature_filter, :]
        snp_depth_matrix_family_filtered = snp_depth_matrix_family[feature_filter, :]

        snp_freq_matrix_sink_filtered = snp_freq_matrix_sink[feature_filter, 0:1]
        snp_depth_matrix_sink_filtered = snp_depth_matrix_sink[feature_filter, 0:1]
        snp_freq_matrix_filtered = snp_freq_matrix[feature_filter, :]
        snp_depth_matrix_filtered = snp_depth_matrix[feature_filter, :]
        snp_freq_matrix_sources_filtered = snp_freq_matrix_filtered[:, source_bit_mask]
        snp_depth_matrix_sources_filtered = snp_depth_matrix_filtered[
            :, source_bit_mask
        ]

        likelihoods = helpers.score_snps(
            snp_freq_matrix_sink_filtered,
            snp_depth_matrix_sink_filtered,
            snp_freq_matrix_sources_filtered,
            snp_depth_matrix_sources_filtered,
        )
        # print("length likelihoods")
        # print(len(likelihoods))

        likelihood_for_cutoff = [l for l in likelihoods if l != -1]

        cutoff_best = 2 * np.std(likelihood_for_cutoff) + np.mean(
            likelihood_for_cutoff
        )  # sede 6 siulation
        print(cutoff_best)
        private_snp_ind = [
            n
            for n, i in enumerate(likelihoods)
            if ((i > cutoff_best and i > 1) or i == -1)
        ]

        # print("SD+2:  length private SNPs")
        # print(len(private_snp_ind))

        select_private_snps_bit_mask = [
            snp in private_snp_ind
            for snp in range(snp_freq_matrix_sources_filtered.shape[0])
        ]
        all_fam_LR[fam] = pd.Series(likelihoods)

        freq_row_snp_names_family = list(
            compress(freq_snp_names_family, select_private_snps_bit_mask)
        )
        snp_freq_matrix_family_filtered_private = snp_freq_matrix_family_filtered[
            select_private_snps_bit_mask, :
        ]
        snp_depth_matrix_family_filtered_private = snp_depth_matrix_family_filtered[
            select_private_snps_bit_mask, :
        ]

        snp_counts_matrix_family_filtered_private = np.round_(
            np.multiply(
                snp_freq_matrix_family_filtered_private,
                snp_depth_matrix_family_filtered_private,
            ),
            decimals=0,
        )
        ref_counts_matrix_family_filtered_private = np.round_(
            np.multiply(
                (1 - snp_freq_matrix_family_filtered_private),
                snp_depth_matrix_family_filtered_private,
            ),
            decimals=0,
        )

        # set to null no longer needed variables
        snp_freq_matrix_family_filtered = None
        snp_depth_matrix_family_filtered = None

        if snp_counts_matrix_family_filtered_private.shape[0] == 0:
            print("no signature SNV")
            continue

        family_private_snps_freq = pd.DataFrame(columns=all_fam, index=private_snp_ind)
        family_private_snps_freq.fillna(0)

        family_private_snps_counts = pd.DataFrame(
            columns=all_fam, index=private_snp_ind
        )
        family_private_snps_counts.fillna(0)
        family_private_refsnps_counts = pd.DataFrame(
            columns=all_fam, index=private_snp_ind
        )
        family_private_refsnps_counts.fillna(0)

        for c in range(len(family_samples_present)):
            family_private_snps_freq[
                family_samples_present[c]
            ] = snp_freq_matrix_family_filtered_private[:, c]
            family_private_snps_counts[
                family_samples_present[c]
            ] = snp_counts_matrix_family_filtered_private[:, c]
            family_private_refsnps_counts[
                family_samples_present[c]
            ] = ref_counts_matrix_family_filtered_private[:, c]

        snp_freq_matrix_family_filtered_private = None
        snp_counts_matrix_family_filtered_private = None
        ref_counts_matrix_family_filtered_private = None
        final_family_private_snps_freq = family_private_snps_freq.values

        final_family_private_snps_counts = family_private_snps_counts.values
        final_family_private_refsnps_counts = family_private_refsnps_counts.values

        final_family_private_allsnps_counts = np.concatenate(
            [final_family_private_snps_counts, final_family_private_refsnps_counts]
        )

        freq_row_snp_names_family1 = ["Alt_" + i for i in freq_row_snp_names_family]
        freq_row_snp_names_family2 = ["Ref_" + i for i in freq_row_snp_names_family]
        freq_row_snp_names_family = (
            freq_row_snp_names_family1 + freq_row_snp_names_family2
        )
        final_family_private_allsnps_counts = np.concatenate(
            (
                np.transpose(np.array([freq_row_snp_names_family])),
                final_family_private_allsnps_counts,
            ),
            axis=1,
        )
        final_family_private_allsnps_counts = pd.DataFrame(
            final_family_private_allsnps_counts
        )

        final_family_private_allsnps_counts.to_csv(
            output_dir + "/" + species + "/" + fam_file_string + "_counts.csv.bz2",
            compression="bz2",
            header=False,
            index=False,
        )
        print("--- %s seconds ---" % (time.time() - start_time))

    if len(all_fam_present) > 0:
        all_fam_presents = pd.concat(all_fam_present, axis=1)
        all_fam_presents.to_csv(
            output_dir
            + "/"
            + species
            + "/"
            + "Summary_Stats_"
            + file_string
            + "_experiments_present_in_metadata.csv",
            sep=",",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", help="species for this job")
    parser.add_argument(
        "--min_reads", help="Minimum read depth for a valid allele frequency", type=int
    )
    parser.add_argument(
        "--start_index",
        help="Starting index for checking the species file because big file",
        type=int,
    )
    parser.add_argument("--end_index", help="End index", type=int)
    parser.add_argument(
        "--config_file_path", help="Path absolute to config file", type=str
    )

    args = parser.parse_args()
    species = args.species
    min_reads = args.min_reads  # minimum read depth
    start_index = args.start_index  # following python system 0 based
    end_index = args.end_index
    config_file_path = args.config_file_path
    main()
