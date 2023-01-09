# python nayfach_main.py --strain Bacteroides_uniformis_57318 --seed 27 --local_bool
# python nayfach_main.py --strain Bacteroides_uniformis_57318 --seed 27

import argparse, sys
import numpy as np
import pandas as pd
import os
from nayfach_helpers import (
    get_column_names,
    read_two_files,
    check_marker_allele,
    get_needed_accessions,
    write_dicts,
)
import time
from itertools import compress
import random


def main():

    if local_bool:
        backhed_input_dir = "/Users/leahbriscoe/Documents/FEASTX/BackhedFiles/"
        backhed_snp_dir = backhed_input_dir + "snps/" + strain + "/"

        simulation_input_dir = "/Users/leahbriscoe/Documents/FEASTX/SimulationDFiles/"
        simulation_snp_dir = simulation_input_dir + "midas_output_snps/" + strain + "/"

    else:
        backhed_input_dir = "/u/project/ngarud/daisyche/mother_infant/data/"
        backhed_snp_dir = backhed_input_dir + "snps/" + strain + "/"

        simulation_input_dir = "/u/project/ngarud/briscoel/FEASTX/SimulationDFiles/"
        simulation_snp_dir = simulation_input_dir + "midas_output_snps/" + strain + "/"

    simulation_marker_snv_dir = simulation_input_dir + "marker_snvs/" + strain + "/"
    if not os.path.exists(simulation_marker_snv_dir):
        os.makedirs(simulation_marker_snv_dir)
    backhed_depth_path = backhed_snp_dir + "snps_depth.txt.bz2"
    backhed_freq_path = backhed_snp_dir + "snps_ref_freq.txt.bz2"

    simulation_depth_path = simulation_snp_dir + "snps_depth.txt.bz2"
    simulation_freq_path = simulation_snp_dir + "snps_ref_freq.txt.bz2"

    inclusion_list = pd.read_csv(
        simulation_input_dir + "Nayfach_analysis/all_inclusion.csv"
    )
    inclusion_list = list(inclusion_list.iloc[:, 0])

    family_record = pd.read_csv(
        simulation_input_dir
        + "sink_source_config/Pipe_ChosenFamily_seed"
        + str(seed)
        + ".csv"
    )

    (
        infant_list,
        mother_list,
        infant_mother_pair_list,
        background_list,
    ) = get_needed_accessions(
        metadata=family_record,
        inclusion_list=inclusion_list,
        data1_path=simulation_depth_path,
        data2_path=backhed_depth_path,
    )

    (
        infant_marker_allele_dict,
        mother_marker_allele_dict,
        pair_marker_allele_dict,
    ) = read_two_files(
        depth_path_1=simulation_depth_path,
        freq_path_1=simulation_freq_path,
        depth_path_2=backhed_depth_path,
        freq_path_2=backhed_freq_path,
        infant_samples=infant_list,
        mother_samples=mother_list,
        infant_mother_pairs=infant_mother_pair_list,
        background_samples=background_list,
    )

    write_dicts(
        infant_marker_allele_dict,
        simulation_marker_snv_dir + "infant_marker_allele_strain__" + strain,
    )
    write_dicts(
        mother_marker_allele_dict,
        simulation_marker_snv_dir + "mother_marker_allele_strain__" + strain,
    )
    write_dicts(
        pair_marker_allele_dict,
        simulation_marker_snv_dir + "mother_infant_marker_allele_strain__" + strain,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--strain", help="Strain for this job")
    parser.add_argument("--seed", help="Seed for random fam draw", type=int)
    parser.add_argument(
        "--local_bool",
        help="Running locally (0 or 1)",
        default=False,
        action="store_true",
    )

    args = parser.parse_args()
    strain = args.strain
    seed = args.seed
    local_bool = args.local_bool
    main()

