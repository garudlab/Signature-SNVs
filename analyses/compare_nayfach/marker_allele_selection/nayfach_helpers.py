from bz2 import BZ2File as bzopen
import numpy as np
import pickle
import pandas as pd

# (samples1_depth_vec >= min_depth).sum() == 0)
# [*filter(lambda x: x >= 30, my_list1)]


def pad_zeroes(the_dict):
    max_len_list = len(max(the_dict.values(), key=len))
    for k in the_dict.keys():
        padding_length = max_len_list - len(the_dict[k])
        if padding_length > 0:
            padding = ["" for p in range(padding_length)]
            the_dict[k].extend(padding)
    return the_dict


def write_dicts(the_dict, path):

    with open(path + ".pickle", "wb") as handle:
        pickle.dump(the_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if len(the_dict) > 0:
        the_dict = pad_zeroes(the_dict)

    df = pd.DataFrame.from_dict(the_dict)
    df.to_csv(path + ".csv", index=False)


def get_column_names(path):

    with bzopen(path, "r") as bzfin:
        """ Handle lines here """

        for i, line in enumerate(bzfin):
            if i > 0:
                break
            return line.decode().rstrip().split("\t")


def get_needed_accessions(metadata, inclusion_list, data1_path, data2_path):

    unique_infants = list(set(metadata.iloc[:, 1:2].values.flatten()))
    unique_mothers = list(
        set(metadata.iloc[:, 3 : (metadata.shape[1])].values.flatten())
    )

    unique_mother_infant_pairs = []
    for r in range(metadata.shape[0]):
        infant_r = metadata.iloc[r, 1:2].values[0]
        for c in range(metadata.shape[1]):
            unique_mother_infant_pairs.append(
                [infant_r, metadata.iloc[r, c : (c + 1)].values[0]]
            )

    # Get column head only to filter down the list
    data1_columns = get_column_names(data1_path)
    data2_columns = get_column_names(data2_path)
    data1_columns = data1_columns[1 : len(data1_columns)]
    data2_columns = data2_columns[1 : len(data2_columns)]

    unique_infants = list(set(unique_infants).intersection(set(data1_columns)))

    unique_mothers = list(set(unique_mothers).intersection(set(data2_columns)))

    background_list = list(set(inclusion_list).intersection(set(data2_columns)))

    new_unique_mother_infant_pairs = []
    for pair in unique_mother_infant_pairs:
        if pair[0] in unique_infants and pair[1] in unique_mothers:
            new_unique_mother_infant_pairs.append(pair)

    return (
        unique_infants,
        unique_mothers,
        new_unique_mother_infant_pairs,
        background_list,
    )


"""
Array operations to return the marker allels
"""


def process_allele_data(
    focal_sample_names,
    focal_depth,
    focal_freq,
    min_depth,
    min_freq,
    background_freq,
    background_depth,
):
    samples_with_marker_allele_ref = []
    samples_with_marker_allele_alt = []
    ## Focal processing
    focal_depth = np.array(focal_depth)
    focal_freq = np.array(focal_freq)

    focal_ref_allele_present_cond1 = focal_depth >= min_depth
    focal_ref_allele_present_cond2 = focal_freq >= min_freq

    focal_alt_allele_present_cond1 = focal_depth >= min_depth
    focal_alt_allele_present_cond2 = focal_freq <= (1 - min_freq)

    # ref allele present in background
    background_freq = np.array(background_freq)
    background_depth = np.array(background_depth)
    background_ref_allele_present_cond1 = background_depth >= min_depth
    background_ref_allele_present_cond2 = background_freq >= min_freq

    background_alt_allele_present_cond1 = background_depth >= min_depth
    background_alt_allele_present_cond2 = background_freq <= (1 - min_freq)

    num_ref_valid_background = sum(
        background_ref_allele_present_cond1 & background_ref_allele_present_cond2
    )
    num_alt_valid_background = sum(
        background_alt_allele_present_cond1 & background_alt_allele_present_cond2
    )

    # THere will be no marker
    if (num_ref_valid_background > 0) and (num_alt_valid_background > 0):
        return [], []
    if num_ref_valid_background == 0:
        focal_ref_conditions = (
            focal_ref_allele_present_cond1 & focal_ref_allele_present_cond2
        )
        samples_with_marker_allele_ref = [
            focal_sample_names[i] for i, x in enumerate(focal_ref_conditions) if x
        ]

    if num_alt_valid_background == 0:
        focal_alt_conditions = (
            focal_alt_allele_present_cond1 & focal_alt_allele_present_cond2
        )
        samples_with_marker_allele_alt = [
            focal_sample_names[i] for i, x in enumerate(focal_alt_conditions) if x
        ]
    return samples_with_marker_allele_ref, samples_with_marker_allele_alt


def check_marker_allele_pairs(
    pair_dict_infant_key,
    focal_samples,
    focal_depth,
    focal_freq,
    background_samples,
    background_depth,
    background_freq,
):
    # filter uptop
    min_depth = 30
    min_freq = 0.10

    samples_with_marker_allele_ref = []
    samples_with_marker_allele_alt = []

    for infant_key in pair_dict_infant_key.keys():
        infant_key_depth = focal_depth[focal_samples.index(infant_key)]
        if infant_key_depth >= min_depth:
            for mother_member in pair_dict_infant_key[infant_key]:
                mother_member_depth = background_depth[
                    background_samples.index(mother_member)
                ]
                if mother_member_depth >= min_depth:

                    mother_member_freq = background_freq[
                        background_samples.index(mother_member)
                    ]
                    infant_key_freq = focal_freq[focal_samples.index(infant_key)]

                    background_freq_fs = np.array(
                        [
                            background_freq[bs]
                            for bs in range(len(background_samples))
                            if background_samples[bs] != mother_member
                        ]
                    )
                    background_depth_fs = np.array(
                        [
                            background_depth[bs]
                            for bs in range(len(background_samples))
                            if background_samples[bs] != mother_member
                        ]
                    )

                    (
                        samples_with_marker_allele_ref_fs,
                        samples_with_marker_allele_alt_fs,
                    ) = process_allele_data(
                        focal_sample_names=[infant_key, mother_member],
                        focal_depth=[infant_key_depth, mother_member_depth],
                        focal_freq=[infant_key_freq, mother_member_freq],
                        min_depth=min_depth,
                        min_freq=min_freq,
                        background_freq=background_freq_fs,
                        background_depth=background_depth_fs,
                    )

                    if len(samples_with_marker_allele_ref_fs) == 2:

                        samples_with_marker_allele_ref.append(
                            infant_key + "_" + mother_member
                        )
                    if len(samples_with_marker_allele_alt_fs) == 2:
                        samples_with_marker_allele_alt.append(
                            infant_key + "_" + mother_member
                        )

    return list(set(samples_with_marker_allele_ref + samples_with_marker_allele_alt))


def check_marker_allele(
    focal_sample_names,
    focal_depth,
    focal_freq,
    background_samples,
    background_depth,
    background_freq,
):

    # filter uptop
    min_depth = 30
    min_freq = 0.10

    check_if_focal_in_background = list(
        set(focal_sample_names).intersection(set(background_samples))
    )
    if len(check_if_focal_in_background) == 0:

        (
            samples_with_marker_allele_ref,
            samples_with_marker_allele_alt,
        ) = process_allele_data(
            focal_sample_names=focal_sample_names,
            focal_depth=focal_depth,
            focal_freq=focal_freq,
            min_depth=min_depth,
            min_freq=min_freq,
            background_freq=background_freq,
            background_depth=background_depth,
        )

    else:
        samples_with_marker_allele_ref = []
        samples_with_marker_allele_alt = []
        for fs in range(len(focal_sample_names)):
            if focal_depth[fs] < min_depth:
                continue
            background_freq_fs = np.array(
                [
                    background_freq[bs]
                    for bs in range(len(background_samples))
                    if background_samples[bs] != focal_sample_names[fs]
                ]
            )
            background_depth_fs = np.array(
                [
                    background_depth[bs]
                    for bs in range(len(background_samples))
                    if background_samples[bs] != focal_sample_names[fs]
                ]
            )

            (
                samples_with_marker_allele_ref_fs,
                samples_with_marker_allele_alt_fs,
            ) = process_allele_data(
                focal_sample_names=[focal_sample_names[fs]],
                focal_depth=[focal_depth[fs]],
                focal_freq=[focal_freq[fs]],
                min_depth=min_depth,
                min_freq=min_freq,
                background_freq=background_freq_fs,
                background_depth=background_depth_fs,
            )
            samples_with_marker_allele_ref.extend(samples_with_marker_allele_ref_fs)
            samples_with_marker_allele_alt.extend(samples_with_marker_allele_alt_fs)
    return list(set(samples_with_marker_allele_ref + samples_with_marker_allele_alt))


"""
Columns 1 is for for simulation (infant),
Columsn 2 fo backhed (mothers)
"""


def read_two_files(
    depth_path_1,
    freq_path_1,
    depth_path_2,
    freq_path_2,
    infant_samples,
    mother_samples,
    infant_mother_pairs,
    background_samples,
):
    min_freq = 0.10
    min_depth = 30
    min_allele_support = 3

    infant_marker_allele_dict = dict()
    infant_not_marker_allele_dict = dict()

    mother_marker_allele_dict = dict()
    mother_not_marker_allele_dict = dict()

    pair_marker_allele_dict = dict()
    pair_not_marker_allele_dict = dict()

    pair_dict_infant_key = dict()
    pair_dict_mother_key = dict()

    for pair in infant_mother_pairs:
        if pair[0] not in pair_dict_infant_key.keys():
            pair_dict_infant_key[pair[0]] = []
        pair_dict_infant_key[pair[0]].append(pair[1])
        if pair[1] not in pair_dict_mother_key.keys():
            pair_dict_mother_key[pair[1]] = []
        pair_dict_mother_key[pair[1]].append(pair[0])

    columns1 = get_column_names(depth_path_1)
    columns2 = get_column_names(depth_path_2)
    columns1 = columns1[1 : len(columns1)]
    columns2 = columns2[1 : len(columns2)]

    infant_indices = [columns1.index(s1) for s1 in infant_samples]
    mother_indices = [columns2.index(s2) for s2 in mother_samples]
    background_indices = [columns2.index(s3) for s3 in background_samples]

    with bzopen(depth_path_1, "r") as bzdepth1, bzopen(
        freq_path_1, "r"
    ) as bzfreq1, bzopen(depth_path_2, "r") as bzdepth2, bzopen(
        freq_path_2, "r"
    ) as bzfreq2:
        for i, (bzdepth1_line, bzfreq1_line, bzdepth2_line, bzfreq2_line) in enumerate(
            zip(bzdepth1, bzfreq1, bzdepth2, bzfreq2)
        ):
            if i == 0:
                continue
            # if i < 31600:
            #     continue

            if (i % 1000) == 0:
                print(i)
            # if i > 31650:
            #     return (
            #         infant_marker_allele_dict,
            #         mother_marker_allele_dict,
            #         pair_marker_allele_dict,
            #     )

            snp_name = bzdepth1_line.decode().rstrip().split("\t")[0]

            ## Preprocessing

            decoded_depth_1 = bzdepth1_line.decode().rstrip().split("\t")
            decoded_freq_1 = bzfreq1_line.decode().rstrip().split("\t")
            decoded_depth_2 = bzdepth2_line.decode().rstrip().split("\t")
            decoded_freq_2 = bzfreq2_line.decode().rstrip().split("\t")

            depth_vec_1 = [int(x) for x in decoded_depth_1[1 : len(decoded_depth_1)]]
            depth_vec_2 = [int(x) for x in decoded_depth_2[1 : len(decoded_depth_2)]]
            freq_vec_1 = [float(x) for x in decoded_freq_1[1 : len(decoded_freq_1)]]
            freq_vec_2 = [float(x) for x in decoded_freq_2[1 : len(decoded_freq_2)]]

            ## Aligned to my samples
            infant_depth = [depth_vec_1[s_i] for s_i in infant_indices]
            mother_depth = [depth_vec_2[s_i] for s_i in mother_indices]
            infant_freq = [freq_vec_1[s_i] for s_i in infant_indices]
            mother_freq = [freq_vec_2[s_i] for s_i in mother_indices]
            background_depth = [depth_vec_2[s_i] for s_i in background_indices]
            background_freq = [freq_vec_2[s_i] for s_i in background_indices]

            num_valid_infants = (np.array(infant_depth) >= min_depth).sum()
            num_valid_mothers = (np.array(mother_depth) >= min_depth).sum()

            ## niether infant or mother has
            if (num_valid_infants == 0) and (num_valid_mothers == 0):
                continue
            ## infant only has

            # print(
            #     "num_valid_infants "
            #     + str(num_valid_infants)
            #     + " num_valid_mothers "
            #     + str(num_valid_mothers)
            # )
            if num_valid_infants > 0:

                samples_with_marker_allele = check_marker_allele(
                    focal_sample_names=infant_samples,
                    focal_depth=infant_depth,
                    focal_freq=infant_freq,
                    background_samples=background_samples,
                    background_depth=background_depth,
                    background_freq=background_freq,
                )
                if len(samples_with_marker_allele) > 0:

                    for sample in samples_with_marker_allele:
                        if sample not in infant_marker_allele_dict.keys():
                            infant_marker_allele_dict[sample] = []
                        infant_marker_allele_dict[sample].append(snp_name)
            if num_valid_mothers > 0:

                samples_with_marker_allele = check_marker_allele(
                    focal_sample_names=mother_samples,
                    focal_depth=mother_depth,
                    focal_freq=mother_freq,
                    background_samples=background_samples,
                    background_depth=background_depth,
                    background_freq=background_freq,
                )
                if len(samples_with_marker_allele) > 0:

                    for sample in samples_with_marker_allele:
                        if sample not in mother_marker_allele_dict.keys():
                            mother_marker_allele_dict[sample] = []
                        mother_marker_allele_dict[sample].append(snp_name)

                # print("before pair")
            if (num_valid_mothers) > 0 and (num_valid_infants > 0):

                pairs_with_marker_allele = check_marker_allele_pairs(
                    pair_dict_infant_key=pair_dict_infant_key,
                    focal_samples=infant_samples,
                    focal_depth=infant_depth,
                    focal_freq=infant_freq,
                    background_samples=background_samples,
                    background_depth=background_depth,
                    background_freq=background_freq,
                )
                if len(pairs_with_marker_allele) > 0:
                    for pair in pairs_with_marker_allele:
                        if pair not in pair_marker_allele_dict.keys():
                            pair_marker_allele_dict[pair] = []
                        pair_marker_allele_dict[pair].append(snp_name)
                # print("after pair")

            # print("done")
    return (
        infant_marker_allele_dict,
        mother_marker_allele_dict,
        pair_marker_allele_dict,
    )

