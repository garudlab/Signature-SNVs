import numpy as np
from scipy.optimize import minimize
from itertools import compress
from bz2 import BZ2File as bzopen


def LL1_fn(n, m, gamma_i):
    alpha_i = 1
    P = alpha_i * gamma_i
    L = np.dot(n, np.log(P)) + np.dot(m, np.log(1 - P))
    return 1 * L


def LL2_fn(alpha_sources, *args):
    sink_n_param, sink_m_param, gamma_mat_param = args[0], args[1], args[2]
    return 1 * (
        np.dot(sink_n_param, np.log(np.dot(gamma_mat_param, alpha_sources)))
        + np.dot(sink_m_param, np.log(np.dot(1 - gamma_mat_param, alpha_sources)))
    )


def LL2_fn_post(alpha_sources, args):
    sink_n_param, sink_m_param, gamma_mat_param = args[0], args[1], args[2]
    return 1 * (
        np.dot(sink_n_param, np.log(np.dot(gamma_mat_param, alpha_sources)))
        + np.dot(sink_m_param, np.log(np.dot(1 - gamma_mat_param, alpha_sources)))
    )


def score_snps(
    filtered_freq_matrix_infant,
    filtered_counts_matrix_infant,
    filtered_freq_matrix_moms,
    filtered_counts_matrix_moms,
):

    ref_frequencies = np.concatenate(
        ((1 - filtered_freq_matrix_infant), (1 - filtered_freq_matrix_moms)), axis=1
    )  # for source
    counts = np.concatenate(
        (filtered_counts_matrix_infant, filtered_counts_matrix_moms), axis=1
    )

    alt_counts = np.round_(np.multiply((1 - ref_frequencies), counts), decimals=0) + 1
    ref_counts = np.round_(np.multiply(ref_frequencies, counts), decimals=0) + 1
    adjusted_frequencies = np.divide(ref_counts, (alt_counts + ref_counts))

    sink_n = ref_counts[:, 0]
    sink_m = alt_counts[:, 0]
    num_samples = counts.shape[1]
    max_likelihoods = []
    num_fail = 0
    for snp_id in range(adjusted_frequencies.shape[0]):

        L_diff = []

        for i in range(1, num_samples):  # each sample as left at some point (see ref A)

            # print("ll1 " + str(i))
            ll1 = LL1_fn(
                sink_n[snp_id], sink_m[snp_id], adjusted_frequencies[snp_id, i]
            )
            # print(ll1)

            if num_samples > 2:

                selector = [
                    x
                    for x in range(adjusted_frequencies.shape[1])
                    if (x != i and x != 0)
                ]  # ref A": i is the sample we are leavin gout

                gamma_mat = adjusted_frequencies[snp_id, selector]

                # print("solve for alpha")
                x0 = np.repeat(float(1) / (len(gamma_mat)), len(gamma_mat))
                # print("X0")
                # print(x0)
                cons = {"type": "eq", "fun": lambda x: 1 - sum(x)}
                bnds = tuple((0, 1) for x in x0)
                # additional = {'sink_n_param':sink_n,'sink_m_param':sink_m}
                additional = (sink_n[snp_id], sink_m[snp_id], gamma_mat)

                # res = minimize(LL2_fn, x0=x0, args =additional,method='SLSQP',
                #              constraints=cons,
                #              bounds=bnds)
                try:
                    res = minimize(
                        LL2_fn,
                        x0=x0,
                        args=additional,
                        method="SLSQP",
                        constraints=cons,
                        bounds=bnds,
                    )
                    # print("suceed")
                except:
                    num_fail += 1
                    # print("Fail")
                    continue

                # print(ll1)
                # print("ll2")
                ll2 = LL2_fn_post(res.x, args=additional)
                # print(ll2)
                L_diff.append(ll1 - ll2)
            else:
                L_diff.append(ll1)

            # print("multiple time")
            # print(np.dot(gamma_mat,  made_up_alpha))

        max_diff = np.max(L_diff)
        max_likelihoods.append(max_diff)
    print("number failed")
    print(num_fail)
    # print(max_likelihoods)

    ## FAVORING SNPS UNIQUE TO INFANT WITH mean depth) and assign highest likelihood in bracket =-> I changed tjhis because it msessed with the SD calc
    zero_freq_alt = filtered_freq_matrix_moms.sum(axis=1) == 0
    mean_depth_infant = np.mean(filtered_counts_matrix_infant[:, 0])

    # unique_infant_filter = np.array(zero_freq_alt) & np.array(filtered_counts_matrix_infant[:,0] > mean_depth_infant)  ## if only want snps with > mean read deptg
    unique_infant_filter = np.array(zero_freq_alt)  # all infant-unqieu snps
    # print("number of zero freq alt")
    # print(sum(zero_freq_alt))
    # print("likelihoods for the snv")
    # zero_likelihoods = []
    # regular_likelihoods = []
    # for i in range(len(max_likelihoods)):
    #   if unique_infant_filter[i]:
    #     zero_likelihoods.append(max_likelihoods[i])
    #   else:
    #     regular_likelihoods.append(max_likelihoods[i])

    # #temp_max_likelihoods = [max_likelihoods[i] if unique_infant_filter[i] else -33 for i in range(len(max_likelihoods)) ]
    # print("zero_likelihoods")
    # if len(zero_likelihoods) > 0:
    #   print(np.percentile(zero_likelihoods, [25, 50 , 75]))
    # print("regular_likelihoods")
    # if len(regular_likelihoods) > 0:
    #   print(np.percentile(regular_likelihoods, [25, 50 , 75]))

    max_likelihoods = [
        -1 if unique_infant_filter[i] else max_likelihoods[i]
        for i in range(len(max_likelihoods))
    ]
    return max_likelihoods


def prevalence_filter(
    freq_matrix,
    freq_rownames,
    freq_colnames,
    depth_matrix,
    depth_rownames,
    depth_colnames,
    species,
    sink_sample_names,
):
    freq_rownames = [species + "|" + i for i in freq_rownames]

    # We don't care about Sites that are not in sink at all
    sink_bit_mask = [s in sink_sample_names for s in freq_colnames]
    freq_matrix_sink = freq_matrix[:, sink_bit_mask]
    print("freq_matrix_sink.shape")
    print(freq_matrix_sink.shape)
    snp_presence = freq_matrix_sink > 0
    snp_prevalence = snp_presence.sum(axis=1) / snp_presence.shape[1]

    prevalence_filter = np.array(snp_prevalence > 0)
    print("sum filter")
    print(sum(prevalence_filter))
    print(prevalence_filter)

    freq_rownames = np.array(freq_rownames)[prevalence_filter]
    depth_rownames = np.array(depth_rownames)[prevalence_filter]

    freq_matrix = freq_matrix[prevalence_filter, :]
    depth_matrix = depth_matrix[prevalence_filter, :]

    return (freq_matrix, freq_rownames, depth_matrix, depth_rownames)


def make_np_array_from_file(
    path,
    start_index,
    end_index,
    all_samples=[],
    custom_delimiter="\t",
    no_headers=False,
    custom_features=[],
):
    # print(no_headers)
    print(path)
    all_samples = list(all_samples)
    # print("all_samples")
    # print(sorted(all_samples))
    # print(type(all_samples))
    # print("LEN " + str(len(all_samples)))

    # end index never implies it's the final SNP matrix
    with bzopen(path, "r") as bzfin:
        """Handle lines here"""
        lines = []
        rownames = []

        for i, line in enumerate(bzfin):
            if i % 10000 == 0:
                print(i)
            if i == end_index:
                break
            if not no_headers:
                if i == 0:

                    colnames = line.decode().rstrip().split(custom_delimiter)
                    # print("colnames")
                    # print(sorted(colnames))
                    # print(type(colnames))
                    # print("LEN " + str(len(colnames)))

                    samplewide_bit_mask = [s in all_samples for s in colnames]

                    # print("INtersection")
                    # print(set(colnames).intersection(set(all_samples)))
                    samplewide_bit_mask[0] = True
                    # print(samplewide_bit_mask)
                    # print("length before filtering out")
                    # print(len(colnames))
                    # print(samplewide_bit_mask)
                    colnames = list(compress(colnames, samplewide_bit_mask))
                    # print(len(newc))

            if i >= start_index:
                # print("CHECK")
                # print(line.decode().rstrip())

                processed_line = line.decode().rstrip().split(custom_delimiter)

                processed_line = list(compress(processed_line, samplewide_bit_mask))
                snp_array = []
                if no_headers:
                    for x in processed_line:
                        try:
                            snp_array.append(float(x))
                        except ValueError:
                            snp_array.append(np.nan)
                            pass
                    lines.append(snp_array)
                else:
                    # print("processed_line")
                    # print(processed_line)
                    if len(custom_features) > 0:
                        if processed_line[0] in custom_features:
                            # print("FOUND")
                            rownames.append(processed_line[0])
                            for x in processed_line[1:]:
                                try:
                                    snp_array.append(float(x))
                                except ValueError:
                                    snp_array.append(np.nan)
                                    pass
                            # print( snp_array)
                            lines.append(snp_array)

                    else:
                        # print(processed_line)

                        rownames.append(processed_line[0])
                        for x in processed_line[1:]:
                            try:
                                snp_array.append(float(x))
                            except ValueError:
                                snp_array.append(np.nan)
                                pass

                        lines.append(snp_array)
    # print(lines)

    snp_matrix = np.array(lines)
    print("dim snp_matrix")
    print(snp_matrix.shape)
    if not no_headers:
        row_snp_names = rownames
        col_sample_names = colnames[1:]
        return (snp_matrix, row_snp_names, col_sample_names)
    else:
        return snp_matrix
