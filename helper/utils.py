import pickle

def add_cluster_from_pickle(andata, pkl_file, new_col, match_col="cell_id"):
    """
    Add clustering results from a pickle dictionary to andata.obs.

    The pickle file should contain:
        key   = original obs name
        value = cluster label
    """

    with open(pkl_file, "rb") as f:
        cluster_dict = pickle.load(f)

    # make sure dictionary keys are strings
    cluster_dict = {str(k): str(v) for k, v in cluster_dict.items()}

    # map cluster labels into AnnData
    andata.obs[new_col] = andata.obs[match_col].map(cluster_dict)

    # report mapping quality
    n_total = andata.n_obs
    n_missing = andata.obs[new_col].isna().sum()

    print(f"{new_col}")
    print(f"  total cells:   {n_total}")
    print(f"  matched cells: {n_total - n_missing}")
    print(f"  missing cells: {n_missing}")

    # convert to category for plotting
    andata.obs[new_col] = andata.obs[new_col].astype("category")

    return andata