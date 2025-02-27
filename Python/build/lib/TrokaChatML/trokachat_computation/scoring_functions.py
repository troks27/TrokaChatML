def calculate_scores_parallel(args):
    row, deg_data, cluster_percent, sample_name, counts_filtered, general_filepath = args

    def calculate_communication_score(row, deg_data, source_pct_total, target_pct_total, source_cellcount, target_cellcount, total_counts_sample):
        ligand = row['Ligand']
        subunits = row['Subunits']
        ligand_data = deg_data[(deg_data['gene'] == ligand) & (deg_data['cluster'] == row['Source'])]
        if not ligand_data.empty:
            ligand_score = (ligand_data['avg_log2FC'].iloc[0] * ((ligand_data['pct.1'].iloc[0] * source_cellcount) / total_counts_sample))
        else:
            return 0
        subunit_scores = []
        for subunit in subunits:
            subunit_data = deg_data[(deg_data['gene'] == subunit) & (deg_data['cluster'] == row['Target'])]
            if not subunit_data.empty:
                subunit_score = (subunit_data['avg_log2FC'].values[0] * ((subunit_data['pct.1'].values[0] * target_cellcount) / total_counts_sample))
                subunit_scores.append(subunit_score)
        average_subunit_score = sum(subunit_scores) / len(subunit_scores)
        communication_score = ligand_score * average_subunit_score
        return communication_score

    def calculate_uniqueness_score(row, deg_data):
        ligand = row['Ligand']
        subunits = row['Subunits']
        ligand_data = deg_data[(deg_data['gene'] == ligand) & (deg_data['cluster'] == row['Source'])]
        ligand_score = ligand_data['pct.1'].values[0] / ligand_data['pct.2'].values[0]
        subunit_scores = []
        for subunit in subunits:
            subunit_data = deg_data[(deg_data['gene'] == subunit) & (deg_data['cluster'] == row['Target'])]
            subunit_score = subunit_data['pct.1'].values[0] / subunit_data['pct.2'].values[0]
            subunit_scores.append(subunit_score)
        average_subunit_score = sum(subunit_scores) / len(subunit_scores)
        uniqueness_score = ligand_score * average_subunit_score
        return uniqueness_score

    source_cluster = str(int(float(row['Source'])))
    target_cluster = str(int(float(row['Target'])))
    
    # Debug: Print the keys of cluster_percent and the sample, source, and target clusters
    #print(f"Calculating scores for row: {row}, sample name: {sample_name}")
    #print(f"Cluster percent keys: {cluster_percent.keys()}")
    #print(f"Attempting to access cluster_percent for sample_name: {sample_name}, Source: {source_cluster}, Target: {target_cluster}")
    
    try:
        source_pct_total = cluster_percent.loc[sample_name, source_cluster]
        target_pct_total = cluster_percent.loc[sample_name, target_cluster]
    except KeyError as e:
        print(f"KeyError: {e}")
        print(f"Available keys in cluster_percent: {cluster_percent.index.tolist()}")
        raise

    source_cellcount = counts_filtered.loc[sample_name, source_cluster]
    target_cellcount = counts_filtered.loc[sample_name, target_cluster]
    total_counts_sample = counts_filtered.loc[sample_name, :].sum()
    communication_score = calculate_communication_score(row, deg_data, source_pct_total, target_pct_total, source_cellcount, target_cellcount, total_counts_sample)
    uniqueness_score = calculate_uniqueness_score(row, deg_data)
    return communication_score, uniqueness_score