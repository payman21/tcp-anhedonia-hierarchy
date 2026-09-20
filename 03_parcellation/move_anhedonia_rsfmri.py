import os
import shutil
import pandas as pd
from pathlib import Path

import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from config import TCP_DATA

# Use absolute path to phenotype file
cluster_df = pd.read_csv(f'{TCP_DATA}/anhedonia/NEW_4_factor_clustering_ipnybV4/anhedonia_polychoric_gmm_clusters.csv')

# Define directories
source_dir = Path(f'{TCP_DATA}/parcellated_rsfmri')
cluster_0_dir = Path(f'{TCP_DATA}/anhedonia/NEW_4_factor_clustering_ipnybV4/cluster_0')
cluster_1_dir = Path(f'{TCP_DATA}/anhedonia/NEW_4_factor_clustering_ipnybV4/cluster_1')

# Create output directories if they don't exist
cluster_0_dir.mkdir(parents=True, exist_ok=True)
cluster_1_dir.mkdir(parents=True, exist_ok=True)

# Track processing - now tracking missing files by cluster too
moved_files = {'cluster_0': [], 'cluster_1': [], 'not_found_cluster_0': [], 'not_found_cluster_1': [], 'no_cluster': []}

# Process each subject in phenotype file
for _, row in cluster_df.iterrows():
    subject_id = row['subjectkey']
    cluster = row['gmm_cluster']
    
    # Build expected file path
    subject_folder = source_dir / subject_id
    expected_filename = f"{subject_id}_concatenated_parcellated_task-rest_bold_Atlas_MSMAll_hp2000_clean.csv"
    source_file = subject_folder / expected_filename
    
    # Check if file exists
    if source_file.exists():
        # Determine destination based on cluster
        if cluster == 0:
            dest_file = cluster_0_dir / expected_filename
            shutil.copy2(source_file, dest_file)
            moved_files['cluster_0'].append(subject_id)
        elif cluster == 1:
            dest_file = cluster_1_dir / expected_filename
            shutil.copy2(source_file, dest_file)
            moved_files['cluster_1'].append(subject_id)
        else:
            moved_files['no_cluster'].append(subject_id)
    else:
        # Track missing files by cluster
        if cluster == 0:
            moved_files['not_found_cluster_0'].append(subject_id)
        elif cluster == 1:
            moved_files['not_found_cluster_1'].append(subject_id)
        else:
            moved_files['no_cluster'].append(subject_id)

# Report results
print(f"Files copied to cluster_0: {len(moved_files['cluster_0'])}")
print(f"Files copied to cluster_1: {len(moved_files['cluster_1'])}")
print(f"Total files not found: {len(moved_files['not_found_cluster_0']) + len(moved_files['not_found_cluster_1'])}")

# Print missing subjects by cluster
if moved_files['not_found_cluster_0']:
    print(f"\nMissing subjects from cluster 0: {len(moved_files['not_found_cluster_0'])} subjects")
    print(f"Subject IDs: {moved_files['not_found_cluster_0']}")

if moved_files['not_found_cluster_1']:
    print(f"\nMissing subjects from cluster 1: {len(moved_files['not_found_cluster_1'])} subjects")
    print(f"Subject IDs: {moved_files['not_found_cluster_1']}")

# Verify balance
print(f"\nCluster distribution in neuroimaging data:")
print(f"Cluster 0 (high anhedonia): {len(moved_files['cluster_0'])} subjects")
print(f"Cluster 1 (low anhedonia): {len(moved_files['cluster_1'])} subjects")

# Save processing log
processing_log = pd.DataFrame({
    'subjectkey': moved_files['cluster_0'] + moved_files['cluster_1'],
    'cluster': [0]*len(moved_files['cluster_0']) + [1]*len(moved_files['cluster_1']),
    'status': ['copied']*len(moved_files['cluster_0'] + moved_files['cluster_1'])
})
