import os
import pandas as pd
import numpy as np
from scipy.io import savemat


def process_cluster_folder(csv_folder_path, output_mat_file):
    """
    Process CSV files from a resting state cluster folder and create a .mat file

    Args:
        csv_folder_path: Path to folder containing CSV files
        output_mat_file: Path to output .mat file

    Returns:
        True if successful, False otherwise
    """
    if not os.path.exists(csv_folder_path):
        print(f"Error: Folder '{csv_folder_path}' does not exist.")
        return False

    print(f"\nProcessing folder: {csv_folder_path}")
    print(f"Output file: {output_mat_file}")

    # Get list of CSV files
    csv_files = [f for f in os.listdir(csv_folder_path) if f.endswith('.csv')]
    num_files = len(csv_files)

    if num_files == 0:
        print(f"Error: No CSV files found in {csv_folder_path}")
        return False

    print(f"Found {num_files} CSV files")

    # Initialize lists for the structure
    data_list = []
    participant_ids = []

    # Process each CSV file
    for i, csv_file in enumerate(csv_files, 1):
        csv_path = os.path.join(csv_folder_path, csv_file)

        print(f"Processing file {i}/{num_files}: {csv_file}")

        try:
            # Read the CSV file
            # Each row contains tab-separated values for all regions at one timepoint
            with open(csv_path, 'r') as f:
                lines = f.readlines()

            print(f"  Found {len(lines)} timepoints (rows)")

            # Parse each line to extract region values
            timepoints_data = []
            for line_idx, line in enumerate(lines):
                line = line.strip()
                if line:  # Skip empty lines
                    # Split by tabs to get values for each region at this timepoint
                    values = line.split('\t')

                    # Convert to float
                    try:
                        float_values = [float(v) for v in values]
                        timepoints_data.append(float_values)
                    except ValueError as e:
                        print(f"  Warning: Could not convert line {line_idx+1} to numbers: {e}")
                        continue

            # Convert to numpy array: timepoints × regions
            timepoints_array = np.array(timepoints_data)
            num_timepoints = timepoints_array.shape[0]
            num_regions = timepoints_array.shape[1]
            print(f"  Parsed data shape: {num_timepoints} timepoints × {num_regions} regions")

            # Transpose to get regions × timepoints (as expected in neuroscience)
            data = timepoints_array.T
            print(f"  Final data shape: {data.shape[0]} regions × {data.shape[1]} timepoints")

            # Store the data
            data_list.append(data)

            # Extract subject key from filename
            # Format: NDAR_INVAG023WG3_parcellated_task-hammerAP_run-01_bold_Atlas_MSMAll_hp2000_clean.csv
            # Subject key is everything before the second underscore (NDAR_INVAG023WG3)
            underscore_positions = [pos for pos, char in enumerate(csv_file) if char == '_']

            if len(underscore_positions) >= 2:
                subject_key = csv_file[:underscore_positions[1]]
            else:
                # Fallback: use filename without extension
                subject_key = os.path.splitext(csv_file)[0]
                print(f"  Warning: Could not extract subject key from {csv_file}, using full filename")

            participant_ids.append(subject_key)
            print(f"  Extracted subject key: {subject_key}\n")

        except Exception as e:
            print(f"Error processing {csv_file}: {str(e)}")
            continue

    if not data_list:
        print(f"Error: No valid CSV files could be processed in {csv_folder_path}")
        return False

    # Create the structure similar to MATLAB format
    # Create a structured array that mimics the MATLAB cell array
    resting_state_data = np.empty((len(data_list), 2), dtype=object)

    for i in range(len(data_list)):
        resting_state_data[i, 0] = data_list[i]
        resting_state_data[i, 1] = participant_ids[i]

    # Save to .mat file
    # Use the filename (without extension) as the variable name
    var_name = os.path.splitext(os.path.basename(output_mat_file))[0]
    try:
        savemat(output_mat_file, {var_name: resting_state_data})
        print(f"\nConversion complete! Saved to: {output_mat_file}")
        print(f"Structure size: {len(data_list)} participants x 2 columns")

        # Display first few entries for verification
        print("\nFirst few entries:")
        for i in range(min(3, len(data_list))):
            print(f"  Row {i+1}: Data size {data_list[i].shape[0]}x{data_list[i].shape[1]}, Subject ID: {participant_ids[i]}")

        print(f"\nExpected format: Each matrix should be regions x timepoints")
        return True

    except Exception as e:
        print(f"Error saving file: {str(e)}")
        return False


def process_cluster_folder_with_count(csv_folder_path, output_mat_file):
    """Wrapper around process_cluster_folder that also returns the participant count."""
    csv_files = [f for f in os.listdir(csv_folder_path) if f.endswith('.csv')] if os.path.exists(csv_folder_path) else []
    success = process_cluster_folder(csv_folder_path, output_mat_file)
    return success, len(csv_files) if success else 0


def main():
    """
    Main function to process both cluster folders and create anhedonia group .mat files
    """
    print("=" * 80)
    print("Resting State fMRI Data Converter")
    print("Converting cluster data to anhedonia group .mat files")
    print("=" * 80)

    # Define the base path where the cluster folders are located
    cluster_base_path = "/Users/proghani/Documents/personal/code_experiments/neuro/phd_stuff/tcp_parcellations/data/anhedonia/NEW_4_factor_clustering_ipnybV4"

    # Define input folders (cluster_0 = high anhedonia, cluster_1 = low anhedonia)
    cluster_0_folder = os.path.join(cluster_base_path, "cluster_0")
    cluster_1_folder = os.path.join(cluster_base_path, "cluster_1")

    # Output files saved inside their respective cluster folders
    output_high_anhedonia = os.path.join(cluster_0_folder, "high_anhedonia.mat")
    output_low_anhedonia = os.path.join(cluster_1_folder, "low_anhedonia.mat")

    # Process cluster_0 (high anhedonia)
    print("\n" + "=" * 80)
    print("Processing Cluster 0 (High Anhedonia)")
    print("=" * 80)
    success_high, count_high = process_cluster_folder_with_count(cluster_0_folder, output_high_anhedonia)

    # Process cluster_1 (low anhedonia)
    print("\n" + "=" * 80)
    print("Processing Cluster 1 (Low Anhedonia)")
    print("=" * 80)
    success_low, count_low = process_cluster_folder_with_count(cluster_1_folder, output_low_anhedonia)

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    if success_high:
        print(f"✓ High anhedonia data saved to: {output_high_anhedonia}")
    else:
        print(f"✗ Failed to process high anhedonia data")

    if success_low:
        print(f"✓ Low anhedonia data saved to: {output_low_anhedonia}")
    else:
        print(f"✗ Failed to process low anhedonia data")

    print(f"\nParticipant counts:")
    print(f"  High anhedonia (cluster_0): {count_high}")
    print(f"  Low anhedonia  (cluster_1): {count_low}")
    print(f"  Total: {count_high + count_low}")

    if success_low and success_high:
        print("\n✓ All conversions completed successfully!")
    else:
        print("\n⚠ Some conversions failed. Please check the errors above.")


if __name__ == "__main__":
    main()
