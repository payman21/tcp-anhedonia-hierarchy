import os
import pandas as pd
import numpy as np
from scipy.io import savemat

def convert_csv_to_mat():
    """
    Convert CSV files to .mat structure with signal data and participant IDs
    CSV structure: Each row = one timepoint with tab-separated values for all regions
    Expected result: 232 regions × timepoints matrix
    """
    
    # Get input folder path from user via terminal
    print("Copy and paste the path to the folder containing your CSV files:")
    csv_folder_path = input().strip()
    
    # Remove quotes if user copied path with quotes
    csv_folder_path = csv_folder_path.strip('"\'')
    
    if not os.path.exists(csv_folder_path):
        print(f"Error: Folder '{csv_folder_path}' does not exist.")
        return
    
    # Get output file path from user via terminal
    print("\nCopy and paste the path where you want to save the .mat file:")
    print("(Include the filename, e.g., /Users/username/Desktop/my_data.mat)")
    output_mat_file = input().strip()
    
    # Remove quotes if user copied path with quotes
    output_mat_file = output_mat_file.strip('"\'')
    
    # Add .mat extension if not provided
    if not output_mat_file.endswith('.mat'):
        output_mat_file += '.mat'
    
    print(f"\nInput folder: {csv_folder_path}")
    print(f"Output file: {output_mat_file}")
    
    # Get list of CSV files
    csv_files = [f for f in os.listdir(csv_folder_path) if f.endswith('.csv')]
    num_files = len(csv_files)
    
    if num_files == 0:
        print("Error: No CSV files found in the specified folder")
        return
    
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
            
            print(f"Found {len(lines)} timepoints (rows)")
            
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
                        print(f"Warning: Could not convert line {line_idx+1} to numbers: {e}")
                        continue
            
            # Convert to numpy array: timepoints × regions
            timepoints_array = np.array(timepoints_data)
            num_timepoints = timepoints_array.shape[0]
            num_regions = timepoints_array.shape[1]
            print(f"Parsed data shape: {num_timepoints} timepoints × {num_regions} regions")
            
            # Transpose to get regions × timepoints (as expected in neuroscience)
            data = timepoints_array.T
            print(f"Final data shape: {data.shape[0]} regions × {data.shape[1]} timepoints")
            
            # Store the data
            data_list.append(data)
            
            # Extract subject key from filename
            # Format: NDAR_INVAG023WG3_concatenated_parcellated_task-rest_bold_Atlas_MSMAll_hp2000_clean.csv
            # Subject key is everything before the second underscore (NDAR_INVAG023WG3)
            underscore_positions = [pos for pos, char in enumerate(csv_file) if char == '_']
            
            if len(underscore_positions) >= 2:
                subject_key = csv_file[:underscore_positions[1]]
            else:
                # Fallback: use filename without extension
                subject_key = os.path.splitext(csv_file)[0]
                print(f"Warning: Could not extract subject key from {csv_file}, using full filename")
            
            participant_ids.append(subject_key)
            print(f"Extracted subject key: {subject_key}\n")
            
        except Exception as e:
            print(f"Error processing {csv_file}: {str(e)}")
            continue
    
    if not data_list:
        print("Error: No valid CSV files could be processed")
        return
    
    # Create the structure similar to MATLAB format
    # Create a structured array that mimics the MATLAB cell array
    psilodep2_before = np.empty((len(data_list), 2), dtype=object)
    
    for i in range(len(data_list)):
        psilodep2_before[i, 0] = data_list[i]
        psilodep2_before[i, 1] = participant_ids[i]
    
    # Save to .mat file
    # Use the filename (without extension) as the variable name
    var_name = os.path.splitext(os.path.basename(output_mat_file))[0]
    try:
        savemat(output_mat_file, {var_name: psilodep2_before})
        print(f"\nConversion complete! Saved to: {output_mat_file}")
        print(f"Structure size: {len(data_list)} participants x 2 columns")
        
        # Display first few entries for verification
        print("\nFirst few entries:")
        for i in range(min(3, len(data_list))):
            print(f"Row {i+1}: Data size {data_list[i].shape[0]}x{data_list[i].shape[1]}, Subject ID: {participant_ids[i]}")
            
        print(f"\nExpected format: Each matrix should be regions x timepoints (dimensions will vary by file)")
            
    except Exception as e:
        print(f"Error saving file: {str(e)}")

if __name__ == "__main__":
    convert_csv_to_mat()