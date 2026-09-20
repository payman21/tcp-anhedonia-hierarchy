import os
import pandas as pd
from pathlib import Path
import csv

def find_results_directories():
    """
    Find all Results directories in NDAR folders, regardless of nesting level
    """
    current_dir = Path('.')
    results_paths = []
    
    # Walk through all directories to find NDAR folders and their Results subdirectories
    for root, dirs, files in os.walk(current_dir):
        root_path = Path(root)
        
        # Check if this is a Results directory
        if root_path.name == 'Results':
            # Find the NDAR folder by walking up the path
            ndar_folder = None
            for parent in root_path.parents:
                if parent.name.startswith('NDAR_'):
                    ndar_folder = parent.name
                    break
            
            if ndar_folder:
                results_paths.append({
                    'ndar_folder': ndar_folder,
                    'results_path': root_path
                })
    
    return results_paths

def scan_folder_structure():
    """
    Scans for NDAR folders and creates a CSV report of their Results subdirectory structure.
    """
    print("Searching for Results directories in NDAR folders...")
    
    # Find all Results directories
    results_data = find_results_directories()
    
    if not results_data:
        print("No Results directories found in NDAR folders.")
        return
    
    print(f"Found Results directories in {len(results_data)} NDAR folders")
    
    # Dictionary to store all unique subfolders and files found across all Results directories
    all_subfolders = set()
    all_files = set()
    
    # Dictionary to store the data for each NDAR folder
    folder_data = {}
    
    # Process each Results directory
    for result_info in results_data:
        ndar_folder = result_info['ndar_folder']
        results_path = result_info['results_path']
        
        if ndar_folder not in folder_data:
            folder_data[ndar_folder] = {'subfolders': set(), 'files': set()}
        
        print(f"Processing: {ndar_folder}")
        
        # Get all subdirectories in Results
        for item in results_path.iterdir():
            if item.is_dir():
                subfolder_name = item.name
                all_subfolders.add(subfolder_name)
                folder_data[ndar_folder]['subfolders'].add(subfolder_name)
                
                # Get all files in this subfolder
                for file_item in item.iterdir():
                    if file_item.is_file():
                        file_key = f"{subfolder_name}/{file_item.name}"
                        all_files.add(file_key)
                        folder_data[ndar_folder]['files'].add(file_key)
    
    # Sort for consistent output
    all_subfolders = sorted(all_subfolders)
    all_files = sorted(all_files)
    
    print(f"\nFound {len(all_subfolders)} unique subfolders across all NDAR folders:")
    for subfolder in all_subfolders:
        print(f"  - {subfolder}")
    
    print(f"\nFound {len(all_files)} unique files across all subfolders")
    
    # Create the CSV data
    csv_data = []
    
    # Create header row
    header = ['NDAR_Folder']
    
    # Add subfolder columns
    for subfolder in all_subfolders:
        header.append(f"Subfolder_{subfolder}")
    
    # Add file columns  
    for file_path in all_files:
        # Clean up the file path for column name
        clean_name = file_path.replace('/', '_').replace(' ', '_')
        header.append(f"File_{clean_name}")
    
    csv_data.append(header)
    
    # Create data rows
    unique_ndar_folders = sorted(set(item['ndar_folder'] for item in results_data))
    
    for ndar_folder in unique_ndar_folders:
        row = [ndar_folder]
        
        # Add subfolder existence data
        for subfolder in all_subfolders:
            if ndar_folder in folder_data:
                exists = "yes" if subfolder in folder_data[ndar_folder]['subfolders'] else "no"
            else:
                exists = "no"
            row.append(exists)
        
        # Add file existence data
        for file_path in all_files:
            if ndar_folder in folder_data:
                exists = "yes" if file_path in folder_data[ndar_folder]['files'] else "no"
            else:
                exists = "no"
            row.append(exists)
        
        csv_data.append(row)
    
    # Write to CSV file
    output_file = 'ndar_results_structure.csv'
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(csv_data)
    
    print(f"\nCSV file created: {output_file}")
    print(f"Total rows: {len(csv_data)} (including header)")
    print(f"Total columns: {len(header)}")
    
    # Print summary statistics
    print(f"\nProcessed {len(unique_ndar_folders)} unique NDAR folders")
    
    print("\nSubfolder presence summary:")
    for subfolder in all_subfolders:
        count = sum(1 for ndar in unique_ndar_folders 
                   if ndar in folder_data and subfolder in folder_data[ndar]['subfolders'])
        percentage = (count / len(unique_ndar_folders)) * 100
        print(f"  {subfolder}: {count}/{len(unique_ndar_folders)} folders ({percentage:.1f}%)")
    
    # Show sample of files found
    print(f"\nSample of files found (showing first 10 of {len(all_files)}):")
    for file_path in list(all_files)[:10]:
        print(f"  {file_path}")
    
    if len(all_files) > 10:
        print(f"  ... and {len(all_files) - 10} more files")

def main():
    print("Scanning NDAR folder structure for Results directories...")
    print("=" * 60)
    
    try:
        scan_folder_structure()
        print("\nScan completed successfully!")
        print("\nThe CSV file shows:")
        print("- One row per NDAR folder")
        print("- One column per task subfolder (e.g., task-stroopPA_run-01_bold)")
        print("- One column per file found in any subfolder")
        print("- 'yes'/'no' values indicating presence of each subfolder/file")
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()