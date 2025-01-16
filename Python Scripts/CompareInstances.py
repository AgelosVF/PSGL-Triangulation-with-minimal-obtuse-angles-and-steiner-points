import json
import os
import shutil
from pathlib import Path

def compare_and_copy_solutions(dir1, dir2, output_dir):
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Get all json files from first directory
    files1 = {f.name: f for f in Path(dir1).glob('*.json')}
    files2 = {f.name: f for f in Path(dir2).glob('*.json')}
    
    # Find common files
    common_files = set(files1.keys()) & set(files2.keys())
    
    for filename in common_files:
        # Read both files
        with open(files1[filename], 'r') as f1:
            data1 = json.load(f1)
        with open(files2[filename], 'r') as f2:
            data2 = json.load(f2)
        
        # Get obtuse counts and convert to integers
        count1 = int(data1['obtuse_count'])
        count2 = int(data2['obtuse_count'])
        
        # Determine which file to copy
        if count1 < count2:
            source_file = files1[filename]
        elif count2 < count1:
            source_file = files2[filename]
        else:
            # If obtuse counts are equal, compare number of Steiner points
            steiner_points1 = len(data1['steiner_points_x'])
            steiner_points2 = len(data2['steiner_points_x'])
            source_file = files1[filename] if steiner_points1 <= steiner_points2 else files2[filename]
        
        # Copy the selected file to output directory
        shutil.copy2(source_file, output_dir)
        
        # Print information about the comparison
        print(f"Processing {filename}:")
        print(f"  Dir1 obtuse count: {count1}, Steiner points: {len(data1['steiner_points_x'])}")
        print(f"  Dir2 obtuse count: {count2}, Steiner points: {len(data2['steiner_points_x'])}")
        print(f"  Copied from: {source_file}")
        print()

# Example usage
if __name__ == "__main__":
    directory1 = "third"
    directory2 = "3_CycleConvex_All_Results_SA_All2"
    output_directory = "Third"
    
    compare_and_copy_solutions(directory1, directory2, output_directory)
