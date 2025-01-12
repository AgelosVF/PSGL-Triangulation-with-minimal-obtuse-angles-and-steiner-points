import os
import subprocess
import shutil

# Create directories if they don't exist
for i in range(1, 6):
    dir_path = f"../{i}"
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

# Get all JSON files in current directory
json_files = [f for f in os.listdir('.') if f.endswith('.json')]

for json_file in json_files:
    # Create temporary output file name
    output_file = "temp_output.json"
    
    # Run the optimize_triangulation command
    try:
        result = subprocess.run(
            ['./optimize_triangulation', '-i', json_file, '-o', output_file, '-preselected_params'],
            capture_output=True,
            text=True
        )
        
        # Extract the type from the output
        for line in result.stdout.split('\n'):
            if "Type:" in line:
                type_num = int(line.split(':')[1].strip())
                
                # Copy the original json file to the appropriate directory
                dest_dir = f"../{type_num}"
                dest_file = os.path.join(dest_dir, json_file)
                shutil.copy2(json_file, dest_file)
                print(f"Copied {json_file} to {dest_dir} (Type {type_num})")
                break
                
    except subprocess.CalledProcessError as e:
        print(f"Error processing {json_file}: {e}")
    except Exception as e:
        print(f"Unexpected error with {json_file}: {e}")
    
    # Clean up temporary output file
    if os.path.exists(output_file):
        os.remove(output_file)
