import os
import subprocess
import json

def process_json_files():
    # Get all immediate subdirectories
    subdirs = [d for d in os.listdir('.') if os.path.isdir(d)]
    
    if not subdirs:
        print("No subdirectories found in current directory")
        return
    
    for subdir in subdirs:
        # Create results directory for this subdirectory
        results_dir = f'{subdir}_results'
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
            print(f"Created directory: {results_dir}")
        
        # Get all .json files in the current subdirectory
        json_files = [f for f in os.listdir(subdir) if f.endswith('.json')]
        
        for json_file in json_files:
            # Create output filename in results directory
            input_file = os.path.join(subdir, json_file)
            output_file = os.path.join(results_dir, f"{json_file}")
            
            # Construct the command
            command = ['./opt_triangulation', '-i', input_file, '-o', output_file, '-preselected_params']
            
            print(f"Processing {input_file}...")
            try:
                # Run the command and capture output
                result = subprocess.run(command, 
                                     capture_output=True, 
                                     text=True, 
                                     check=True)
                print(f"Successfully processed {json_file}")
                print(f"Output saved to: {output_file}")
                
            except subprocess.CalledProcessError as e:
                print(f"Error processing {json_file}:")
                print(f"Error output: {e.stderr}")
                
            except Exception as e:
                print(f"Unexpected error processing {json_file}: {str(e)}")

if __name__ == "__main__":
    process_json_files()
   #count_obtuse_zeros()
