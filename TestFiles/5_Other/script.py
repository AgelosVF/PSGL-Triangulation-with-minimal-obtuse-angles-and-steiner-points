import os
import subprocess

def process_json_files():
    # Create results directory if it doesn't exist
    results_dir = 'results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        print(f"Created directory: {results_dir}")
    
    # Get all .json files in the current directory
    json_files = [f for f in os.listdir('.') if f.endswith('.json')]
    
    for json_file in json_files:
        # Create output filename in results directory
        output_file = os.path.join(results_dir, f"output_{json_file}")
        
        # Construct the command
        command = ['./Triangulation', '-i', json_file, '-o', output_file]
        
        print(f"Processing {json_file}...")
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
