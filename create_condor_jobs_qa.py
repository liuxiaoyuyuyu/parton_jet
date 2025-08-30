#!/usr/bin/env python3

import os
import sys

def create_condor_jobs(filelist_path, output_sub_file="condor_qa.sub"):
    """
    Create a condor submission file for QA analysis jobs.
    
    Args:
        filelist_path: Path to the file containing the list of input files
        output_sub_file: Output condor submission file name
    """
    
    # Read the file list
    with open(filelist_path, 'r') as f:
        files = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(files)} files to process")
    
    # Create the condor submission file
    with open(output_sub_file, 'w') as f:
        f.write("Universe        = vanilla\n")
        f.write("Executable      = exec.sh\n\n")
        
        f.write("Should_Transfer_Files = YES\n")
        f.write("Transfer_Input_Files  = parton_qa.C\n\n")
        
        f.write("# Transfer_Output_Files = parton_qa_$(Process).root\n\n")
        
        f.write("Log              = logs/condor_qa_$(Process).log\n")
        f.write("Output           = logs/condor_qa_$(Process).out\n")
        f.write("Error            = logs/condor_qa_$(Process).err\n")
        f.write("+MaxRuntime =1800\n\n")
        
        # Add Arguments for each file
        for i, filename in enumerate(files):
            f.write(f"Arguments_{i} = \"{filename}\"\n")
        
        f.write("\nArguments = $(Arguments_$(Process))\n")
        f.write(f"Queue {len(files)}\n")
    
    print(f"Created condor submission file: {output_sub_file}")
    print(f"To submit jobs: condor_submit {output_sub_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_condor_jobs_qa.py <filelist>")
        print("Example: python create_condor_jobs_qa.py filelist.txt")
        sys.exit(1)
    
    filelist_path = sys.argv[1]
    if not os.path.exists(filelist_path):
        print(f"Error: File list {filelist_path} not found")
        sys.exit(1)
    
    create_condor_jobs(filelist_path)
