#!/usr/bin/env python3

import os
import sys

def create_condor_jobs_batch(list_dir, analysis_type="qa", output_sub_file=None):
    """
    Create condor submission files for batch processing where each job handles multiple input files.
    
    Args:
        list_dir: Directory containing the file lists (one list per job)
        analysis_type: Type of analysis ("qa" or "v2")
        output_sub_file: Output condor submission file name
    """
    
    if not os.path.exists(list_dir):
        print(f"Error: Directory {list_dir} not found")
        return
    
    # Transfer all necessary files for both analysis types
    transfer_files = "parton_qa.C, parton_v2.C, coordinateTools.h, binning.h"
    
    if output_sub_file is None:
        output_sub_file = f"condor_{analysis_type}.sub"
    
    # Get all list files in the directory
    list_files = []
    for filename in os.listdir(list_dir):
        file_path = os.path.join(list_dir, filename)
        if os.path.isfile(file_path) and not filename.startswith('.'):
            list_files.append(file_path)
    
    list_files.sort()  # Sort for consistent ordering
    
    print(f"Found {len(list_files)} list files in {list_dir}")
    
    # Create the condor submission file
    with open(output_sub_file, 'w') as f:
        f.write("Universe        = vanilla\n")
        f.write("GetEnv          = True\n")
        f.write(f"Executable      = {executable}\n")
        f.write(f"Log             = logs/{log_prefix}_$(Process).log\n")
        f.write(f"Output          = logs/{log_prefix}_$(Process).out\n")
        f.write(f"Error           = logs/{log_prefix}_$(Process).err\n")
        f.write("+MaxRuntime     = 20000\n\n")
        
        f.write("Should_Transfer_Files = YES\n")
        f.write(f"Transfer_Input_Files  = {transfer_files}\n\n")
        
        f.write("# Transfer_Output_Files = parton_qa_$(Process).root\n\n")
        
        # Add Arguments for each list file with analysis type
        for i, list_file in enumerate(list_files):
            f.write(f"Arguments_{i} = \"{list_file} {analysis_type}\"\n")
        
        f.write("\nArguments = $(Arguments_$(Process))\n")
        f.write(f"Queue {len(list_files)}\n")
    
    print(f"Created condor submission file: {output_sub_file}")
    print(f"To submit jobs: condor_submit {output_sub_file}")
    
    return output_sub_file

def create_condor_jobs_single(filelist_path, analysis_type="qa", output_sub_file="condor.sub"):
    """
    Create a condor submission file for single file processing (original functionality).
    
    Args:
        filelist_path: Path to the file containing the list of input files or a single .root file
        analysis_type: Type of analysis ("qa" or "v2")
        output_sub_file: Output condor submission file name
    """
    
    # Check if it's a single .root file or a file list
    if filelist_path.endswith('.root'):
        # Single .root file
        files = [filelist_path]
        print(f"Processing single .root file: {filelist_path}")
    else:
        # Read the file list
        with open(filelist_path, 'r') as f:
            files = [line.strip() for line in f if line.strip()]
        print(f"Found {len(files)} files to process from list: {filelist_path}")
    
    # Transfer all necessary files for both analysis types
    transfer_files = "parton_qa.C, parton_v2.C, coordinateTools.h, binning.h"
    
    # Create the condor submission file
    with open(output_sub_file, 'w') as f:
        f.write("Universe        = vanilla\n")
        f.write("Executable      = exec.sh\n\n")
        
        f.write("Should_Transfer_Files = YES\n")
        f.write(f"Transfer_Input_Files  = {transfer_files}\n\n")
        
        f.write("# Transfer_Output_Files = parton_qa_$(Process).root\n\n")
        
        f.write("Log              = logs/condor_$(Process).log\n")
        f.write("Output           = logs/condor_$(Process).out\n")
        f.write("Error            = logs/condor_$(Process).err\n")
        f.write("+MaxRuntime = 20000\n\n")
        
        # Add Arguments for each file with analysis type
        for i, filename in enumerate(files):
            f.write(f"Arguments_{i} = \"{filename} {analysis_type}\"\n")
        
        f.write("\nArguments = $(Arguments_$(Process))\n")
        f.write(f"Queue {len(files)}\n")
    
    print(f"Created condor submission file: {output_sub_file}")
    print(f"To submit jobs: condor_submit {output_sub_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python create_condor_jobs.py <list_directory_or_file> [analysis_type] [mode]")
        print("Examples:")
        print("  # Batch processing (multiple files per job)")
        print("  python create_condor_jobs.py Run2_tree_list/list_25/ qa batch")
        print("  python create_condor_jobs.py Run2_tree_list/list_25/ v2 batch")
        print("  # Single file processing (one file per job)")
        print("  python create_condor_jobs.py filelist.txt qa single")
        print("  python create_condor_jobs.py filelist.txt v2 single")
        print("  # Quick testing with single .root file")
        print("  python create_condor_jobs.py test.root qa single")
        print("  python create_condor_jobs.py test.root v2 single")
        sys.exit(1)
    
    input_path = sys.argv[1]
    analysis_type = sys.argv[2] if len(sys.argv) > 2 else "qa"
    mode = sys.argv[3] if len(sys.argv) > 3 else "batch"
    
    if mode == "batch":
        if not os.path.exists(input_path):
            print(f"Error: Directory {input_path} not found")
            sys.exit(1)
        create_condor_jobs_batch(input_path, analysis_type)
    elif mode == "single":
        if not os.path.exists(input_path):
            print(f"Error: File {input_path} not found")
            sys.exit(1)
        create_condor_jobs_single(input_path, analysis_type)
    else:
        print(f"Error: Unknown mode {mode}. Use 'batch' or 'single'")
        sys.exit(1)
