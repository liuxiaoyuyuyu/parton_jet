#!/usr/bin/env python3

import sys
import os

def create_condor_submission(filename_list, output_sub_file="condor_parton_v2.sub"):
    """
    Create a condor submission file based on a list of filenames.
    
    Args:
        filename_list: List of input filenames
        output_sub_file: Output condor submission file name
    """
    
    with open(output_sub_file, 'w') as f:
        f.write("Universe        = vanilla\n")
        f.write("Executable      = exec.sh\n\n")
        
        f.write("Should_Transfer_Files = YES\n")
        f.write("Transfer_Input_Files  = parton_v2.C, binning.h, coordinateTools.h\n\n")
        
        #f.write("# Transfer output files back\n")
        #f.write("Transfer_Output_Files = parton_v2_output_$(Process).root\n\n")
        
        f.write("Log              = logs/condor_parton_v2_$(Process).log\n")
        f.write("Output           = logs/condor_parton_v2_$(Process).out\n")
        f.write("Error            = logs/condor_parton_v2_$(Process).err\n")
        f.write("+MaxRuntime = 3600\n\n")
        
        # Write individual job entries
        for i, filename in enumerate(filename_list):
            f.write(f"Arguments = \"{filename}\"\n")
            f.write(f"Queue 1\n\n")
    
    print(f"Created condor submission file: {output_sub_file}")
    print(f"Total jobs: {len(filename_list)}")

def read_filename_list(list_file):
    """
    Read filenames from a text file.
    
    Args:
        list_file: Path to file containing list of filenames
    
    Returns:
        List of filenames
    """
    filenames = []
    with open(list_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                filenames.append(line)
    return filenames

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_condor_jobs.py <filename_list.txt>")
        print("The filename_list.txt should contain one filename per line")
        sys.exit(1)
    
    list_file = sys.argv[1]
    if not os.path.exists(list_file):
        print(f"Error: File {list_file} does not exist")
        sys.exit(1)
    
    filenames = read_filename_list(list_file)
    if not filenames:
        print("Error: No filenames found in the list file")
        sys.exit(1)
    
    create_condor_submission(filenames)
    print(f"Processed {len(filenames)} files") 