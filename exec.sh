#!/bin/bash

# Check if argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <file_or_filelist>"
    echo "Example: $0 /path/to/single_file.root"
    echo "Example: $0 /path/to/filelist.txt"
    exit 1
fi

INPUT="$1"

# Check if input is a file list (ends with .txt) or a single file
if [[ "$INPUT" == *.txt ]]; then
    echo "Processing file list: $INPUT"
    # Check if file list exists
    if [ ! -f "$INPUT" ]; then
        echo "Error: File list $INPUT not found"
        exit 1
    fi
    
    # Determine analysis type based on the condor job name or other logic
    # For now, default to QA analysis for file lists
    echo "Running QA analysis on file list"
    root -l -b -q "parton_qa.C(\"$INPUT\", \"/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/anaOutput/qa/\")"
else
    echo "Processing single file: $INPUT"
    # Check if file exists
    if [ ! -f "$INPUT" ]; then
        echo "Error: File $INPUT not found"
        exit 1
    fi
    
    # For single files, run v2 analysis (original behavior)
    echo "Running v2 analysis on single file"
    root -l -b -q "parton_v2.C(\"$INPUT\", \"/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/anaOutput/round5/\")"
fi