#!/bin/bash

# Check if arguments are provided
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <file_or_filelist> [analysis_type]"
    echo "Example: $0 /path/to/single_file.root qa"
    echo "Example: $0 /path/to/filelist.txt v2"
    echo "Example: $0 /path/to/single_file.root (defaults to v2 for single files)"
    exit 1
fi

INPUT="$1"
ANALYSIS_TYPE="$2"

# Check if input file exists
if [ ! -f "$INPUT" ]; then
    echo "Error: Input file $INPUT not found"
    exit 1
fi

# Determine analysis type
if [ -z "$ANALYSIS_TYPE" ]; then
    # Default analysis type based on file type
    if [[ "$INPUT" == *.root ]]; then
        ANALYSIS_TYPE="v2"  # Default to v2 for single .root files
    else
        ANALYSIS_TYPE="qa"  # Default to qa for file lists
    fi
fi

echo "Processing: $INPUT"
echo "Analysis type: $ANALYSIS_TYPE"

# Run the appropriate analysis
if [ "$ANALYSIS_TYPE" = "qa" ]; then
    echo "Running QA analysis"
    root -l -b -q "parton_qa.C(\"$INPUT\", \"/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/anaOutput/qa/\")"
elif [ "$ANALYSIS_TYPE" = "v2" ]; then
    echo "Running v2 analysis"
    root -l -b -q "parton_v2.C(\"$INPUT\", \"/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/anaOutput/round5/\")"
else
    echo "Error: Unknown analysis type '$ANALYSIS_TYPE'. Use 'qa' or 'v2'"
    exit 1
fi