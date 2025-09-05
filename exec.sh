#!/bin/bash

# Check if arguments are provided
if [ $# -lt 1 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <file_or_filelist> [analysis_type] [process_id]"
    echo "Example: $0 /path/to/single_file.root qa 0"
    echo "Example: $0 /path/to/filelist.txt v2 1"
    echo "Example: $0 /path/to/single_file.root (defaults to v2 for single files, process_id=0)"
    exit 1
fi

INPUT="$1"
ANALYSIS_TYPE="$2"
PROCESS_ID="$3"

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

# Set default process ID if not provided
if [ -z "$PROCESS_ID" ]; then
    PROCESS_ID="0"
fi

echo "Processing: $INPUT"
echo "Analysis type: $ANALYSIS_TYPE"

# Run the appropriate analysis
if [ "$ANALYSIS_TYPE" = "qa" ]; then
    echo "Running QA analysis"
    root -l -b -q "parton_qa.C(\"$INPUT\", \"/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/anaOutput/qa/parton_qa_output_${PROCESS_ID}.root\")"
elif [ "$ANALYSIS_TYPE" = "v2" ]; then
    echo "Running v2 analysis"
    root -l -b -q "parton_v2.C(\"$INPUT\", \"/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/anaOutput/round5/parton_v2_output_${PROCESS_ID}.root\")"
else
    echo "Error: Unknown analysis type '$ANALYSIS_TYPE'. Use 'qa' or 'v2'"
    exit 1
fi