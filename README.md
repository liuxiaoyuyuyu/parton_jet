# Parton Jet Analysis

This repository contains ROOT macros for analyzing parton jet data with Condor job submission support.

## Files

- `parton_qa.C` - Quality assurance analysis
- `parton_v2.C` - v2/eccentricity analysis  
- `create_condor_jobs.py` - Condor job submission script
- `exec.sh` - Execution script for Condor jobs
- `coordinateTools.h` - Coordinate transformation utilities
- `binning.h` - Binning definitions

## Quick Start

### 1. Prepare File Lists

Create file lists containing multiple `.root` files:

```bash
find /eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/sample/ -name "*batch4*">>Treelist_zpc.list
find /eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/sample/ -name "*batch5*">>Treelist_zpc.list
...

mkdir -p Model_tree_list/list_25/
cd Model_tree_list/list_25/
split -l25 -d -a 3 ../../Treelist_zpc list_job
```

### 2. Generate Condor Jobs

```bash
# Batch processing (multiple files per job)
python create_condor_jobs.py Model_tree_list/list_25/ qa batch
python create_condor_jobs.py Model_tree_list/list_25/ v2 batch

# Single file processing (one file per job)
python create_condor_jobs.py filelist.txt qa single
python create_condor_jobs.py filelist.txt v2 single

# Quick testing with single .root file
python create_condor_jobs.py test.root qa single
python create_condor_jobs.py test.root v2 single
```

### 3. Submit Jobs

```bash
mkdir -p logs
condor_submit condor_qa.sub    # For QA analysis
condor_submit condor_v2.sub    # For v2 analysis
condor_submit condor.sub       # For single mode
```

### 4. Monitor Jobs

```bash
condor_q                       # Check job status
tail -f logs/condor_0.out      # Check job output
```

## Analysis Types

- **qa**: Quality assurance analysis with jet multiplicity, collision correlations, etc.
- **v2**: v2/eccentricity analysis with proper time evolution

## Output Files

- **QA**: `/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/anaOutput/qa/`
- **v2**: `/eos/cms/store/group/phys_heavyions/xiaoyul/wenbin/anaOutput/round5/`

## File List Format

File lists can have any name and contain one `.root` file path per line:

```
/path/to/file1.root
/path/to/file2.root
/path/to/file3.root
```

Comments (lines starting with `#`) and empty lines are ignored.
