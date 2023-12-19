#!/bin/bash

# I will add here the file matching rather than path input to allow for the input to be in multiple folders
# I would also probably add input and output as $1 and $2
# files=$(ls /path/to/your/root/folder/*fastq.gz)
# folder_path="/path/to/your/folder"
# folder_path="/Users/pavel/Desktop/PROJECTS/hooman-2/customcageq/assets/mock_fq"

# Check that there are exactly two positional arguments
if [ $# -ne 2 ]; then
    echo "Error: Provide exactly 2 positional arguments: a directory with FASTQ files and the file name of the CSV table to create."
    exit 1
fi

# Check that the first argument is an existing directory
if ! [ -d $1 ]; then
    echo "Error: The first argument must be an existing directory."
    exit 1
fi

# Assign positional argument values to variables for further use
folder_path="$1"
output_csv="$2"

# Write the header to the CSV
echo "sample,fastq_1,fastq_2,single_end" > "${output_csv}"

# Loop through each R1 file
for r1 in "${folder_path}"/*_R1*; do
    
    # Extract sample name (assumes format: sampleName_R1...)
    sample_name=$(basename "$r1" | cut -d'_' -f1)

    # Check for corresponding R2 file
    full_sample_name=$(basename "$r1" | awk -F'_R[1-2]' '{print $1}')
    r2="${folder_path}/${full_sample_name}_R2*"

    if [ -f $r2 ]; then
        echo "${sample_name},$r1,$(ls $r2),False" >> "${output_csv}"
    else
        echo "${sample_name},$r1,,True" >> "${output_csv}"
    fi
done

echo "CSV file created at: ${output_csv}"
