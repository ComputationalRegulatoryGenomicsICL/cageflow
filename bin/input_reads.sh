#!/bin/bash

# I will add here the file matching rather than path input to allow for the input to be in multiple folders
# I would also probably add input and output as $1 and $2
# files=$(ls /path/to/your/root/folder/*fastq.gz)
#folder_path="/path/to/your/folder"
folder_path="./customcageq/assets/mock_fq"
#output_csv="./output.csv"
output_csv="${folder_path}/output.csv"

# Write the header to the CSV
echo "sample,fastq_1,fastq_2,single_end" > "${output_csv}"

# Loop through each R1 file
for r1 in "${folder_path}"/*_R1*; do
    # Extract sample name (assumes format: sampleName_R1...)
    sample_name=$(basename "$r1" | cut -d'_' -f1)

    # Check for corresponding R2 file
    full_sample_name=$(basename "$r1" | cut -d"R" -f1)
    r2="${folder_path}/${full_sample_name}R2*"

    # Do not forget to remove these commented lines
    #r2_ending=$(basename "$r1" | cut -d"R" -f2 | cut -c2-)
    #r2="${folder_path}/${sample_name}*R2*"

    if [ -f $r2 ]; then

        # Do not forget to remove these commented lines
        #echo "${sample_name},$(basename "$r1"),$(echo $(basename "$r2")$r2_ending | tr -d '*'),False" >> "${output_csv}"

        echo "${sample_name},$(basename "$r1"),$(basename $(ls $r2)),False" >> "${output_csv}"
    else
        echo "${sample_name},$(basename "$r1"),,True" >> "${output_csv}"
    fi
done

echo "CSV file created at: ${output_csv}"