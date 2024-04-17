#!/bin/bash
cd ./temp &&
wait
# Get a list of reference files
reference_files=(*_pr_on_rl.txt)


# Loop through each reference file
for reference_file in "${reference_files[@]}"; do
  output_file="${reference_file%.txt}_Output.txt"

  # Extract the prefix of the reference file
  suffix="_pr_on_rl.txt"
  prefix="${reference_file%"$suffix"}"

  # Get the files associated with the current reference file
  files=("${prefix}_pr_on_perfish_I_rl.txt" "${prefix}_pr_sg1indel_II_rl.txt" "${prefix}_pr_R1wIII_rl.txt" "${prefix}_pr_R2wIII_rl.txt" "${prefix}_pr_sg4indel_IV_rl.txt" "${prefix}_pr_ldel1to3_V_rl.txt" "${prefix}_pr_ldel3to4_VI_rl.txt" "${prefix}_pr_ldel1to4R1_VIIr1_rl.txt" "${prefix}_pr_ldel1to4R2_VIIr2_rl.txt")

  # Read the reference file line by line and store the items in an array
  reference_items=()
  while IFS= read -r line; do
    reference_items+=("$line")
  done < "$reference_file"

  # Create the output file
  > "$output_file"

  # Loop through each item in the reference file
  for item in "${reference_items[@]}"; do
    line="$item "

    # Check if the item is present in each file
    for file in "${files[@]}"; do
      if grep -Fxq "$item" "$file"; then
        line+="Y"
      else
        line+="N"
      fi
    done

    # Append the line to the output file
    echo "$line" >> "$output_file"
  done
done
cd ../ 