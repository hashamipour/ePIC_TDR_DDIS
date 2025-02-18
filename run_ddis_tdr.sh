#!/bin/bash

# Check if two arguments are provided
if [ $# -ne 1]; then
  echo "Usage: $0 <list_file>"
  exit 1  # Exit with an error code
fi

# Assign the arguments to variables
list_file="$1"

# Now use the variables in your ROOT command
root -b -q "./DDIS_TDR.cpp(\"${list_file}\")"

# Or, if you have compiled your code:
# ./myprogram "$rec_file" "$outputfile"

# Optional: Print the filenames for confirmation (good for debugging)
#echo "Reconstructed file: $rec_file"
#echo "Output file: $outputfile"
