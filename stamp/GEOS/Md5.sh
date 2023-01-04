#!/bin/bash

# create an empty file called "output.csv"
touch output.csv

# loop through all files in the current directory
for file in *; do
  # get the md5 hash of the file
  hash=$(md5sum "$file" | cut -d ' ' -f 1)
  # append the filename and hash to "output.csv"
  echo "$file,$hash" >> output.csv
done

# print a message to let the user know that the script is finished
echo "Done!"
