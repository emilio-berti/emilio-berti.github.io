#!/bin/bash

# Check if a file argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 filename"
    exit 1
fi

file="$1"

# Check if the file exists
if [ ! -f "$file" ]; then
    echo "File '$file' not found!"
    exit 1
fi

# Search for the string "install"
if grep -q "install" "$file"; then
    # Print red warning
    echo -e "Install check: \033[31mWARNING\033[0m"
else
    echo -e "Install check: \033[32mOK\033[0m"
fi
