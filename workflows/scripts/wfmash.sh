#!/bin/bash

ref="$1"
query_files="$2"
output_dir="$3"

count=0

while IFS= read -r filename; do
	basename="${filename%.*}"
	query=$filename
	order_file="order.txt"
	wfmash_out=".txt"
	final_fasta=$output_dir/out_${count}.fa
	
	wfmash $ref $query -m -p 95 -N -t 16 > $wfmash_out
	
	sort -k9,9n $wfmash_out | awk '{print $1}' > "$order_file"
	
	echo "Segment names sorted and saved to $order_file"
	
	concat_fasta() {
	    local fasta_file="$1"
	    local order_file="$2"
	    local output_file="$3"
	
	    declare -A sequences
	
	    local current_seq=""
	    while IFS= read -r line || [[ -n $line ]]; do
		    # Remove '\n' from the line
		    line=$(echo "$line" | tr -d '\n')
		    
		    if [[ $line == ">"* ]]; then
		        # Start of a new sequence header
		        current_seq=${line#>}
		        echo "$current_seq"
		        sequences["$current_seq"]=""
		    else
		        # Append the line to the current sequence
		        sequences["$current_seq"]+="$line"
		    fi
	    done < "$fasta_file"
	
	    > "$output_file"
	    echo ">$basename" >> "$output_file"
	    output_data=""
	    while IFS= read -r seq_name || [[ -n $seq_name ]]; do
	        if [[ -n ${sequences["$seq_name"]} ]]; then
	            output_data+="${sequences["$seq_name"]}"$'\n'
	        else
	            echo "Warning: Sequence '$seq_name' not found in FASTA file" >&2
	        fi
	    done < "$order_file"
	    echo -n "$output_data" >> "$output_file"
	    
	    echo "Concatenated sequences written to $output_file"
	}
	
	concat_fasta $query $order_file $final_fasta
	count=$((count + 1))

done < "$query_files"
