#!/bin/bash

#Navigate to the folder containing contamination report files from NCBI
#Store the name of the folder in a variable
contamination_folder="/Users/shreyas/Documents/bmc_seq_submission/contamination_reports"
#Navigate to the folder
cd $contamination_folder

#loop through the files
for report in *
	do
		echo $report
		#Each contamination report file has three parts: base info, reads to be excluded because of contamination, and the reads to be trimmed
		#The goal is to make two files for each report: one for excluding reads, and another for trimming
		#Store name of the file without the extension into a variable
		name="$(cut -d'.' -f1 <<<"$report")"
		echo $name
		#Remove everything above the line "Exclude:"
		sed -e '1,/Exclude:/d' $report > ${name}_exclude_trim.txt
		#Break down the file created above into two: one for exclude and one for trim
		exclude=${name}_exclude.txt
		trim=${name}_trim.txt
		sed -e '/Trim:/,$d' ${name}_exclude_trim.txt > $exclude  
		sed -e '1,/Trim:/d' ${name}_exclude_trim.txt > $trim
	done
