#!/bin/bash

# Check if output filename provided
if [ -z $1 ]; then
	out=bam_statistics.txt
else
	out=$1
fi

dup_files="_duplicate_metrics.txt"
hs_files=".HS_Metrics.txt"

rm -f temp1_$out temp2_$out

# Output samples and duplicate statistics
for f in *${dup_files}; do
    sample=${f%${dup_files}}
    dup=`awk '{if (NR==8) {gsub ( /\./,"," ) ; print $9}}' $f`
    # echo -e "sample\tdup" >> temp1_$out
    echo -e "${sample}\t${dup}" >> temp1_$out
done

# Get HS Metrics
for f in *${hs_files}; do
    # awk '{if (NR==7) {print $0}}' $f >> temp2_$out
    awk '{if (NR==8) {gsub ( /\./,"," ) ; print $0}}' $f >> temp2_$out
done

paste -d "\t" temp1_$out temp2_$out > $out
rm -f temp1_$out temp2_$out

# for d in */; do echo -e "${d::-1}" >> $file; done
#
# # Get HS Metrics
# echo -e "\n\n# Get HS Metrics (sample order as above)" >> $file
# for d in */; do cat ${d}${d::-1}_HS_Metrics.txt | grep ^MedExome_hg19_capture_targets | sed '/[0-9]\./s/\./,/g' >> $file; done
#
# # Get duplicate statistics
# echo -e "\n\n# Duplicate statistics (sample order as above)" >> $file
# for d in */; do cat ${d}duplicate_metrics.txt | grep ^NimbleGen | awk '{gsub ( /\./,"," ) ; print $9}' >> $file; done
#
