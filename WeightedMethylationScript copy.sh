#This is the script I used to generate the total number of methylated reads for every cytosine in the genome for each of the files
for i in *; do awk '{s+=$4} END {print s}' $i; done
#This is the script I used to generate the total number of reads (methylated + unmethylated) for every cytosine in the genome for each of the files
for i in *; do awk '{s+=$3} END {print s}' $i; done

##I put these values into Excel for each sample and divided methylated reads by total number of reads to get percent methylation
##What I'd like to do is subset the text files by chromosome, which is in column 1 of the text file, and listed as chr1, chr2, chr3,...,chr8
##And then calculate weighted percent methylation for all individuals for each methylation context and chromosome.