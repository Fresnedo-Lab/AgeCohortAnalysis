# Make output folder. I keep these in scratch to save user drive space.

mkdir /fs/project/PAS0471/osu9453/AgeCohortProject/FASTQC_PreTrim_Out

# Make and travel to the batch script folder. The loop will generate and submit one script per sample within this folder.
	# I like to store loop scripts in:	$HOME/<path>/scripts
	# and batch scripts in: 			$HOME/<path>/batchScripts/<script name>
	
mkdir /fs/project/PAS0471/osu9453/AgeCohortProject/batchScripts/FASTQC_Pre
cd /fs/project/PAS0471/osu9453/AgeCohortProject/batchScripts/FASTQC_Pre

# Make sample lists - one for forward reads, and one for reverse reads.
	# Paired reads need to be listed in the same order. You can pipe sort just incase to address this.
	# These lists will be fed into the loop.
	
ls /fs/project/PAS0471/osu9453/AgeCohortProject/RawData/*R1*gz | sort > list1.txt
ls /fs/project/PAS0471/osu9453/AgeCohortProject/RawData/*R2*gz | sort > list2.txt


# Loop
	# n script files are created where n = the number of lines in list1.txt = the number of lines in list2.txt
	# list1.txt lines are defined as R1, list2.txt lines as R2
	
while IFS= read R1 && IFS= read R2 <&3
do

 ID1=`echo $R1 | xargs basename | sed -E 's/\.fastq\.gz//'`	#use sed and regex to define basename for output.
 ID2=`echo $R2 | xargs basename | sed -E 's/\.fastq\.gz//'`
 
 {
 echo '#PBS -l walltime=00:30:00'
 echo '#PBS -l nodes=1:ppn=2'
 echo '#PBS -A PAS0471'
 echo '#PBS -l mem=8GB'
 
 echo 'module load fastqc/0.11.7'
 
 echo 'cd /fs/project/PAS0471/osu9453/AgeCohortProject/FASTQC_PreTrim_Out'
 
 echo "fastqc -t 2 -o /fs/project/PAS0471/osu9453/AgeCohortProject/FASTQC_PreTrim_Out \
 $R1 $R2"
 } > fastqc_pre_$ID1.sh 	#echo'd text in brackets are written to a script. 1 script per listed sample.
 
 qsub fastqc_pre_$ID1.sh 	#submit the script
 
done < list1.txt 3< list2.txt