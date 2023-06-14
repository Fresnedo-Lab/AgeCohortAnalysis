# Make output folder. I keep these in scratch to save user drive space.

mkdir /fs/project/PAS0471/osu9453/AgeCohortProject/MergedLanes

# Make and travel to the batch script folder. The loop will generate and submit one script per sample within this folder.
	# I like to store loop scripts in:	$HOME/<path>/scripts
	# and batch scripts in: 			$HOME/<path>/batchScripts/<script name>
	
mkdir /fs/project/PAS0471/osu9453/AgeCohortProject/batchScripts/MergeLanes
cd /fs/project/PAS0471/osu9453/AgeCohortProject/batchScripts/MergeLanes

# Make sample lists - one for forward reads, and one for reverse reads.
	# Paired reads need to be listed in the same order. You can pipe sort just incase to address this.
	# These lists will be fed into the loop.
	
ls /fs/project/PAS0471/osu9453/AgeCohortProject/RawData/ | sed -E 's/_L.*fastq.*//' | uniq > list.txt

# Loop
	# n script files are created where n = the number of lines in list1.txt = the number of lines in list2.txt
	# list1.txt lines are defined as R1, list2.txt lines as R2
	
while IFS= read FILE
do
 
 OUTFILE=`echo $FILE | sed -E 's/R(.*)/R\1\.fq/'`	#use sed and regex to define basename for output.
 
 {
 echo '#PBS -l walltime=01:00:00'
 echo '#PBS -l nodes=1:ppn=1'
 echo '#PBS -A PAS0471'
 echo '#PBS -l mem=4GB'
 
 
 echo 'cd /fs/project/PAS0471/osu9453/AgeCohortProject/RawData/'
 
 echo "zcat $FILE* > ../MergedLanes/$OUTFILE"
 echo "gzip -f ../MergedLanes/$OUTFILE"
 } > mergeLanes_$OUTFILE.sh 	#echo'd text in brackets are written to a script. 1 script per listed sample.
 
 qsub mergeLanes_$OUTFILE.sh 	#submit the script
 
done < list.txt