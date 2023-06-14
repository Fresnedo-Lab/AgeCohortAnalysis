# Make output folder. I keep these in scratch to save user drive space.

mkdir /fs/project/PAS0471/osu9453/AgeCohortProject/Bismark_Alignment_Out

# Make and travel to the batch script folder. The loop will generate and submit one script per sample within this folder.
	# I like to store loop scripts in:	$HOME/<path>/scripts
	# and batch scripts in: 			$HOME/<path>/batchScripts/<script name>
	
mkdir /fs/project/PAS0471/osu9453/AgeCohortProject/batchScripts/Bismark_Alignment
cd /fs/project/PAS0471/osu9453/AgeCohortProject/batchScripts/Bismark_Alignment

# Make sample lists - one for forward reads, and one for reverse reads.
	# Paired reads need to be listed in the same order. You can pipe sort just incase to address this.
	# These lists will be fed into the loop.
	
ls /fs/project/PAS0471/osu9453/AgeCohortProject/trimgaloreOutPE/*R1_trimmed*gz | sort > list1.txt
ls /fs/project/PAS0471/osu9453/AgeCohortProject/trimgaloreOutPE/*R2_trimmed*gz | sort > list2.txt


# Loop
	# n script files are created where n = the number of lines in list1.txt = the number of lines in list2.txt
	# list1.txt lines are defined as R1, list2.txt lines as R2
	
while IFS= read R1 && IFS= read R2 <&3
do

 ID1=`echo $R1 | xargs basename | sed -E 's/\_trimmed\.fq\.gz//'`	#use sed and regex to define basename for output.
 ID2=`echo $R2 | xargs basename | sed -E 's/\_trimmed\.fq\.gz//'`
 
 {
 echo '#PBS -l walltime=10:00:00'
 echo '#PBS -l nodes=1:ppn=28'
 echo '#PBS -A PAS0471'
 echo '#PBS -l mem=8GB'
 
 echo 'module use $HOME/local/share/lmodfiles'
 echo 'module load module load bowtie2/2.4.1'
 echo 'module load module load bismark/0.22.3'
 echo 'module load samtools/1.9'
 
 echo 'cd /fs/project/PAS0471/osu9453/AgeCohortProject/Bismark_Alignment_Out'
 
 echo "bismark --genome /fs/project/PAS0471/osu9453/AgeCohortProject/genome \
  --path_to_bowtie2 /usr/local/bowtie2/2.4.1 --samtools_path /usr/local/samtools/1.9/bin --parallel 28 \
  -1 $R1 -2 $R2 
"
 } > bismark_alignment_$ID1.sh 	#echo'd text in brackets are written to a script. 1 script per listed sample.
 
 qsub bismark_alignment_$ID1.sh 	#submit the script
 
done < list1.txt 3< list2.txt