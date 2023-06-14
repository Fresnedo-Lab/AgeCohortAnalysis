#text

mkdir /fs/project/PAS0471/osu9453/AgeCohortProject/batchScripts/trimgalorePE
cd /fs/project/PAS0471/osu9453/AgeCohortProject/batchScripts/trimgalorePE

mkdir -p /fs/project/PAS0471/osu9453/AgeCohortProject/trimgaloreOutPE

ls /fs/project/PAS0471/osu9453/AgeCohortProject/MergedLanes/*R1* > listPE_R1.txt
ls /fs/project/PAS0471/osu9453/AgeCohortProject/MergedLanes/*R2* > listPE_R2.txt

while IFS= read R1 && IFS= read R2 <&3
do
SAMPID=`echo $R1 | xargs basename | sed -E 's/_R*//'`
 {
 
 echo '#PBS -l walltime=02:00:00'
 echo '#PBS -l nodes=1:ppn=1'
 echo '#PBS -A PAS0471'
 echo '#PBS -l mem=4315MB'
 
 echo 'module use $HOME/local/share/lmodfiles'
 echo 'module load TrimGalore/0.6.5'
 echo 'module load Cutadapt/2.10'
 
 echo "trim_galore --paired --retain_unpaired --phred33 --illumina --path_to_cutadapt /users/PAS0471/osu9453/.local/bin/cutadapt \
 --output_dir /fs/project/PAS0471/osu9453/AgeCohortProject/trimgaloreOutPE \
 --length 20 -q 20 -e 0.1 \
 $R1 \
 $R2"
 
 } > trimgalore_$SAMPID.sh
 
 qsub trimgalore_$SAMPID.sh

done < listPE_R1.txt 3< listPE_R2.txt
