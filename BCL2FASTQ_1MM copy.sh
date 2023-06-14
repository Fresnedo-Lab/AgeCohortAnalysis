#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=28
#PBS -A PAS0471
#PBS -l mem=128GB
 
module load singularity/current

 
cd /fs/project/PAS0471/osu9453/AgeCohortProject/Demultiplexing/L005_RawData/200714_K00400_0176_AHGVGVBBXY_Pearlly_L005
 
singularity exec ../../bcl2fastq-centos.sif bcl2fastq 