#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=12
#PBS -A PAS0471
#PBS -l mem=48GB

cd /fs/project/PAS0471/osu9453/AgeCohortProject

module load bowtie2/2.4.1

software/Bismark-0.22.3/bismark_genome_preparation --verbose --path_to_aligner /usr/local/bowtie2/2.4.1 \
--parallel 6 /fs/project/PAS0471/osu9453/AgeCohortProject/genome
