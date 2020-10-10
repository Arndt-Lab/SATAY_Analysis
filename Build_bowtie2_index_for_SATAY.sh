#initiate srun session
srun -n1 -t04:00:00 --cpus-per-task=4 --mem=64g --pty bash

#MODULES
module load bowtie2/2.4.1

yeast_chromosomes_dir=/bgfs/karndt/Annotations/Saccharomyces_cerevisiae_Ensembl_R64-1-1/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/Chromosomes
satay_plasmid_dir=/bgfs/karndt/Annotations

#navigate to the appropriate directory
cd /bgfs/karndt

#make a directory for the index
mkdir satay_index
#save path
satay_index=/bgfs/karndt/satay_index


#BUILD BOWTIE2 INDEX
bowtie2-build --threads 4 -f ${yeast_chromosomes_dir}/I.fa,${yeast_chromosomes_dir}/II.fa,${yeast_chromosomes_dir}/III.fa,${yeast_chromosomes_dir}/IV.fa,${yeast_chromosomes_dir}/IX.fa,${yeast_chromosomes_dir}/MT.fa,${yeast_chromosomes_dir}/V.fa,${yeast_chromosomes_dir}/VI.fa,${yeast_chromosomes_dir}/VII.fa,${yeast_chromosomes_dir}/VIII.fa,${yeast_chromosomes_dir}/X.fa,${yeast_chromosomes_dir}/XI.fa,${yeast_chromosomes_dir}/XII.fa,${yeast_chromosomes_dir}/XIII.fa,${yeast_chromosomes_dir}/XIV.fa,${yeast_chromosomes_dir}/XV.fa,${yeast_chromosomes_dir}/XVI.fa,${satay_plasmid_dir}/SATAY_plasmid.fa ${satay_index}/satay_index > ${satay_index}/build_info.txt

