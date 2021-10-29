#!/bin/bash
#
# Runs jackhmmer search with bitscore thresholds
#
### #SBATCH --cluster=<clustername> 
#SBATCH --partition=pbatch
### #SBATCH --account=<accountname>
#SBATCH --job-name=jackhmmer
#SBATCH --output=jackhmmer.out
### #SBATCH --gres=gpu:0              # Number of GPU(s) per node.
#SBATCH --cpus-per-task=8         # CPU cores/threads
### #SBATCH --mem=48000M              # memory per node
#SBATCH --time=0-24:00            # Max time (DD-HH:MM)
#SBATCH --ntasks=1                # Only set to >1 if you want to use multi-threading

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## USAGE
## Create a directory (prot_dir), and put the WT sequence in the file wt.fa in that dir
## sbatch jackhmmer.sh <bin_dir> <prot_dir> <bitscore_threshold> <niter> <seqdb>

bindir=$1
protdir=$2
bitscore=$3   # e.g. 0.5
niter=$4      # 5 typically
seqdb=$5      # db location for e.g. uniref100 or uniref90 fasta files

jckhmm="$bindir/jackhmmer"
query="$protdir/wt.fa"
tblout="$protdir/targets.tblout"
alignmentfile="$protdir/alignment.sto"
hmmprefix="$protdir/iter"
aliprefix="$protdir/iter"

cpu=1
if [[ ! -z "$SLURM_CPUS_PER_TASK" ]]
then
    cpu="$SLURM_CPUS_PER_TASK"
fi

wtseq=$(sed 1d $query)
seqlen=${#wtseq}
bitscore=$(echo "$seqlen*$bitscore" | bc)   # scale bitscore by seqlen
echo "INFO: The scaled bitscore for jackhmmer is: $bitscore"

#EVcouplings defaults
$jckhmm -N $niter \
    --incT $bitscore --incdomT $bitscore -T $bitscore --domT $bitscore \
    --popen 0.02 --pextend 0.4 --mx BLOSUM62 \
    --tblout $tblout -A $alignmentfile --noali --notextw\
    --chkhmm $hmmprefix --chkali $aliprefix \
    --cpu $cpu \
    $query $seqdb
