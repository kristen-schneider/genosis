#SBATCH -p short
#SBATCH --job-name=slurm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32gb
#SBATCH --time=23:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=/home/name/precision-medicine/workflow/out/slurm.out
#SBATCH --error=/home/name/precision-medicine/workflow/err/slurm.err

sbatch SLURM_slice-encode.sh
sbatch SLURM_embed.sh
sbatch SLURM_index.sh
sbatch SLURM_search-aggregate-evaluate.sh
