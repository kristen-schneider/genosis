#SBATCH -p short
#SBATCH --job-name=index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32gb
#SBATCH --time=02:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=/home/name/precision-medicine/workflow/out/index.out
#SBATCH --error=/home/name/precision-medicine/workflow/err/index.err

pmed_dir="/home/name/precision-medicine/"
out_dir="/home/name/out_dir/"
config="/home/name/out_dir/config.yml"
singularity="/home/name/singularity.sif"

singularity run $singularity \
	bash $pmed_dir"index.sh" \
		$pmed_dir \
		$out_dir \
		$config \
