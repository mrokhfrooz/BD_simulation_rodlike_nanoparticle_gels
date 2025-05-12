#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=35
#SBATCH --mem=14G
#SBATCH --job-name=LJ3_AR3
#SBATCH --account=???
#SBATCH --mail-type=ALL
#SBATCH --mail-user=???@gmail.com
#SBATCH --output=??_output.txt
#SBATCH --error=??_error.txt

################################################################################

module load gcc

cd $SLURM_SUBMIT_DIR

# Run the job times in parallel
for i in {1..40}; do
    # Create a unique folder for each run using a timestamp
    timestamp=$(date +%Y%m%d%H%M%S)
    run_folder="demo_run_${timestamp}_${i}"

    mkdir -p "$run_folder"
    cd "$run_folder"

    # Run your command in the background
    ../output_filename &

    cd $SLURM_SUBMIT_DIR  # Return to the original working directory
done

# Wait for all background processes to finish
wait