import subprocess

# Function to download and rename FASTQ files using fastq-dump
def download_and_rename_fastq(exp_name, sra_id):
    # Command to download FASTQ files in gzip format and split paired-end files
    fastq_dump_cmd = f"fastq-dump --gzip --split-files {sra_id}"
    subprocess.run(fastq_dump_cmd, shell=True)

    # Rename the downloaded FASTQ files
    rename_r1_cmd = f"mv {sra_id}_1.fastq.gz {sra_id}.{exp_name}_R1.fastq.gz"
    rename_r2_cmd = f"mv {sra_id}_2.fastq.gz {sra_id}.{exp_name}_R2.fastq.gz"
    subprocess.run(rename_r1_cmd, shell=True)
    subprocess.run(rename_r2_cmd, shell=True)

# Read table.txt file (assuming each line is "ExperimentName SRA_ID")
with open("sraids.txt", "r") as f:
    lines = f.readlines()

# Loop through each line in the table and download + rename FASTQ files
for line in lines:
    exp_name, sra_id = line.strip().split(" ")
    print(f"Downloading and renaming files for experiment {exp_name} with SRA ID {sra_id}")
    download_and_rename_fastq(exp_name, sra_id)

