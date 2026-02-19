import os
import pandas as pd
import random

# Base directory for dummy project
BASE_DIR = "Deliverables_Dummy"

def create_structure(base_dir):
    structure = [
        "01_Raw_Data",
        "02_reference_genome_and_gff",
        "03_transcript_assembly_gtf",
        "04_transcript_sequences_fasta",
        "05_differential_expression_analysis",
        "06_Significant_DGE_GO",
        "07_Significant_DGE_pathways",
        "08_Significant_DGE_Enrichment"
    ]
    for folder in structure:
        os.makedirs(os.path.join(base_dir, folder), exist_ok=True)

def create_raw_data(base_dir, samples):
    path = os.path.join(base_dir, "01_Raw_Data")
    # Write fake fastq.gz to prompt sample detection
    stats_data = []
    mapping_data = []

    for sample in samples:
        open(os.path.join(path, f"{sample}_R1_001.fastq.gz"), 'w').close()
        open(os.path.join(path, f"{sample}_R2_001.fastq.gz"), 'w').close()
        
        # Raw Stats Data
        reads = random.randint(20000000, 30000000)
        bases = reads * 150 * 2
        gb = bases / 1e9
        stats_data.append({
            "Sample": sample.replace("Sample", "NGS_Dummy_Sample"), 
            "Total Reads in R1": reads,
            "Total Reads in R2": reads,
            "Total Reads(R1+R2)": reads * 2,
            "Total Bases (R1+R2)": bases,
            "Total Data(GB)": f"{gb:.2f}"
        })

        # Mapping Stats Data
        mapped = int(reads * 2 * 0.9) # 90% mapping
        unique = int(mapped * 0.8)    # 80% unique
        mapping_data.append({
            "Sample Name": sample.replace("Sample", "NGS_Dummy_Sample"),
            "Total Reads": reads * 2,
            "No. of Mapped reads": mapped,
            "% of mapped reads": "90.00",
            "# uniquely mapped reads": unique,
            "% uniquely mapped reads": "72.00"
        })

    # Create Raw Stats File (NGS_Dummy_Raw_Stats_with_GB.txt)
    df_stats = pd.DataFrame(stats_data)
    df_stats.to_csv(os.path.join(path, "NGS_Dummy_Raw_Stats_with_GB.txt"), sep='\t', index=False)

    # Create Mapping Stats File (Mapping.txt)
    df_mapping = pd.DataFrame(mapping_data)
    df_mapping.to_csv(os.path.join(path, "Mapping.txt"), sep='\t', index=False)

def create_ref_genome(base_dir):
    path = os.path.join(base_dir, "02_reference_genome_and_gff")
    with open(os.path.join(path, "dummy_genome.fa"), "w") as f:
        f.write(">chr1\n" + "A" * 1000 + "\n")
        f.write(">chr2\n" + "T" * 2000 + "\n")
    
    with open(os.path.join(path, "dummy.gff"), "w") as f:
        # Minimal GFF-like lines
        for i in range(10):
            f.write(f"chr1\t.\tgene\t{i*100}\t{i*100+50}\t.\t+\t.\tID=gene{i}\n")

def create_assembly_stats(base_dir, samples):
    # Mocking as if we scan GTF or Fasta
    path = os.path.join(base_dir, "04_transcript_sequences_fasta")
    with open(os.path.join(path, "Merged.fasta"), "w") as f:
        for i in range(100):
            f.write(f">transcript_{i}\n" + "AGCT" * random.randint(50, 200) + "\n")

    for sample in samples:
        with open(os.path.join(path, f"NGS_{sample}Aligned_transcript.fasta"), "w") as f:
            for i in range(50):
                 f.write(f">tr_{i}\n" + "AGCT" * random.randint(50, 100) + "\n")

def create_dge_data(base_dir):
    path = os.path.join(base_dir, "05_differential_expression_analysis")
    comparisons = ["Comp1_WT_vs_Mutatnt", "Comp2_Severe_vs_Mild"]
    
    for comp in comparisons:
        # Create Excel with 'logFC' and 'PValue' columns
        data = {
            "GeneID": [f"Gene_{i}" for i in range(100)],
            "logFC": [random.uniform(-5, 5) for _ in range(100)],
            "FDR": [random.uniform(0, 0.1) for _ in range(100)]
        }
        df = pd.DataFrame(data)
        df.to_excel(os.path.join(path, f"{comp}_DGE.xlsx"), index=False)

def create_pathway_data(base_dir):
    path = os.path.join(base_dir, "07_Significant_DGE_pathways")
    # Mock data for pathway counts
    data = {
        "Pathway": ["Glycolysis", "TCA Cycle", "Photosynthesis"],
        "Count": [15, 10, 25],
        "Category": ["Metabolism", "Metabolism", "Metabolism"]
    }
    df = pd.DataFrame(data)
    df.to_excel(os.path.join(path, "Comp1_Significant_DGE_Pathways.xlsx"), index=False)

def main():
    samples = ["SampleA", "SampleB", "SampleC", "SampleD"]
    create_structure(BASE_DIR)
    create_raw_data(BASE_DIR, samples)
    create_ref_genome(BASE_DIR)
    create_assembly_stats(BASE_DIR, samples)
    create_dge_data(BASE_DIR)
    create_pathway_data(BASE_DIR)
    print(f"Dummy data created in {BASE_DIR}")

if __name__ == "__main__":
    main()
