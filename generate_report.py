import os
import argparse
import pandas as pd
from jinja2 import Environment, FileSystemLoader
import glob
import re
from datetime import datetime


def _choose_file(files):
    if not files:
        return None
    existing = [f for f in files if os.path.exists(f)]
    if not existing:
        return None
    existing.sort(key=lambda p: (os.path.getmtime(p), os.path.getsize(p)), reverse=True)
    return existing[0]


def _find_first_existing(input_dir, relative_globs):
    candidates = []
    for rel_glob in relative_globs:
        candidates.extend(glob.glob(os.path.join(input_dir, rel_glob)))
    return _choose_file(candidates)


def _parse_project_readme(input_dir):
    readme_path = _find_first_existing(input_dir, ["Readme.txt", "README.txt", "readme.txt"])
    details = {}
    if not readme_path:
        return details

    try:
        with open(readme_path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()

        patterns = {
            "project_id": r"Project\s*ID\s*:\s*(.+)",
            "project_pi": r"Project\s*PI\s*:\s*(.+)",
            "application": r"Application\s*:\s*(.+)",
            "no_of_samples": r"No\s*of\s*Samples\s*:\s*(.+)",
        }
        for key, pat in patterns.items():
            m = re.search(pat, content, flags=re.IGNORECASE)
            if m:
                details[key] = m.group(1).strip()
    except Exception:
        return details

    return details


def _validate_deliverables_structure(input_dir):
    warnings = []
    required_dirs = [
        "01_Raw_Data",
        "02_reference_genome_and_gff",
        "03_transcript_assembly_gtf",
        "04_transcript_sequences_fasta",
        "05_differential_expression_analysis",
    ]
    optional_dirs = [
        "06_Significant_DGE_GO",
        "07_Significant_DGE_pathways",
        "08_Significant_DGE_Enrichment",
    ]

    for d in required_dirs:
        if not os.path.isdir(os.path.join(input_dir, d)):
            warnings.append(f"Missing required directory: {d}")
    for d in optional_dirs:
        if not os.path.isdir(os.path.join(input_dir, d)):
            warnings.append(f"Missing optional directory: {d}")

    raw_dir = os.path.join(input_dir, "01_Raw_Data")
    if os.path.isdir(raw_dir):
        fastqs = glob.glob(os.path.join(raw_dir, "*.fastq.gz"))
        if not fastqs:
            warnings.append("No .fastq.gz files found in 01_Raw_Data (sample detection may rely on Mapping/Stats).")

    return warnings

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate NGS Report")
    parser.add_argument("--input", required=True, help="Path to the project directory")
    parser.add_argument("--output", default="report.html", help="Output HTML file path")
    parser.add_argument("--strict", action="store_true", help="Fail if required inputs are missing")
    parser.add_argument("--client-name", default=None, help="Client name to display on the cover page")
    parser.add_argument("--client-org", default=None, help="Client organization to display on the cover page")
    parser.add_argument("--reference-organism", default=None, help="Reference organism name used in the report")
    parser.add_argument("--logo-path", default=None, help="Logo path to use in the report")
    return parser.parse_args()

def get_project_details(input_dir):
    # Try to infer project ID from directory name or a config
    project_id = os.path.basename(os.path.normpath(input_dir))
    parsed = _parse_project_readme(input_dir)
    if parsed.get("project_id"):
        project_id = parsed["project_id"]
    return {
        "project_id": project_id,
        "report_date": datetime.now().strftime("%d-%m-%Y"),
        "client_name": "Client Name", # Placeholder
        "client_org": "Client Organization", # Placeholder
        "project_pi": parsed.get("project_pi", ""),
        "application": parsed.get("application", ""),
        "no_of_samples": parsed.get("no_of_samples", "")
    }

def get_samples(input_dir):
    raw_data_dir = os.path.join(input_dir, "01_Raw_Data")
    if not os.path.exists(raw_data_dir):
        return []
    
    samples = set()
    
    # Priority 1: Try to look for fastq files
    files = glob.glob(os.path.join(raw_data_dir, "*.fastq.gz"))
    for f in files:
        basename = os.path.basename(f)
        match = re.match(r"(.*)_R[12].*\.fastq\.gz", basename)
        if match:
            samples.add(match.group(1))

    # Priority 2: If no fastq, extract from Mapping file (txt or xlsx)
    if not samples:
        mapping_files = glob.glob(os.path.join(raw_data_dir, "Mapping.*"))
        mapping_files = sorted(mapping_files)
        for mapping_file in mapping_files:
            try:
                if mapping_file.endswith('.xlsx'):
                    df = pd.read_excel(mapping_file)
                else:
                    df = pd.read_csv(mapping_file, sep=None, engine='python')
                
                df.columns = df.columns.str.strip()
                sample_col = next((c for c in df.columns if 'Sample' in c), None)
                if sample_col:
                    samples.update(df[sample_col].dropna().astype(str).unique())
            except Exception as e:
                print(f"Error reading {mapping_file} in get_samples: {e}")
    
    if not samples:
        stats_files = glob.glob(os.path.join(raw_data_dir, "*Raw*Stats*"))
        # Also check for Rawstats.xlsx specifically if glob didn't catch it
        if not stats_files:
            stats_files = glob.glob(os.path.join(raw_data_dir, "*Rawstats*"))
            
        for f in stats_files:
            try:
                if f.endswith('.xlsx'):
                    df = pd.read_excel(f)
                else:
                    df = pd.read_csv(f, sep="\t")
                
                df.columns = df.columns.str.strip()
                sample_col = next((c for c in df.columns if 'Sample' in c), None) # 'Sample' or 'SampleID'
                if sample_col:
                    samples.update(df[sample_col].dropna().astype(str).unique())
            except Exception as e:
                print(f"Error reading Stats file {f} in get_samples: {e}")

    return sorted(list(samples))

def get_sequencing_stats(input_dir, samples):
    # Parse *Raw_Stats_with_GB.txt or .xlsx in 01_Raw_Data
    stats = []
    raw_data_dir = os.path.join(input_dir, "01_Raw_Data")
    
    # Try txt then fallback to xlsx or wildcard
    stats_files = glob.glob(os.path.join(raw_data_dir, "*Raw*Stats*"))
    if not stats_files:
        stats_files = glob.glob(os.path.join(raw_data_dir, "*Rawstats*"))
    stats_files = sorted(stats_files)
    
    if stats_files:
        try:
            f = _choose_file(stats_files)
            if f.endswith('.xlsx'):
                df = pd.read_excel(f)
            else:
                df = pd.read_csv(f, sep="\t")

            # Normalize columns
            df.columns = df.columns.str.strip()
            
            for _, row in df.iterrows():
                # Flexible matching
                sample_col = next((c for c in df.columns if 'sample' in c.lower()), None)
                # Specific columns from user's file: "Total Reads(R1+R2)", "Total Bases (R1+R2)", "Total Data(GB)"
                reads_col = next((c for c in df.columns if '(R1+R2)' in c and 'Reads' in c), None) 
                if not reads_col: reads_col = next((c for c in df.columns if 'ReadCount' in c), None)
                
                bases_col = next((c for c in df.columns if '(R1+R2)' in c and 'Bases' in c), None)
                if not bases_col: bases_col = next((c for c in df.columns if 'BaseCount' in c), None)

                gb_col = next((c for c in df.columns if 'Total Data' in c), None)

                if sample_col:
                    reads_val = str(row[reads_col]).replace(',', '') if reads_col and pd.notna(row[reads_col]) else "0"
                    bases_val = str(row[bases_col]).replace(',', '') if bases_col and pd.notna(row[bases_col]) else "0"
                    gb_val = str(row[gb_col]).replace(',', '') if gb_col and pd.notna(row[gb_col]) else "0"

                    try:
                        reads_fmt = f"{int(float(reads_val)):,}"
                    except: reads_fmt = "N/A"
                    
                    try:
                        bases_fmt = f"{int(float(bases_val)):,}"
                    except: bases_fmt = "N/A"

                    stats.append({
                        "sample": str(row[sample_col]),
                        "total_reads": reads_fmt,
                        "total_bases": bases_fmt,
                        "data_gb": gb_val
                    })
        except Exception as e:
            print(f"Error parsing raw stats file: {e}")
    
    # If no stats found
    if not stats:
        for s in samples:
            stats.append({
                "sample": s,
                "total_reads": "N/A",
                "total_bases": "N/A",
                "data_gb": "N/A"
            })
    return stats

def get_mapping_stats(input_dir):
    raw_data_dir = os.path.join(input_dir, "01_Raw_Data")
    mapping_stats = []
    mapping_files = glob.glob(os.path.join(raw_data_dir, "Mapping.*"))
    mapping_files = sorted(mapping_files)
    
    if mapping_files:
        try:
            f = _choose_file(mapping_files)
            # Fix: sep=None with engine='python' allows sniffing, 
            # and avoid passing regex characters that conflict
            if f.endswith('.xlsx'):
                df = pd.read_excel(f)
            else:
                df = pd.read_csv(f, sep=r'\s+|\t', engine='python')
                
            df.columns = df.columns.str.strip()

            for _, row in df.iterrows():
                # Columns: Sample Name, Total Reads, No. of Mapped reads, % of mapped reads, # uniquely mapped reads, % uniquely mapped reads
                sample_col = next((c for c in df.columns if 'Sample' in c), None)
                mapped_col = next((c for c in df.columns if 'No. of Mapped reads' in c), None)
                pct_mapped_col = next((c for c in df.columns if '% of mapped reads' in c), None)
                unique_col = next((c for c in df.columns if 'uniquely mapped reads' in c and '%' not in c), None) # Avoid % col
                pct_unique_col = next((c for c in df.columns if '% uniquely mapped reads' in c or '% of uniquely mapped reads' in c), None)

                total_reads_col = next((c for c in df.columns if 'Total Reads' in c), None)

                if sample_col:
                     mapping_stats.append({
                        "sample_name": row[sample_col],
                        "total_reads": f"{row[total_reads_col]:,}" if total_reads_col and pd.notna(row[total_reads_col]) else "N/A",
                        "mapped_reads": f"{row[mapped_col]:,}" if mapped_col and pd.notna(row[mapped_col]) else "N/A",
                        "percent_mapped": str(row[pct_mapped_col]) if pct_mapped_col and pd.notna(row[pct_mapped_col]) else "N/A",
                        "unique_reads": f"{row[unique_col]:,}" if unique_col and pd.notna(row[unique_col]) else "N/A",
                        "percent_unique": str(row[pct_unique_col]) if pct_unique_col and pd.notna(row[pct_unique_col]) else "N/A"
                    })
        except Exception as e:
            print(f"Error parsing Mapping.txt: {e}")
            
    return mapping_stats

def get_ref_stats(input_dir):
    ref_dir = os.path.join(input_dir, "02_reference_genome_and_gff")
    stats = {
        "total_scaffolds": 0,
        "genome_length": 0,
        "mean_scaffold_size": 0,
        "max_scaffold_size": 0
    }
    total_genes = 0
    
    if os.path.exists(ref_dir):
        # Parse Fasta for genome stats
        fasta_files = glob.glob(os.path.join(ref_dir, "*.fa")) + glob.glob(os.path.join(ref_dir, "*.fasta"))
        if fasta_files:
            scaffold_lengths = []
            fasta_path = _choose_file(fasta_files)
            with open(fasta_path, 'r') as f:
                curr_len = 0
                for line in f:
                    if line.startswith(">"):
                        if curr_len > 0:
                            scaffold_lengths.append(curr_len)
                        curr_len = 0
                    else:
                        curr_len += len(line.strip())
                if curr_len > 0:
                    scaffold_lengths.append(curr_len)
            
            if scaffold_lengths:
                stats["total_scaffolds"] = len(scaffold_lengths)
                stats["genome_length"] = f"{sum(scaffold_lengths):,}"
                stats["mean_scaffold_size"] = f"{int(sum(scaffold_lengths)/len(scaffold_lengths)):,}"
                stats["max_scaffold_size"] = f"{max(scaffold_lengths):,}"

        # Parse GFF for gene count
        gff_files = glob.glob(os.path.join(ref_dir, "*.gff")) + glob.glob(os.path.join(ref_dir, "*.gff3"))
        if gff_files:
            gff_path = _choose_file(gff_files)
            with open(gff_path, 'r') as f:
                for line in f:
                    if "\tgene\t" in line:
                        total_genes += 1
    
    return stats, f"{total_genes:,}"

def get_assembly_stats(input_dir, samples):
    assembly_stats = []
    # Logic: Read 04_transcript_sequences_fasta for merged and sample-wise fastas
    fasta_dir = os.path.join(input_dir, "04_transcript_sequences_fasta")
    
    def calc_fasta_stats(filepath):
        lengths = []
        with open(filepath, 'r') as f:
            curr = 0
            for line in f:
                if line.startswith(">"):
                    if curr > 0: lengths.append(curr)
                    curr = 0
                else:
                    curr += len(line.strip())
            if curr > 0: lengths.append(curr)
        
        if not lengths: return 0, 0, 0
        return len(lengths), sum(lengths), int(sum(lengths)/len(lengths))

    if os.path.exists(fasta_dir):
        # Merged
        merged_files = glob.glob(os.path.join(fasta_dir, "Merged.fasta"))
        if merged_files:
             merged_path = _choose_file(merged_files)
             count, total, mean = calc_fasta_stats(merged_path)
             assembly_stats.append({
                 "sample": "merged.fasta",
                 "num_transcripts": f"{count:,}",
                 "total_bp": f"{total:,}",
                 "mean_size": f"{mean:,}"
             })
        
        # Per Sample (trying to match sample names)
        for s in samples:
             # Try to find a file that looks like it belongs to the sample
             # Pattern from user: NGS_240540_PM-CONTROL-2Aligned_transcript.fasta
             # It's fuzzy, so we look for files containing the sample name
             files = glob.glob(os.path.join(fasta_dir, f"*{s}*.fasta"))
             for f in files:
                if "Merged" in f: continue
                count, total, mean = calc_fasta_stats(f)
                assembly_stats.append({
                 "sample": os.path.basename(f),
                 "num_transcripts": f"{count:,}",
                 "total_bp": f"{total:,}",
                 "mean_size": f"{mean:,}"
             })

    return assembly_stats

def get_dge_stats(input_dir):
    dge_dir = os.path.join(input_dir, "05_differential_expression_analysis")
    stats = []
    labels = []
    up_counts = []
    down_counts = []
    
    if os.path.exists(dge_dir):
        files = glob.glob(os.path.join(dge_dir, "*DGE.xlsx"))
        for f in files:
            try:
                df = pd.read_excel(f)
                # Assuming standard edgeR output or similar with logFC and FDR/PValue
                # Adjust column names as per actual files
                # User's file implies: Comparison1_DGE.xlsx
                comparison_name = os.path.basename(f).replace("_DGE.xlsx", "")
                
                # Check for likely column names
                fc_col = next((c for c in df.columns if "logFC" in c or "log2FoldChange" in c), None)
                p_col = next((c for c in df.columns if "FDR" in c or "adj.P.Val" in c or "padj" in c), None)
                
                if fc_col and p_col:
                    sig_up = len(df[(df[fc_col] > 1) & (df[p_col] < 0.05)])
                    sig_down = len(df[(df[fc_col] < -1) & (df[p_col] < 0.05)])
                    total_sig = sig_up + sig_down
                    
                    stats.append({
                        "comparison": comparison_name,
                        "total_degs": len(df), # Or maybe total genes tested? usually DGE file has all genes
                        "sig_down": sig_down,
                        "sig_up": sig_up,
                        "total_sig": total_sig
                    })
                    
                    labels.append(comparison_name)
                    up_counts.append(sig_up)
                    down_counts.append(sig_down)
            except Exception as e:
                print(f"Error reading DGE file {f}: {e}")

    return stats, labels, up_counts, down_counts

def get_pathway_stats(input_dir):
    # This usually parses the KEGG/GO result files
    # User's path: 07_Significant_DGE_pathways/Comparison1_Significant_DGE_Pathways.xlsx
    path_dir = os.path.join(input_dir, "07_Significant_DGE_pathways")
    stats = []
    if os.path.exists(path_dir):
        files = glob.glob(os.path.join(path_dir, "*.xlsx"))
        if files:
            try:
                f = _choose_file(files)
                df = pd.read_excel(f)
                # Expecting columns like Pathway, Count, Category
                # Mapping likely columns
                cat_col = next((c for c in df.columns if "Category" in c or "Class" in c), "Category")
                sub_col = next((c for c in df.columns if "Sub" in c or "Pathway" in c), "Pathway")
                count_col = next((c for c in df.columns if "Count" in c or "Gene" in c), "Count")

                # Aggregate if needed, or just take top N
                # Here we just take the first 10 rows
                for _, row in df.head(10).iterrows():
                    stats.append({
                        "category": row.get(cat_col, "N/A"),
                        "subcategory": row.get(sub_col, "N/A"),
                        "counts": row.get(count_col, 0)
                    })
            except Exception as e:
                 print(f"Error reading pathway file: {e}")
    return stats

def main():
    args = parse_arguments()
    input_dir = args.input
    output_file = args.output
    
    print(f"Processing project: {input_dir}")
    
    warnings = _validate_deliverables_structure(input_dir)
    if args.strict and any(w.startswith("Missing required directory") for w in warnings):
        raise SystemExit("\n".join(warnings))

    # 1. Get Project Details
    project_details = get_project_details(input_dir)
    
    # 2. Get Samples
    samples = get_samples(input_dir)
    print(f"Found {len(samples)} samples: {samples}")
    
    # 3. Get Sequencing Stats
    seq_stats = get_sequencing_stats(input_dir, samples)
    
    # 4. Reference Genome Stats
    ref_stats, total_genes = get_ref_stats(input_dir)
    
    # 5. Mapping Stats
    mapping_stats = get_mapping_stats(input_dir)
    # Fallback if empty logic (optional)
    if not mapping_stats:
        # Try to infer or leave empty
        for s in samples:
            mapping_stats.append({
                "sample_name": s,
                "mapped_reads": "N/A",
                "percent_mapped": "N/A",
                "unique_reads": "N/A",
                "percent_unique": "N/A"
            })

    # 6. Assembly Stats
    assembly_stats = get_assembly_stats(input_dir, samples)
    
    # 7. DGE Stats
    dge_stats, dge_labels, dge_up, dge_down = get_dge_stats(input_dir)
    
    # 8. Component Stats (Pathway)
    pathway_stats = get_pathway_stats(input_dir)

    client_name = args.client_name if args.client_name is not None else project_details["client_name"]
    client_org = args.client_org if args.client_org is not None else project_details["client_org"]
    reference_organism = args.reference_organism if args.reference_organism is not None else "Organism Name"
    logo_path = args.logo_path if args.logo_path is not None else "assets/logo.png"

    # Prepare Context
    context = {
        "project_id": project_details["project_id"],
        "report_date": project_details["report_date"],
        "client_name": client_name,
        "client_org": client_org,
        "project_pi": project_details.get("project_pi", ""),
        "application": project_details.get("application", ""),
        "no_of_samples": project_details.get("no_of_samples", ""),
        "sample_count": len(samples),
        "samples": samples,
        "qubit_data": [], # Placeholder
        "library_sizes": [450] * len(samples), # Placeholder
        "sequencing_stats": seq_stats,
        "reference_organism": reference_organism,
        "total_genes": total_genes,
        "ref_stats": ref_stats,
        "mapping_stats": mapping_stats,
        "total_transcripts": sum(int(x['num_transcripts'].replace(',', '')) for x in assembly_stats if x['sample'] == 'merged.fasta') if assembly_stats else 0,
        "mean_transcript_size": next((x['mean_size'] for x in assembly_stats if x['sample'] == 'merged.fasta'), 0),
        "assembly_stats": assembly_stats,
        "diff_expr_stats": dge_stats,
        "dge_chart_labels": dge_labels,
        "dge_chart_up": dge_up,
        "dge_chart_down": dge_down,
        "pathway_stats": pathway_stats,
        "pathway_image_src": "", # Logic to find image if exists
        "logo_path": logo_path,
        "warnings": warnings
    }
    
    # Render
    script_dir = os.path.dirname(os.path.abspath(__file__))
    templates_dir = os.path.join(script_dir, 'templates')
    env = Environment(loader=FileSystemLoader(templates_dir))
    template = env.get_template('report_template.html')
    
    with open(output_file, 'w') as f:
        f.write(template.render(context))
        
    print(f"Report generated: {output_file}")

if __name__ == "__main__":
    main()
