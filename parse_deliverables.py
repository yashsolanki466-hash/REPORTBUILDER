import os
import re
import json
import pandas as pd
from pathlib import Path
from collections import defaultdict


class DeliverablesParser:

    def __init__(self, root):
        self.root = Path(root)
        self.data = defaultdict(dict)

    # -----------------------------
    # PROJECT METADATA
    # -----------------------------
    def parse_readme(self):
        readme = self.root / "Readme.txt"

        if not readme.exists():
            return

        text = readme.read_text(errors="ignore")

        patterns = {
            "project_id": r"Project\s*ID\s*:\s*(.*)",
            "pi": r"Project\s*PI\s*:\s*(.*)",
            "application": r"Application\s*:\s*(.*)",
            "samples": r"No\s*of\s*Samples\s*:\s*(.*)"
        }

        for k, p in patterns.items():
            m = re.search(p, text, re.IGNORECASE)
            if m:
                self.data["project"][k] = m.group(1).strip()

    # -----------------------------
    # SAMPLE DETECTION
    # -----------------------------
    def parse_samples(self):

        raw = self.root / "01_Raw_Data"
        samples = set()

        for f in raw.glob("*.fastq.gz"):

            name = f.name
            m = re.match(r"(.*?)_S\d+_L\d+_R[12]", name)

            if m:
                samples.add(m.group(1))

        self.data["samples"] = sorted(samples)

    # -----------------------------
    # RAW STATS
    # -----------------------------
    def parse_raw_stats(self):

        raw = self.root / "01_Raw_Data"

        stats_file = None

        for f in raw.glob("*Raw_stats*"):
            stats_file = f

        if stats_file is None:
            return

        try:

            df = pd.read_csv(stats_file, sep="\t")

            rows = []

            for _, r in df.iterrows():

                rows.append({
                    "sample": str(r.get("Sample", "")),
                    "total_reads": int(r.get("Total Reads", 0)),
                    "total_bases": int(r.get("Total Bases", 0)),
                    "data_gb": round(int(r.get("Total Bases", 0)) / 1e9, 2)
                })

            self.data["raw_stats"] = rows

        except:
            pass

    # -----------------------------
    # REFERENCE GENOME
    # -----------------------------
    def parse_reference(self):

        ref = self.root / "02_reference_genome_and_gff"

        fna = list(ref.glob("*.fna"))
        gff = list(ref.glob("*.gff"))

        self.data["reference"] = {
            "genome_fasta": str(fna[0]) if fna else None,
            "genome_gff": str(gff[0]) if gff else None
        }

    # -----------------------------
    # TRANSCRIPT ASSEMBLY
    # -----------------------------
    def parse_gtf(self):

        gtf_dir = self.root / "03_transcript_assembly_gtf"

        merged = gtf_dir / "merged_transcripts.gtf"

        if not merged.exists():
            return

        transcripts = set()
        total_len = 0

        with open(merged) as f:

            for line in f:

                if line.startswith("#"):
                    continue

                cols = line.split("\t")

                if len(cols) < 9:
                    continue

                if cols[2] != "exon":
                    continue

                start = int(cols[3])
                end = int(cols[4])

                total_len += (end - start + 1)

                attr = cols[8]

                m = re.search(r'transcript_id "([^"]+)"', attr)

                if m:
                    transcripts.add(m.group(1))

        self.data["assembly"] = {
            "transcripts": len(transcripts),
            "total_length": total_len,
            "mean_length": int(total_len / max(len(transcripts), 1))
        }

    # -----------------------------
    # DIFFERENTIAL EXPRESSION
    # -----------------------------
    def parse_dge(self):

        dge_dir = self.root / "05_differential_expression_analysis"

        xlsx = list(dge_dir.glob("*DGE*.xlsx"))

        if not xlsx:
            return

        df = pd.read_excel(xlsx[0])

        up = 0
        down = 0

        if "log2FoldChange" in df.columns:

            up = (df["log2FoldChange"] > 0).sum()
            down = (df["log2FoldChange"] < 0).sum()

        self.data["dge"] = {
            "total": len(df),
            "upregulated": int(up),
            "downregulated": int(down)
        }

    # -----------------------------
    # GO RESULTS
    # -----------------------------
    def parse_go(self):

        go_dir = self.root / "06_Significant_DGE_GO"

        files = list(go_dir.glob("*.xlsx"))

        if not files:
            return

        df = pd.read_excel(files[0])

        self.data["go"] = {
            "terms": len(df)
        }

    # -----------------------------
    # KEGG
    # -----------------------------
    def parse_kegg(self):

        kegg = self.root / "07_Significant_DGE_pathways"

        res = list(kegg.glob("*Results*.xlsx"))

        imgs = list((kegg / "kegg_pathway_images").glob("*.png"))

        self.data["kegg"] = {
            "pathways": len(res),
            "images": len(imgs)
        }

    # -----------------------------
    # ENRICHMENT
    # -----------------------------
    def parse_enrichment(self):

        enr = self.root / "08_Significant_DGE_Enrichment"

        go = list(enr.glob("*GO_enrichment_results*.xlsx"))
        kegg = list(enr.glob("*KEGG_enrichment_results*.xlsx"))

        self.data["enrichment"] = {
            "go_files": len(go),
            "kegg_files": len(kegg)
        }

    # -----------------------------
    # PCA
    # -----------------------------
    def parse_pca(self):

        pca = self.root / "09_PCA_plot"

        files = list(pca.glob("*"))

        self.data["pca"] = [str(x) for x in files]

    # -----------------------------
    # PCOA
    # -----------------------------
    def parse_pcoa(self):

        pcoa = self.root / "10_PCOA_plot"

        files = list(pcoa.glob("*"))

        self.data["pcoa"] = [str(x) for x in files]

    # -----------------------------
    # RUN ALL
    # -----------------------------
    def run(self):

        self.parse_readme()
        self.parse_samples()
        self.parse_raw_stats()
        self.parse_reference()
        self.parse_gtf()
        self.parse_dge()
        self.parse_go()
        self.parse_kegg()
        self.parse_enrichment()
        self.parse_pca()
        self.parse_pcoa()

        return dict(self.data)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="report_data.json")

    args = parser.parse_args()

    p = DeliverablesParser(args.input)

    result = p.run()

    with open(args.output, "w") as f:
        json.dump(result, f, indent=4)

    print("Report data extracted to", args.output)