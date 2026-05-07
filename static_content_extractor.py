import json
import os
import re
from datetime import datetime, timezone


_DEFAULT_CONSTANTS_TS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "constants.ts")


def _extract_default_template_html(constants_ts_text: str) -> str:
    m = re.search(r"export\s+const\s+DEFAULT_TEMPLATE\s*=\s*`", constants_ts_text)
    if not m:
        raise ValueError("DEFAULT_TEMPLATE start not found")

    start = m.end()
    end = constants_ts_text.rfind("`;")
    if end == -1 or end <= start:
        raise ValueError("DEFAULT_TEMPLATE end not found")

    return constants_ts_text[start:end]


def _strip_tags(s: str) -> str:
    s = re.sub(r"<\s*br\s*/?\s*>", "\n", s, flags=re.IGNORECASE)
    s = re.sub(r"<[^>]+>", "", s)
    return re.sub(r"\s+", " ", s).strip()


def _contains_placeholder(s: str) -> bool:
    return "{{" in s or "}}" in s


def _extract_tables(section_html: str):
    tables = []
    for t in re.finditer(r"<table(?P<attrs>[^>]*)>(?P<body>.*?)</table>", section_html, flags=re.DOTALL | re.IGNORECASE):
        table_html = t.group(0)
        if _contains_placeholder(table_html):
            continue

        headers = []
        thead = re.search(r"<thead>(.*?)</thead>", table_html, flags=re.DOTALL | re.IGNORECASE)
        if thead:
            headers = [_strip_tags(x) for x in re.findall(r"<th[^>]*>(.*?)</th>", thead.group(1), flags=re.DOTALL | re.IGNORECASE)]

        rows = []
        tbody = re.search(r"<tbody>(.*?)</tbody>", table_html, flags=re.DOTALL | re.IGNORECASE)
        body_region = tbody.group(1) if tbody else table_html
        for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", body_region, flags=re.DOTALL | re.IGNORECASE):
            cols = [_strip_tags(td) for td in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", tr, flags=re.DOTALL | re.IGNORECASE)]
            if cols and any(c for c in cols):
                rows.append(cols)

        if headers or rows:
            tables.append({"headers": headers, "rows": rows})

    return tables


def _extract_static_text_blocks(section_html: str):
    blocks = []
    for tag in ["p", "h2", "h3", "h4", "h5", "li"]:
        for inner in re.findall(rf"<{tag}[^>]*>(.*?)</{tag}>", section_html, flags=re.DOTALL | re.IGNORECASE):
            if _contains_placeholder(inner):
                continue
            text = _strip_tags(inner)
            if not text:
                continue
            # Skip obvious UI-only artifacts.
            if text.lower() in {"png", "pdf"}:
                continue
            blocks.append({"tag": tag.lower(), "text": text})

    # de-dup while preserving order
    seen = set()
    uniq = []
    for b in blocks:
        k = (b["tag"], b["text"])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(b)
    return uniq


def extract_static_content_from_constants(constants_ts_path: str = _DEFAULT_CONSTANTS_TS) -> dict:
    with open(constants_ts_path, "r", encoding="utf-8", errors="ignore") as f:
        constants_ts_text = f.read()

    html = _extract_default_template_html(constants_ts_text)

    sections = {}
    for m in re.finditer(
        r"<section\s+[^>]*id=\"(?P<id>[^\"]+)\"[^>]*>(?P<body>.*?)</section>",
        html,
        flags=re.DOTALL | re.IGNORECASE,
    ):
        section_id = m.group("id").strip()
        body = m.group("body")

        title = ""
        subtitle = ""
        mt = re.search(r"<h1[^>]*class=\"page-title\"[^>]*>(.*?)</h1>", body, flags=re.DOTALL | re.IGNORECASE)
        if mt and not _contains_placeholder(mt.group(1)):
            title = _strip_tags(mt.group(1))

        ms = re.search(r"<p[^>]*class=\"page-subtitle\"[^>]*>(.*?)</p>", body, flags=re.DOTALL | re.IGNORECASE)
        if ms and not _contains_placeholder(ms.group(1)):
            subtitle = _strip_tags(ms.group(1))

        sections[section_id] = {
            "title": title,
            "subtitle": subtitle,
            "text_blocks": _extract_static_text_blocks(body),
            "tables": _extract_tables(body),
        }

    snippets = {
        "wet_lab": {},
        "mapping": {},
        "reference": {},
        "assembly": {},
        "dge": {},
        "go": {},
        "workflow": {},
        "bioinformatics": {},
        "data_stats": {},
    }

    def _find_first_text(section_id: str, startswith: str) -> str:
        sec = sections.get(section_id, {})
        for b in sec.get("text_blocks", []):
            txt = b.get("text", "")
            if txt.startswith(startswith):
                return txt
        return ""

    def _find_first_text_any(startswith: str) -> str:
        for sec in sections.values():
            for b in sec.get("text_blocks", []):
                txt = b.get("text", "")
                if txt.startswith(startswith):
                    return txt
        return ""

    snippets["wet_lab"]["rna_isolation_qc"] = _find_first_text(
        "wet-lab",
        "Total RNA was extracted",
    )
    snippets["wet_lab"]["library_preparation"] = _find_first_text(
        "wet-lab",
        "Sequencing libraries were prepared",
    )

    snippets["mapping"]["methodology"] = _find_first_text(
        "mapping",
        "Indexing of the reference genome was carried out",
    )

    snippets["reference"]["genome_info"] = _find_first_text(
        "reference",
        "Based on information received from client",
    )
    snippets["reference"]["gff_info"] = _find_first_text(
        "reference",
        "The reference GFF file contains genome annotation",
    )

    snippets["assembly"]["stringtie_methodology_1"] = _find_first_text(
        "assembly",
        "StringTie assembles transcripts",
    )
    snippets["assembly"]["stringtie_methodology_2"] = _find_first_text(
        "assembly",
        "After assembling each sample",
    )
    snippets["assembly"]["output_notes"] = _find_first_text(
        "assembly",
        "The output consists of assembled gene/transcript GTF file",
    )
    snippets["assembly"]["gffcompare_method"] = _find_first_text(
        "assembly",
        "For comparison of assembled transcripts with reference transcripts",
    )
    snippets["assembly"]["mstrg_note"] = _find_first_text(
        "assembly",
        "\u201cMSTRG",
    )

    wf_tables = sections.get("workflow", {}).get("tables", [])
    if wf_tables:
        snippets["workflow"]["software_table"] = wf_tables[0]
    else:
        # Hardcoded software table from original report
        snippets["workflow"]["software_table"] = {
            "headers": ["Software", "Version", "Application"],
            "rows": [
                ["Trim Galore", "0.6.4", "Adapter and low-quality data removal"],
                ["STAR", "2.7.10", "Mapping of reads to reference genome"],
                ["Stringtie", "2.2.1", "Reference based assembly"],
                ["gffcompare", "0.12.6", "Comparing reference genome gff with merged gff file generated from stringtie"],
                ["gffread", "0.12.7", "Fetching fasta sequences corresponding to gff file generated after reference-based assembly by stringtie"],
                ["Blastx", "2.30.0+", "Similarity search against NCBI NR database"],
                ["Blast2go cli", "1.5", "GO mapping and annotation"],
                ["DESeq2", "1.44.0", "Differential expression analysis and its visualization"],
                ["KEGG KAAS", "Webserver", "Pathway analysis"],
                ["clusterProfiler", "4.12.6", "Over-representation (or enrichment) analysis"],
            ]
        }

    snippets["bioinformatics"]["analysis_intro"] = _find_first_text_any(
        "Bioinformatics analysis includes"
    )
    
    # Hardcoded workflow narrative from original report
    workflow_text = (
        "Workflow overview: Raw data was filtered to remove adapters and low-quality bases. Clean reads thus "
        "produced were used for analysis. Reference guided transcript assembly was performed for samples first "
        "by mapping clean reads on reference genome using STAR (v2.7.10a) aligner and then performing "
        "transcript assembly by StringTie (v 2.2.1). A consensus set of transcripts was obtained using SringTie "
        "merge function which merges together all the gene structures found in any of the samples. Transcript "
        "abundance was then estimated for individual sample using merged transcript consensus again using "
        "StringTie and read counts thus obtained for each transcript were taken as input for differential expression "
        "analysis using DESeq2 package. Gene Ontology and pathway analysis of the significantly differentially "
        "expressed transcripts were performed using Blast2go cli and KEGG-KAAS server respectively. Overall "
        "bioinformatics workflow is graphically represented in figure above and software version along with its "
        "use in table below."
    )
    snippets["bioinformatics"]["workflow_narrative"] = workflow_text

    snippets["data_stats"]["intro"] = _find_first_text_any(
        "Next-generation sequencing statistics"
    )
    # Hardcoded data stats intro from original report
    data_stats_text = (
        "The next generation sequencing for these samples were performed on the Illumina Novaseq X plus "
        "platform. Statistics of raw data is shown in Table below:"
    )
    snippets["data_stats"]["illumina_intro"] = data_stats_text

    # Hardcoded reference genome intro from original report
    ref_genome_text = (
        "The reference genome of Solanum lycopersicum (tomato) and its corresponding GFF was downloaded "
        "from NCBI (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/512/215/GCF_036512215.1_SLM_r2.1/)."
    )
    snippets["reference"]["genome_source"] = ref_genome_text

    # Hardcoded GFF info from original report
    gff_info_text = (
        "The reference GFF file contains genome annotation. It includes information (locus and description) "
        "regarding genes and proteins. There are a total of 40021 genes and 44391 mRNAs as per the GFF file "
        "of NCBI reference genome. The fasta sequence of the reference genome along with the corresponding "
        "GFF file is provided with the deliverables in the folder named '02_reference_genome_and_gff' for "
        "client's reference."
    )
    snippets["reference"]["gff_annotation"] = gff_info_text

    # Hardcoded mapping methodology from original report
    mapping_text = (
        "Indexing of the reference genome was carried out using STAR genomeGenerate mode option. Then the "
        "input reads, in FASTQ format, along with the indexed reference genome generated in previous step was "
        "given to the STAR aligner. The analysis begins by mapping HQ reads against the reference genome to "
        "identify the positions from where the reads originated. This mapping information allows us to collect "
        "subsets of the reads corresponding to each gene, and then to assemble and quantify transcripts represented "
        "by those reads. Hence, the high-quality reads were mapped to reference genome using STAR aligner to "
        "create alignment in BAM format for each sample with default parameters. The mapping statistics are "
        "provided in the table below:"
    )
    snippets["mapping"]["star_methodology"] = mapping_text

    # Hardcoded transcript assembly methodology
    stringtie_text = (
        "StringTie assembles transcripts from RNA-seq reads that have been aligned to the genome, first grouping "
        "the reads into distinct gene loci and then assembling each locus into as many isoforms as are needed to "
        "explain the data. Following this, StringTie simultaneously assembles and quantify the final transcripts by "
        "using network flow algorithm and starting from most highly abundant transcripts. The reference genome "
        "GFF annotation files, containing exon structures of 'known' genes, are then used to annotate the "
        "assembled transcripts and quantify the expression of known genes as well derive clues if a novel "
        "transcript has been found in the sample."
    )
    snippets["assembly"]["stringtie_methodology"] = stringtie_text

    # Hardcoded StringTie merge explanation
    stringtie_merge_text = (
        "After assembling each sample, the full set of assemblies is passed to StringTie's merge function, which "
        "merges together all the gene structures found in any of the samples. This step is required because "
        "transcripts in some of the samples might only be partially covered by reads, and as a consequence only "
        "partial versions of them will be assembled in the initial StringTie run. The merge step creates a set of "
        "transcripts that is consistent across all samples, so that the transcripts can be compared in subsequent "
        "steps in figure below:"
    )
    snippets["assembly"]["stringtie_merge"] = stringtie_merge_text

    stringtie_merge_caption = (
        "Explanation of merging transcript assemblies using StringTie’s merge function. In this "
        "example, four partial assemblies from four different samples are merged into two transcripts A and "
        "B. Samples 1 and 2 are both consistent with the reference annotation, which is used here to merge "
        "and extend them to create transcript A. Samples 3 and 4 are consistent with each other but not with "
        "the annotation, and these are merged to create transcript B. Please note that the 4 samples shown "
        "here are for example purpose only (Ref: Micheal P et al., 2016)."
    )
    snippets["assembly"]["stringtie_merge_caption"] = stringtie_merge_caption

    # Hardcoded gffcompare methodology
    gffcompare_text = (
        "For comparison of assembled transcripts with reference transcripts, Gffcompare utility was run taking "
        "the reference GTF and the string-tie merged GTF file. This produces an output file, which adds to each "
        "transcript a 'class code' and the name of the transcript from the reference annotation file to check how "
        "the predicted transcripts relate to an annotation file. Meaning of each class code is depicted in figure below:"
    )
    snippets["assembly"]["gffcompare_methodology"] = gffcompare_text

    # Hardcoded novel isoforms text
    novel_isoforms_text = (
        "To identify novel isoform transcripts not present in the reference GTF file, gffcompare utility was ran "
        "taking the reference GTF and the string-tie merged GTF file. This produces an output file, which adds to "
        "each transcript a 'class code' and the name of the transcript from the reference annotation file to check "
        "how the predicted transcripts relate to an annotation file. Class code 'j' means that predicted transcript "
        "is a potential novel isoform that shares at least one splice junction with a reference transcript. A total of "
        "8,934 novel isoforms with the class code 'j' were extracted and are provided with the deliverables in a "
        "file named novel.isoforms.gtf in the folder '03_transcript_assembly_gtf'. A brief description of the "
        "column names in the GTF file is provided in the same folder in the file named ReadMe.txt."
    )
    snippets["assembly"]["novel_isoforms"] = novel_isoforms_text

    # Differential expression (static narrative + tables)
    dge_intro_text = (
        "Abundances of the transcripts in all samples were estimated using StringTie with the help of merged "
        "transcripts generated from merge stringtie step. A python program (prepDE.py) was used to extract the "
        "read count information directly from the files generated by StringTie. For differential expression analysis, "
        "sample comparison was made according to grouping information provided by client:"
    )
    snippets["dge"]["intro"] = dge_intro_text

    snippets["dge"]["comparison_info_table"] = {
        "headers": ["Comparison", "Description"],
        "rows": [
            ["ST v/s SC", "Comparison 1 [ST (ST1+ST2) Vs SC (SC1+SC2)]"],
            ["RC v/s SC", "Comparison 2 [RC (RC1+RC2) Vs SC (SC1+SC2)]"],
            ["RT v/s SC", "Comparison 3 [RT (RT1+RT2) Vs SC (SC1+SC2)]"],
            ["RT v/s ST", "Comparison 4 [RT (RT1+RT2) Vs ST (ST1+ST2)]"],
            ["RC v/s ST", "Comparison 5 [RC (RC1+RC2) Vs ST (ST1+ST2)]"],
            ["RT v/s RC", "Comparison 6 [RT (RT1+RT2) Vs RC (RC1+RC2)]"],
        ],
    }

    snippets["dge"]["group_sample_table"] = {
        "headers": ["Group Name", "Sample Name"],
        "rows": [
            ["SC", "SC1"],
            ["SC", "SC2"],
            ["ST", "ST1"],
            ["ST", "ST2"],
            ["RC", "RC1"],
            ["RC", "RC2"],
            ["RT", "RT1"],
            ["RT", "RT2"],
        ],
    }

    deseq2_text = (
        "Differential gene expression was inferred between samples by applying the R package DESeq2. It is "
        "bioconductor package based on negative binomial distribution method. The analysis provides tabular "
        "result along with the normalized i.e. \"baseMean\" for each involved sample. The description of selected "
        "result columns has been given in below mentioned table:"
    )
    snippets["dge"]["deseq2_intro"] = deseq2_text

    snippets["dge"]["deseq2_columns_table"] = {
        "headers": ["Column Name", "Description"],
        "rows": [
            ["", "Feature identifier"],
            ["baseMean", "The average of the normalized count values, dividing by size factors, taken over all samples"],
            ["log2FoldChange", "The effect size estimate. This value indicates how much the gene or transcript's expression seems to have changed between the comparison and control groups. This value is reported on a logarithmic scale to base 2."],
            ["lfcSE", "The standard error estimate for the log2 fold change estimate."],
            ["stat", "The value of the test statistic for the gene or transcript."],
            ["pvalue", "P-value of the test for the gene or transcript."],
            ["padj", "Adjusted P-value for multiple testing for the gene or transcript."],
        ],
    }

    criteria_intro = (
        "The criterion used to identify upregulated and downregulated transcripts along with the significance is "
        "provided in the below table:"
    )
    snippets["dge"]["criteria_intro"] = criteria_intro

    snippets["dge"]["criteria_table"] = {
        "headers": ["Condition", "Status"],
        "rows": [
            ["log2FoldChange > 0", "Up regulated"],
            ["log2FoldChange < 0", "Down regulated"],
            ["log2FoldChange > 0 and padj <0.05", "Significantly up regulated"],
            ["log2FoldChange < 0 and padj <0.05", "Significantly down regulated"],
        ],
    }

    dge_stats_note = (
        "Differential expression analysis statistics is provided in table below and the data deliverable are given in "
        "the folder \"05_differential_expression_analysis\" which contains the Differential gene expression (DGE) "
        "analysis results."
    )
    snippets["dge"]["stats_note"] = dge_stats_note

    snippets["dge"]["annotation_note"] = (
        "Annotation for differentially expressed transcripts was obatined based on blastX result by similarity search "
        "against NCBI’s NR database in order to obtain gene ontology analysis using blast2GO."
    )

    snippets["dge"]["ma_plot_text"] = (
        "The plot visualizes the differences between measurements taken in two samples, by transforming the data "
        "onto M (log ratio) and A (mean average) scales, then plotting these values. For illustration purpose, MA "
        "plot for Comparison1 has been depicted below:"
    )
    snippets["dge"]["ma_plot_caption"] = (
        "Figure 10.1: MA plot showing differentially expressed transcripts in Comparison-1. On X-axis normalized "
        "counts for all the samples and on Y-axis log2foldchange are plotted. Points colored are with red if the "
        "q-value is less than 0.05 and black if the q-value is greater than 0.05."
    )
    snippets["dge"]["ma_plot_note"] = (
        "Note: MA plot for only one comparison has been shown above, however plots for all other combinations have "
        "been provided as a part of data deliverables \"05_differential_expression_analysis\"."
    )

    snippets["dge"]["volcano_plot_text"] = (
        "The \"volcano plot\" arranges expressed genes along dimensions of biological as well as statistical "
        "significance. For illustration purpose, Volcano plot for Comparison1 has been depicted below:"
    )
    snippets["dge"]["volcano_plot_caption"] = (
        "Figure 10.2: Volcano plot showing differentially expressed transcripts in Comparison-1. X-axis represents "
        "the log2fold change values and Y-axis represents log10 of q-value (FDR corrected p-value. Points are colored "
        "red if the q-value is less than 0.05 and black if the q-value is greater than 0.05."
    )
    snippets["dge"]["volcano_plot_note"] = (
        "Note: Volcano plot for only one comparison has been shown above, however plots for all other combinations "
        "has been provided as a part of data deliverables \"05_differential_expression_analysis\"."
    )

    # Gene Ontology (GO) analysis for significant DEG (static narrative + table)
    go_text = (
        "The Gene Ontology project provides controlled vocabularies of defined terms representing gene product "
        "properties. These cover three domains: Cellular Component, the parts of a cell or its extracellular "
        "environment; Molecular Function, the elemental activities of a gene product at the molecular level, such "
        "as binding or catalysis; and Biological Process, operations or sets of molecular events with a defined "
        "beginning and end, pertinent to the functioning of integrated living units: cells, tissues, organs, and "
        "organisms. GO was assigned to significant differentially expressed transcripts using Blast2go cli. Single "
        "gene can be assigned with multiple GO categories and hence multiple GO terms. The GO domain "
        "distribution is shown in table below:"
    )
    snippets["go"]["intro"] = go_text

    snippets["go"]["distribution_table"] = {
        "headers": [
            "Sample Name",
            "Significant DGE",
            "# Seq with GO",
            "Biological Process",
            "Cellular Component",
            "Molecular Function",
        ],
        "rows": [
            ["Comparison1", "1357", "498", "398", "351", "412"],
            ["Comparison2", "4616", "1492", "1178", "1022", "1189"],
            ["Comparison3", "5601", "1910", "1475", "1311", "1534"],
            ["Comparison4", "4916", "1551", "1194", "1036", "1214"],
            ["Comparison5", "5331", "1805", "1410", "1229", "1433"],
            ["Comparison6", "2793", "1215", "982", "879", "981"],
        ],
    }

    snippets["go"]["deliverables_note"] = (
        "Assigned Gene Ontology are also provided in deliverables \"06_Significant_DGE_GO\"."
    )

    return {
        "source": {
            "file": os.path.abspath(constants_ts_path),
            "generated_at": datetime.now(timezone.utc).isoformat(),
        },
        "snippets": snippets,
        "sections": sections,
    }


def write_static_content_json(output_path: str, constants_ts_path: str = _DEFAULT_CONSTANTS_TS) -> None:
    data = extract_static_content_from_constants(constants_ts_path)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(script_dir, "report_static_content.json")
    write_static_content_json(out)
    print(f"Wrote static content JSON: {out}")
