export const LATEX_TEMPLATE = `\\documentclass[a4paper,11pt]{article}

% --- Packages ---
\\usepackage[utf8]{inputenc}
\\usepackage[T1]{fontenc}
\\usepackage{lmodern}
\\usepackage[margin=2.5cm, top=4cm, bottom=3.5cm, headheight=80pt]{geometry}
\\usepackage{graphicx}
\\usepackage{booktabs}
\\usepackage[table]{xcolor}
\\usepackage{hyperref}
\\usepackage{fancyhdr}
\\usepackage{titlesec}
\\usepackage{longtable}
\\usepackage{array}
\\usepackage{float}
\\usepackage{caption}
\\usepackage{tabularx}
\\usepackage{enumitem}
\\usepackage{tikz}
\\usepackage{mdframed}

% --- Custom Colors ---
\\definecolor{uniorange}{HTML}{ED8936}
\\definecolor{uniblue}{HTML}{2563EB}
\\definecolor{unidarkblue}{HTML}{1E3A8A}
\\definecolor{unigray}{HTML}{718096}
\\definecolor{unilight}{HTML}{F7FAFC}

% --- Title & Section Styling ---
\\hypersetup{
    colorlinks=true,
    linkcolor=uniblue,
    filecolor=uniblue,      
    urlcolor=uniblue,
}

\\titleformat{\\section}{\\color{unidarkblue}\\normalfont\\Large\\bfseries}{\\thesection}{1em}{}
\\titleformat{\\subsection}{\\color{unidarkblue}\\normalfont\\large\\bfseries}{\\thesubsection}{1em}{}

% --- Custom Header & Footer ---
\\pagestyle{fancy}
\\fancyhf{}

\\renewcommand{\\headrulewidth}{0pt}
\\renewcommand{\\footrulewidth}{0.5pt}

% Header Implementation
\\fancyhead[C]{
    \\begin{tikzpicture}[remember picture, overlay]
        % Orange side bars similar to image
        \\fill[uniorange] ([xshift=0pt, yshift=0pt]current page.north west) rectangle ([xshift=15pt, yshift=-100pt]current page.north west);
        \\fill[uniorange] ([xshift=0pt, yshift=0pt]current page.north east) rectangle ([xshift=-15pt, yshift=-100pt]current page.north east);
    \\end{tikzpicture}
    \\vspace{-1cm}
    \\begin{flushright}
        \\includegraphics[width=4cm]{unigenome_logo} \\\\
        {\\tiny \\color{unigray} Lab Facility: 2A,3A,3B PASL -House, Beside Sahjanand college, Opposite Kamdenu Complex, Panjarapole, Ambawadi, Ahmedabad-380015, Gujarat \\\\
        Ph: +91-79-49006800 | WhatsApp: 6356005900 | Email: info@unigenome.in | Website: www.unigenome.in \\\\
        \\textbf{{{PROJECT_ID}}} | Date: {{DATE}} }
    \\end{flushright}
}

% Footer Implementation
\\fancyfoot[L]{\\color{unigray}\\small Unigenome | Confidential Restricted use only}
\\fancyfoot[R]{
    \\begin{tikzpicture}[baseline=(pagebox.base)]
        \\node[fill=uniorange, text=white, font=\\bfseries, inner sep=4pt, minimum width=1cm] (pagebox) {\\thepage};
    \\end{tikzpicture}
}

% --- Document Commands ---
\\newcommand{\\bluecaption}[1]{
    \\captionsetup{font={color=unidarkblue, it, small}}
    \\caption{#1}
}

% --- Document ---
\\begin{document}

\\begin{titlepage}
    \\centering
    \\vspace*{2cm}
    \\includegraphics[width=0.5\\textwidth]{unigenome_logo} \\\\
    \\vspace{2cm}
    {\\Huge\\bfseries\\color{unidarkblue} RNA-SEQUENCING ANALYSIS REPORT \\par}
    \\vspace{1cm}
    {\\Large\\color{uniorange} Comprehensive Transcriptome Profiling \\par}
    \\vspace{2cm}
    
    \\begin{table}[H]
        \\centering
        \\renewcommand{\\arraystretch}{1.5}
        \\rowcolors{2}{unilight}{white}
        \\begin{tabularx}{0.8\\textwidth}{lX}
            \\hline
            \\textbf{Project ID} & {{PROJECT_ID}} \\\\
            \\textbf{Client Name} & {{CLIENT}} \\\\
            \\textbf{Institute} & {{INSTITUTE}} \\\\
            \\textbf{Organism} & {{ORGANISM}} \\\\
            \\textbf{Genome Build} & {{GENOME_BUILD}} \\\\
            \\textbf{Platform} & {{PLATFORM}} \\\\
            \\textbf{Service Type} & {{SERVICE_TYPE}} \\\\
            \\textbf{Date} & {{DATE}} \\\\
            \\hline
        \\end{tabularx}
    \\end{table}

    \\vfill
    {\\large \\today \\par}
\\end{titlepage}

\\tableofcontents
\\newpage

\\section{Executive Summary}
This report provides an in-depth establishment of the RNA-Seq transcriptome analysis for project \\textbf{{{PROJECT_ID}}}.

\\subsection{Global Statistics}
\\begin{table}[H]
\\centering
\\caption{High-Level Project Statistics}
\\rowcolors{2}{unilight}{white}
\\begin{tabularx}{\\textwidth}{Xl}
\\toprule
\\textbf{Metric} & \\textbf{Value} \\\\
\\midrule
Total Samples & {{TOTAL_SAMPLES}} \\\\
Total Data & {{TOTAL_DATA_GB}} GB \\\\
Mean Mapping Rate & {{MAPPING_RATE}} \\\\
Total Transcripts & {{MERGED_TRANSCRIPTS}} \\\\
Novel Isoforms & {{NOVEL_ISOFORMS}} \\\\
\\bottomrule
\\end{tabularx}
\\end{table}

\\section{Methods \\& Workflow}
\\subsection{Sequence Processing}
\\begin{enumerate}
    \\item \\textbf{Quality Control}: Raw reads were assessed for quality using FastQC.
    \\item \\textbf{Trimming}: Adapters and low quality bases were removed.
    \\item \\textbf{Alignment}: Reads were mapped to the reference genome \\textbf{{{GENOME_BUILD}}} using STAR.
    \\item \\textbf{Quantification}: Transcript assembly and quantification were performed using StringTie.
\\end{enumerate}

\\newpage
{{TRANSCRIPT_STATS}}

\\section{Quality Control establishment}
\\subsection{Sequencing Yield Stats}
{{DATA_STATS_TABLE}}

\\subsection{Mapping Performance}
{{MAPPING_STATS_TABLE}}

\\newpage
\\section{Differential Expression establishment}
\\subsection{Summary of Results}
{{DGE_SUMMARY_TABLE}}

\\newpage
\\section{Functional Enrichment}
\\subsection{Gene Ontology (GO) Detail}
{{GO_DETAIL}}

\\newpage
\\subsection{KEGG Pathway Detail}
{{KEGG_DETAIL}}

\\newpage
\\section{Deliverables}
\\begin{verbatim}
{{DELIVERABLES_TREE}}
\\end{verbatim}

\\end{document}`;
