import fs from 'fs-extra';
import path from 'path';
import { fileURLToPath } from 'url';
import XLSX from 'xlsx';
import chalk from 'chalk';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

export interface ProjectData {
  project_id: string;
  report_date: string;
  client_name: string;
  client_org: string;
  project_pi: string;
  application: string;
  no_of_samples: string;
  sample_count: number;
  samples: string[];
  qubit_data: any[];
  library_sizes: number[];
  sequencing_stats: any[];
  reference_organism: string;
  total_genes: number;
  ref_stats: any;
  mapping_stats: any[];
  total_transcripts: number;
  mean_transcript_size: number;
  assembly_stats: any[];
  diff_expr_stats: any[];
  dge_chart_labels: string[];
  dge_chart_up: number[];
  dge_chart_down: number[];
  pathway_stats: any[];
  go_distribution: any[];
  pathway_image_src: string;
  dge_figures: any[];
  dge_comparison_table: any;
  dge_group_table: any;
  func_assets: any;
  deliverables_tree: string;
  gffcompare_codes_src: string;
  workflow_figure_src: string;
  stringtie_merge_figure_src: string;
  isoforms_figure_src: string;
  pathway_ex_figure_src: string;
  logo_path: string;
  warnings: string[];
  static_content: any;
  static_snippets: any;
}

export async function parseProjectData(inputDir: string, metadataOverride?: any): Promise<ProjectData> {
  console.log(chalk.gray('Parsing project structure...'));
  
  const scriptDir = path.join(__dirname, '..');
  
  // Load static content
  const static_content = await loadStaticContent(scriptDir);
  const static_snippets = static_content?.snippets || {};
  
  // Parse README for metadata
  const readmeData = await parseReadme(inputDir);
  
  // Get project details
  const project_details = await getProjectDetails(inputDir);
  
  // Get warnings
  const warnings = await validateDeliverablesStructure(inputDir);
  
  // Get samples
  const samples = await getSamples(inputDir);
  const sample_count = samples.length;
  
  // Get sequencing stats
  const sequencing_stats = await getSequencingStats(inputDir, samples);
  
  // Get reference stats
  const ref_stats = await getReferenceStats(inputDir);
  const total_genes = ref_stats?.total_genes || 0;
  
  // Get mapping stats
  const mapping_stats = await getMappingStats(inputDir, samples);
  
  // Get assembly stats
  const assembly_stats = await getAssemblyStats(inputDir, samples);
  const total_transcripts = assembly_stats.find((s: any) => s.sample === 'merged.fasta')?.num_transcripts || 0;
  const mean_transcript_size = assembly_stats.find((s: any) => s.sample === 'merged.fasta')?.mean_size || 0;
  
  // Get DGE stats
  const { stats: diff_expr_stats, labels: dge_chart_labels, up_counts: dge_chart_up, down_counts: dge_chart_down } = await getDGEStats(inputDir);
  
  // Get pathway stats
  const pathway_stats = await getPathwayStats(inputDir);
  
  // Get GO distribution
  const go_distribution = await getGODistribution(inputDir);
  
  // Get DGE figures
  const dge_figures = await getDGEFigures(inputDir);
  
  // Get DGE comparison/group tables
  const { dge_comparison_table, dge_group_table } = await getDGEComparisonTables(inputDir);
  
  // Get functional assets
  const func_assets = await extractFunctionalAssets(inputDir);
  
  // Get deliverables tree
  const deliverables_tree = await getDeliverablesTree(inputDir);
  
  // Get component assets
  const gffcompare_codes_src = await findComponentAsset(scriptDir, 'gffcompare_codes.png', inputDir);
  const workflow_figure_src = await findComponentAsset(scriptDir, ['workflow.png', 'workflow.svg'], inputDir);
  const stringtie_merge_figure_src = await findComponentAsset(scriptDir, 'stringtie_merge_illustration.svg', inputDir);
  const isoforms_figure_src = await findComponentAsset(scriptDir, 'isoforms.png', inputDir);
  const pathway_ex_figure_src = await findComponentAsset(scriptDir, 'pathway_ex.png', inputDir);
  
  // Logo path
  const logo_path = metadataOverride?.logo || path.join(scriptDir, 'assets', 'logo.png');
  
  const reference_organism = metadataOverride?.reference_organism || ref_stats?.organism || 'Organism Name';
  
  return {
    project_id: project_details.project_id || metadataOverride?.projectId || readmeData.project_id || path.basename(inputDir),
    report_date: project_details.report_date || new Date().toLocaleDateString(),
    client_name: metadataOverride?.client_name || project_details.client_name || 'Unknown',
    client_org: metadataOverride?.client_org || project_details.client_org || 'Unknown',
    project_pi: metadataOverride?.piName || project_details.project_pi || readmeData.project_pi || '',
    application: project_details.application || readmeData.application || '',
    no_of_samples: project_details.no_of_samples || readmeData.no_of_samples || '',
    sample_count,
    samples,
    qubit_data: [],
    library_sizes: Array(sample_count).fill(450),
    sequencing_stats,
    reference_organism,
    total_genes,
    ref_stats,
    mapping_stats,
    total_transcripts,
    mean_transcript_size,
    assembly_stats,
    diff_expr_stats,
    dge_chart_labels,
    dge_chart_up,
    dge_chart_down,
    pathway_stats,
    go_distribution,
    pathway_image_src: '',
    dge_figures,
    dge_comparison_table,
    dge_group_table,
    func_assets,
    deliverables_tree,
    gffcompare_codes_src,
    workflow_figure_src,
    stringtie_merge_figure_src,
    isoforms_figure_src,
    pathway_ex_figure_src,
    logo_path,
    warnings,
    static_content,
    static_snippets
  };
}

async function loadStaticContent(scriptDir: string): Promise<any> {
  const staticContentPath = path.join(scriptDir, 'report_static_content.json');
  if (await fs.pathExists(staticContentPath)) {
    try {
      const content = await fs.readFile(staticContentPath, 'utf-8');
      return JSON.parse(content);
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not load static content JSON'));
    }
  }
  return {};
}

async function parseReadme(inputDir: string): Promise<any> {
  const readmePaths = ['Readme.txt', 'README.txt', 'readme.txt'];
  
  for (const readmeFile of readmePaths) {
    const readmePath = path.join(inputDir, readmeFile);
    if (await fs.pathExists(readmePath)) {
      try {
        const content = await fs.readFile(readmePath, 'utf-8');
        const details: any = {};
        
        const patterns = {
          project_id: /Project\s*ID\s*:\s*(.+)/i,
          project_pi: /Project\s*PI\s*:\s*(.+)/i,
          application: /Application\s*:\s*(.+)/i,
          no_of_samples: /No\s*of\s*Samples\s*:\s*(.+)/i,
        };
        
        for (const [key, pattern] of Object.entries(patterns)) {
          const match = content.match(pattern);
          if (match) {
            details[key] = match[1].trim();
          }
        }
        
        return details;
      } catch (e) {
        console.warn(chalk.yellow('Warning: Could not parse README'));
      }
    }
  }
  
  return {};
}

async function getProjectDetails(inputDir: string): Promise<any> {
  const metadataPath = path.join(inputDir, 'metadata.json');
  let metadata: any = {};
  
  if (await fs.pathExists(metadataPath)) {
    try {
      metadata = JSON.parse(await fs.readFile(metadataPath, 'utf-8'));
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not parse metadata.json'));
    }
  }
  
  return {
    project_id: metadata.project_id || '',
    project_pi: metadata.project_pi || '',
    client_name: metadata.client_name || '',
    client_org: metadata.client_org || '',
    application: metadata.application || '',
    no_of_samples: metadata.no_of_samples || '',
    report_date: new Date().toLocaleDateString()
  };
}

async function validateDeliverablesStructure(inputDir: string): Promise<string[]> {
  const warnings: string[] = [];
  const requiredDirs = [
    '01_Raw_Data',
    '02_reference_genome_and_gff',
    '03_Mapping',
    '04_transcript_assembly_gtf',
    '05_differential_expression_analysis',
    '06_Significant_DGE_GO',
    '07_Significant_DGE_pathways'
  ];
  
  for (const dir of requiredDirs) {
    const dirPath = path.join(inputDir, dir);
    if (!await fs.pathExists(dirPath)) {
      warnings.push(`Missing required directory: ${dir}`);
    }
  }
  
  return warnings;
}

async function getSamples(inputDir: string): Promise<string[]> {
  const samples: string[] = [];
  const rawDataDir = path.join(inputDir, '01_Raw_Data');
  
  if (await fs.pathExists(rawDataDir)) {
    try {
      const files = await fs.readdir(rawDataDir);
      for (const file of files) {
        if ((file.includes('Raw_stats') || file.includes('Stats')) && file.endsWith('.xlsx')) {
          const filePath = path.join(rawDataDir, file);
          const workbook = XLSX.readFile(filePath);
          const sheetName = workbook.SheetNames[0];
          const worksheet = workbook.Sheets[sheetName];
          const rawData = XLSX.utils.sheet_to_json(worksheet, { header: 1 }) as any[][];
          
          if (rawData.length > 1) {
            const headers = rawData[0] as string[];
            const fileIdx = headers.indexOf('file');
            
            const sampleNames = new Set<string>();
            for (let i = 1; i < rawData.length; i++) {
              const row = rawData[i];
              if (row[fileIdx]) {
                const filename = String(row[fileIdx]);
                const match = filename.match(/^([^_]+)/);
                if (match) {
                  sampleNames.add(match[1]);
                }
              }
            }
            
            samples.push(...Array.from(sampleNames));
          }
        }
      }
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not parse samples from Raw_Data'));
    }
  }
  
  return samples.length > 0 ? samples : [];
}

async function getSequencingStats(inputDir: string, samples: string[]): Promise<any[]> {
  const stats: any[] = [];
  const rawDataDir = path.join(inputDir, '01_Raw_Data');
  
  if (await fs.pathExists(rawDataDir)) {
    try {
      const files = await fs.readdir(rawDataDir);
      for (const file of files) {
        if ((file.includes('Raw_stats') || file.includes('Stats')) && file.endsWith('.xlsx')) {
          const filePath = path.join(rawDataDir, file);
          const workbook = XLSX.readFile(filePath);
          const sheetName = workbook.SheetNames[0];
          const worksheet = workbook.Sheets[sheetName];
          const data = XLSX.utils.sheet_to_json(worksheet);
          
          for (const row of data) {
            const r = row as any;
            stats.push({
              sample: r.file || r.Sample || 'Unknown',
              raw_reads: r['Total Reads'] || r.num_seqs || 'N/A',
              clean_reads: r['Total Reads'] || r.num_seqs || 'N/A',
              q30: r.Q30 || 'N/A',
              gc_content: r['GC%'] || 'N/A'
            });
          }
        }
      }
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not parse sequencing stats'));
    }
  }
  
  return stats;
}

async function getReferenceStats(inputDir: string): Promise<any> {
  const refDir = path.join(inputDir, '02_reference_genome_and_gff');
  
  if (await fs.pathExists(refDir)) {
    try {
      const files = await fs.readdir(refDir);
      for (const file of files) {
        if (file.endsWith('.fasta') || file.endsWith('.fa') || file.endsWith('.fna')) {
          const filePath = path.join(refDir, file);
          const content = await fs.readFile(filePath, 'utf-8');
          const lines = content.split('\n');
          const geneCount = lines.filter(line => line.startsWith('>')).length;
          
          return {
            organism: 'Organism Name',
            total_genes: geneCount,
            source: file
          };
        }
      }
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not parse reference stats'));
    }
  }
  
  return { organism: 'Organism Name', total_genes: 0 };
}

async function getMappingStats(inputDir: string, samples: string[]): Promise<any[]> {
  const stats: any[] = [];
  const mappingDir = path.join(inputDir, '03_Mapping');
  
  if (await fs.pathExists(mappingDir)) {
    try {
      const files = await fs.readdir(mappingDir);
      for (const file of files) {
        if (file.endsWith('.xlsx') || file.endsWith('.xls')) {
          const filePath = path.join(mappingDir, file);
          const workbook = XLSX.readFile(filePath);
          const sheetName = workbook.SheetNames[0];
          const worksheet = workbook.Sheets[sheetName];
          const data = XLSX.utils.sheet_to_json(worksheet);
          
          for (const row of data) {
            const r = row as any;
            stats.push({
              sample_name: r.Sample || r['Sample Name'] || 'Unknown',
              mapped_reads: r['Mapped Reads'] || 'N/A',
              percent_mapped: r['% Mapped'] || 'N/A',
              unique_reads: r['Uniquely Mapped'] || 'N/A',
              percent_unique: r['% Unique'] || 'N/A'
            });
          }
        }
      }
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not parse mapping stats'));
    }
  }
  
  return stats;
}

async function getAssemblyStats(inputDir: string, samples: string[]): Promise<any[]> {
  const stats: any[] = [];
  const assemblyDir = path.join(inputDir, '04_transcript_assembly_gtf');
  
  if (await fs.pathExists(assemblyDir)) {
    try {
      const files = await fs.readdir(assemblyDir);
      for (const file of files) {
        if (file.endsWith('.xlsx') || file.endsWith('.xls')) {
          const filePath = path.join(assemblyDir, file);
          const workbook = XLSX.readFile(filePath);
          const sheetName = workbook.SheetNames[0];
          const worksheet = workbook.Sheets[sheetName];
          const data = XLSX.utils.sheet_to_json(worksheet);
          
          for (const row of data) {
            const r = row as any;
            stats.push({
              sample: r.Sample || r['Sample Name'] || 'Unknown',
              num_transcripts: r['Number of Transcripts'] || r.Transcripts || '0',
              mean_size: r['Mean Size'] || r['Mean Transcript Size'] || '0'
            });
          }
        }
      }
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not parse assembly stats'));
    }
  }
  
  return stats;
}

async function getDGEStats(inputDir: string): Promise<{ stats: any[]; labels: string[]; up_counts: number[]; down_counts: number[] }> {
  const stats: any[] = [];
  const labels: string[] = [];
  const up_counts: number[] = [];
  const down_counts: number[] = [];
  const dgeDir = path.join(inputDir, '05_differential_expression_analysis');
  
  if (await fs.pathExists(dgeDir)) {
    try {
      const files = await fs.readdir(dgeDir);
      for (const file of files) {
        if (file.toLowerCase().includes('dge') && file.endsWith('.xlsx')) {
          const base = file.replace(/_?dge.*$/i, '').trim();
          const comparison = base || file;
          
          const filePath = path.join(dgeDir, file);
          const workbook = XLSX.readFile(filePath);
          const sheetName = workbook.SheetNames[0];
          const worksheet = workbook.Sheets[sheetName];
          const data = XLSX.utils.sheet_to_json(worksheet);
          
          let total = 0;
          let sig_up = 0;
          let sig_down = 0;
          
          for (const row of data) {
            const r = row as any;
            total++;
            const logFC = parseFloat(r.logFC || r['log2FoldChange'] || 0);
            const fdr = parseFloat(r.FDR || r.padj || 1);
            
            if (logFC > 0 && fdr < 0.05) sig_up++;
            if (logFC < 0 && fdr < 0.05) sig_down++;
          }
          
          stats.push({
            comparison,
            total_genes: total,
            sig_up,
            sig_down,
            total_sig: sig_up + sig_down
          });
          
          labels.push(comparison);
          up_counts.push(sig_up);
          down_counts.push(sig_down);
        }
      }
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not parse DGE stats'));
    }
  }
  
  return { stats, labels, up_counts, down_counts };
}

async function getPathwayStats(inputDir: string): Promise<any[]> {
  const stats: any[] = [];
  const pathwayDir = path.join(inputDir, '07_Significant_DGE_pathways');
  
  if (await fs.pathExists(pathwayDir)) {
    try {
      const files = await fs.readdir(pathwayDir);
      for (const file of files) {
        if (file.toLowerCase().includes('kegg') && file.endsWith('.xlsx')) {
          const filePath = path.join(pathwayDir, file);
          const workbook = XLSX.readFile(filePath);
          const sheetName = workbook.SheetNames[0];
          const worksheet = workbook.Sheets[sheetName];
          const data = XLSX.utils.sheet_to_json(worksheet);
          
          for (const row of data.slice(0, 20)) {
            const r = row as any;
            stats.push({
              pathway: r.Pathway || r['Pathway Name'] || 'Unknown',
              count: r.Count || r['Gene Count'] || 0
            });
          }
        }
      }
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not parse pathway stats'));
    }
  }
  
  return stats;
}

async function getGODistribution(inputDir: string): Promise<any[]> {
  const stats: any[] = [];
  const goDir = path.join(inputDir, '06_Significant_DGE_GO');
  
  if (await fs.pathExists(goDir)) {
    try {
      const files = await fs.readdir(goDir);
      for (const file of files) {
        if (file.toLowerCase().includes('go') && file.endsWith('.xlsx')) {
          const filePath = path.join(goDir, file);
          const workbook = XLSX.readFile(filePath);
          const sheetName = workbook.SheetNames[0];
          const worksheet = workbook.Sheets[sheetName];
          const data = XLSX.utils.sheet_to_json(worksheet);
          
          for (const row of data.slice(0, 30)) {
            const r = row as any;
            stats.push({
              term: r.Term || r['GO Term'] || 'Unknown',
              count: r.Count || r['Gene Count'] || 0,
              ontology: r.Namespace || r.Ontology || r.Category || 'Unknown'
            });
          }
        }
      }
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not parse GO distribution'));
    }
  }
  
  return stats;
}

async function getDGEFigures(inputDir: string): Promise<any[]> {
  const dgeDir = path.join(inputDir, '05_differential_expression_analysis');
  const figures: any[] = [];

  if (!await fs.pathExists(dgeDir)) {
    return figures;
  }

  try {
    const files = await fs.readdir(dgeDir);

    for (const file of files) {
      if (file.toLowerCase().includes('dge') && file.endsWith('.xlsx')) {
        const base = file.replace(/_?dge.*$/i, '').trim();
        const comparison = base || file;

        // Find heatmap
        let heatmap = '';
        for (const pattern of [
          `${comparison}*Heatmap*.png`,
          `${comparison}*heatmap*.png`,
          '*Heatmap*.png'
        ]) {
          const matches = files.filter(f => f.match(new RegExp(pattern.replace(/\*/g, '.*'), 'i')));
          if (matches.length > 0) {
            heatmap = path.join(dgeDir, matches[0]);
            break;
          }
        }

        // Find MA/Volcano PDF
        let maVolcano = '';
        for (const pattern of [
          `${comparison}*MA*Volcano*.pdf`,
          `${comparison}*Volcano*.pdf`,
          '*Volcano*.pdf'
        ]) {
          const matches = files.filter(f => f.match(new RegExp(pattern.replace(/\*/g, '.*'), 'i')));
          if (matches.length > 0) {
            maVolcano = path.join(dgeDir, matches[0]);
            break;
          }
        }

        // Compute plot data from DGE Excel file
        let plotData: { ma: { [key: string]: any[] }, volcano: { [key: string]: any[] } } = { ma: { up: [], down: [], ns: [] }, volcano: { up: [], down: [], ns: [] } };
        let hasMAData = false;
        let hasVolcanoData = false;

        try {
          const xlsxPath = path.join(dgeDir, file);
          console.log(chalk.blue(`Computing plot data from: ${xlsxPath}`));
          const workbook = XLSX.readFile(xlsxPath);
          const sheetName = workbook.SheetNames[0];
          const worksheet = workbook.Sheets[sheetName];
          const data = XLSX.utils.sheet_to_json(worksheet) as any[];
          console.log(chalk.blue(`Loaded ${data.length} rows from DGE file`));

          // Find columns with fallback logic
          let fcCol = findColumn(data, ['logfc', 'log2fc', 'log2foldchange', 'logfoldchange', 'lfc']);
          if (!fcCol) {
            fcCol = Object.keys(data[0]).find((c: string) => normCol(c).includes('logfc')) || null;
            if (!fcCol) {
              fcCol = Object.keys(data[0]).find((c: string) => normCol(c).includes('log2fold')) || null;
            }
          }

          let cpmCol = findColumn(data, ['logcpm', 'log_cpm']);
          if (!cpmCol) {
            cpmCol = Object.keys(data[0]).find((c: string) => normCol(c).includes('logcpm')) || null;
          }

          let meanCol = findColumn(data, ['basemean', 'aveexpr']);
          if (!meanCol) {
            meanCol = Object.keys(data[0]).find((c: string) => normCol(c).includes('basemean')) || null;
          }

          let sigCol = findColumn(data, ['padj', 'adjp', 'adjpval', 'fdr', 'qvalue', 'qval', 'q_value', 'pvalue', 'pval', 'p_value']);
          if (!sigCol) {
            sigCol = Object.keys(data[0]).find((c: string) => {
              const norm = normCol(c);
              return ['padj', 'fdr', 'qvalue', 'qval'].some(tok => norm.includes(tok));
            }) || null;
            if (!sigCol) {
              sigCol = Object.keys(data[0]).find((c: string) => {
                const norm = normCol(c);
                return ['pvalue', 'pval'].some(tok => norm.includes(tok));
              }) || null;
            }
          }

          console.log(chalk.blue(`Columns found - FC: ${fcCol}, CPM: ${cpmCol}, Mean: ${meanCol}, Sig: ${sigCol}`));

          if (fcCol) {
            for (const row of data) {
              const fc = asFloat(row[fcCol]);
              if (fc === null) continue;

              const cpm = cpmCol ? asFloat(row[cpmCol]) : null;
              const meanExpr = meanCol ? asFloat(row[meanCol]) : null;
              const pval = sigCol ? asFloat(row[sigCol]) : null;

              const isSig = (pval !== null && pval > 0 && pval < 0.05);
              let bucket = 'ns';
              if (isSig && fc >= 1) bucket = 'up';
              else if (isSig && fc <= -1) bucket = 'down';

              // MA plot data
              let maX = null;
              if (cpm !== null) {
                maX = cpm;
              } else if (meanExpr !== null && meanExpr >= 0) {
                maX = Math.log10(meanExpr + 1.0);
              }
              if (maX !== null) {
                plotData.ma[bucket].push({ x: maX, y: fc });
              }

              // Volcano plot data
              if (pval !== null && pval > 0) {
                plotData.volcano[bucket].push({ x: fc, y: -Math.log10(pval) });
              }
            }

            // Downsample data
            const maxPoints = 8000;
            for (const k in plotData.ma) {
              if (plotData.ma[k].length > maxPoints) {
                const step = Math.ceil(plotData.ma[k].length / maxPoints);
                plotData.ma[k] = plotData.ma[k].filter((_, i) => i % step === 0);
              }
            }
            for (const k in plotData.volcano) {
              if (plotData.volcano[k].length > maxPoints) {
                const step = Math.ceil(plotData.volcano[k].length / maxPoints);
                plotData.volcano[k] = plotData.volcano[k].filter((_, i) => i % step === 0);
              }
            }

            hasMAData = Object.values(plotData.ma).some(arr => arr.length > 0);
            hasVolcanoData = Object.values(plotData.volcano).some(arr => arr.length > 0);

            console.log(chalk.blue(`MA data points: ${Object.values(plotData.ma).flat().length}, Volcano data points: ${Object.values(plotData.volcano).flat().length}`));
          } else {
            console.warn(chalk.yellow(`Could not find fold change column in ${file}`));
          }
        } catch (e) {
          console.warn(chalk.yellow(`Warning: Could not compute plot data for ${file}: ${e}`));
        }

        figures.push({
          comparison,
          heatmap_png_src: heatmap,
          ma_volcano_pdf_src: maVolcano,
          dge_xlsx_src: path.join(dgeDir, file),
          plot_data_json: JSON.stringify(plotData),
          ma_plot_data: hasMAData,
          volcano_plot_data: hasVolcanoData
        });
      }
    }
  } catch (e) {
    console.warn(chalk.yellow('Warning: Could not parse DGE figures'));
  }

  return figures;
}

function normCol(c: string): string {
  return c.toLowerCase().replace(/[^a-z0-9]+/g, '');
}

function asFloat(v: any): number | null {
  if (v === null || v === undefined) return null;
  if (typeof v === 'number') return v;
  try {
    const s = String(v).trim().replace(/,/g, '');
    if (s === '' || s.toLowerCase() === 'na' || s.toLowerCase() === 'nan' || s.toLowerCase() === 'none') return null;
    const f = parseFloat(s);
    return isNaN(f) ? null : f;
  } catch {
    return null;
  }
}

function findColumn(data: any[], candidates: string[]): string | null {
  if (!data || data.length === 0) return null;
  const columns = Object.keys(data[0]);
  const normMap: { [key: string]: string } = {};
  for (const col of columns) {
    normMap[normCol(col)] = col;
  }

  for (const cand of candidates) {
    const nc = normCol(cand);
    if (nc in normMap) {
      return normMap[nc];
    }
  }

  return null;
}

async function getDGEComparisonTables(inputDir: string): Promise<{ dge_comparison_table: any; dge_group_table: any }> {
  // Return static tables from static content or defaults
  return {
    dge_comparison_table: {
      headers: ['Comparison', 'Description'],
      rows: []
    },
    dge_group_table: {
      headers: ['Group Name', 'Sample Name'],
      rows: []
    }
  };
}

async function extractFunctionalAssets(inputDir: string): Promise<any> {
  const assets: any = {
    go_results_xlsx: '',
    kegg_results_xlsx: '',
    go_results_preview: [],
    kegg_results_preview: [],
    kegg_pathway_image_src: '',
    enrichment_image_src: ''
  };

  // Recursively search for GO and KEGG files
  async function findFiles(pattern: string, extensions: string[]): Promise<string[]> {
    const found: string[] = [];
    async function search(dir: string) {
      if (await fs.pathExists(dir)) {
        const entries = await fs.readdir(dir, { withFileTypes: true });
        for (const entry of entries) {
          const fullPath = path.join(dir, entry.name);
          if (entry.isDirectory()) {
            await search(fullPath);
          } else if (entry.isFile()) {
            const nameLower = entry.name.toLowerCase();
            if (nameLower.includes(pattern) && extensions.some(ext => nameLower.endsWith(ext))) {
              found.push(fullPath);
            }
          }
        }
      }
    }
    await search(inputDir);
    return found;
  }

  // Find GO files
  const goFiles = await findFiles('go', ['.xlsx', '.xls', '.csv', '.tsv']);
  if (goFiles.length > 0) {
    assets.go_results_xlsx = goFiles[0]; // Absolute path
    // Read preview
    try {
      const workbook = XLSX.readFile(goFiles[0]);
      const sheetName = workbook.SheetNames[0];
      const worksheet = workbook.Sheets[sheetName];
      const data = XLSX.utils.sheet_to_json(worksheet, { header: 1 }) as any[][];
      assets.go_results_preview = data.slice(0, 11).map((row, idx) => {
        if (idx === 0) return row; // Keep header
        return row.map((cell: any) => cell ?? '');
      });
    } catch (e) {
      console.warn('Could not read GO file preview');
    }
  }

  // Find KEGG files
  const keggFiles = await findFiles('kegg', ['.xlsx', '.xls', '.csv', '.tsv']);
  if (keggFiles.length > 0) {
    assets.kegg_results_xlsx = keggFiles[0]; // Absolute path
    // Read preview
    try {
      const workbook = XLSX.readFile(keggFiles[0]);
      const sheetName = workbook.SheetNames[0];
      const worksheet = workbook.Sheets[sheetName];
      const data = XLSX.utils.sheet_to_json(worksheet, { header: 1 }) as any[][];
      assets.kegg_results_preview = data.slice(0, 11).map((row, idx) => {
        if (idx === 0) return row; // Keep header
        return row.map((cell: any) => cell ?? '');
      });
    } catch (e) {
      console.warn('Could not read KEGG file preview');
    }
  }

  // Search for images
  async function findImages(pattern: string): Promise<string[]> {
    const found: string[] = [];
    async function search(dir: string) {
      if (await fs.pathExists(dir)) {
        const entries = await fs.readdir(dir, { withFileTypes: true });
        for (const entry of entries) {
          const fullPath = path.join(dir, entry.name);
          if (entry.isDirectory()) {
            await search(fullPath);
          } else if (entry.isFile()) {
            const nameLower = entry.name.toLowerCase();
            if (nameLower.includes(pattern) && ['.png', '.jpg', '.jpeg', '.svg'].some(ext => nameLower.endsWith(ext))) {
              found.push(fullPath);
            }
          }
        }
      }
    }
    await search(inputDir);
    return found;
  }

  // Find KEGG pathway images
  const keggPathImages = await findImages('kegg');
  const pathwayImages = keggPathImages.filter(f => path.basename(f).toLowerCase().includes('path'));
  if (pathwayImages.length > 0) {
    assets.kegg_pathway_image_src = pathwayImages[0]; // Absolute path
  }

  // Find enrichment images
  const enrichImages = await findImages('enrich');
  const dotplotImages = await findImages('dotplot');
  const allEnrichImages = [...enrichImages, ...dotplotImages];
  if (allEnrichImages.length > 0) {
    assets.enrichment_image_src = allEnrichImages[0]; // Absolute path
  }

  
  return assets;
}

async function getDeliverablesTree(inputDir: string): Promise<string> {
  const readmePath = path.join(inputDir, 'Readme.txt');
  if (await fs.pathExists(readmePath)) {
    try {
      const content = await fs.readFile(readmePath, 'utf-8');
      const lines = content.split('\n');
      const treeLines = lines.filter(line => 
        line.includes('|--') || line.includes('├──') || line.includes('└──') || line.includes('│   ') || line.includes('--')
      );
      if (treeLines.length > 0) {
        return treeLines.join('\n').trim();
      }
    } catch (e) {
      console.warn(chalk.yellow('Warning: Could not extract tree from README'));
    }
  }
  
  // Build tree from filesystem
  return await buildDeliverablesTreeFromFS(inputDir);
}

async function buildDeliverablesTreeFromFS(inputDir: string, maxDepth: number = 4): Promise<string> {
  const rootDepth = inputDir.split(path.sep).length;
  const out: string[] = [path.basename(inputDir) + '/'];
  
  async function walk(dir: string, depth: number) {
    if (depth >= maxDepth) return;
    
    try {
      const entries = await fs.readdir(dir, { withFileTypes: true });
      const dirs = entries.filter(e => e.isDirectory()).sort((a, b) => a.name.localeCompare(b.name));
      const files = entries.filter(e => e.isFile() && !e.name.startsWith('.')).sort((a, b) => a.name.localeCompare(b.name));
      
      for (const d of dirs) {
        const fullPath = path.join(dir, d.name);
        const relPath = path.relative(inputDir, fullPath);
        const currentDepth = relPath.split(path.sep).length;
        const indent = '  '.repeat(currentDepth);
        out.push(`${indent}${d.name}/`);
        await walk(fullPath, currentDepth + 1);
      }
      
      for (let i = 0; i < Math.min(files.length, 25); i++) {
        const f = files[i];
        const relPath = path.relative(inputDir, path.join(dir, f.name));
        const currentDepth = relPath.split(path.sep).length;
        const indent = '  '.repeat(currentDepth);
        out.push(`${indent}${f.name}`);
      }
    } catch (e) {
      // Ignore errors
    }
  }
  
  await walk(inputDir, 0);
  return out.join('\n');
}

async function findComponentAsset(scriptDir: string, assetName: string | string[], inputDir: string): Promise<string> {
  const names = Array.isArray(assetName) ? assetName : [assetName];
  const componentsDir = path.join(scriptDir, 'components');
  
  for (const name of names) {
    const assetPath = path.join(componentsDir, name);
    if (await fs.pathExists(assetPath)) {
      // Return absolute path for the template to use
      return assetPath;
    }
  }
  
  return '';
}

async function getAllFilesRecursive(dir: string): Promise<string[]> {
  const files: string[] = [];
  
  try {
    const entries = await fs.readdir(dir, { withFileTypes: true });
    
    for (const entry of entries) {
      const fullPath = path.join(dir, entry.name);
      if (entry.isDirectory()) {
        const subFiles = await getAllFilesRecursive(fullPath);
        files.push(...subFiles);
      } else {
        files.push(fullPath);
      }
    }
  } catch (e) {
    // Ignore errors
  }
  
  return files;
}
