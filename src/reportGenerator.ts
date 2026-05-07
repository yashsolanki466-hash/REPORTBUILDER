import { launch, Browser } from 'puppeteer';
import ejs from 'ejs';
import fs from 'fs-extra';
import path from 'path';
import { Document, Paragraph, Packer, AlignmentType, HeadingLevel } from 'docx';
import chalk from 'chalk';
import ora from 'ora';
import { parseProjectData, ProjectData } from './dataParser.js';

interface GeneratorOptions {
  templateDir: string;
  headless: boolean;
}

interface GenerateOptions {
  inputDir: string;
  outputName: string;
  formats: string[];
  metadata?: {
    projectId?: string;
    piName?: string;
    clientName?: string;
    logo?: string;
  };
  template: string;
}

interface GenerateResult {
  files: Array<{
    path: string;
    size: number;
  }>;
  duration: number;
}

export class ReportGenerator {
  private options: GeneratorOptions;
  private browser: Browser | null = null;

  constructor(options: GeneratorOptions) {
    this.options = options;
  }

  async init(): Promise<void> {
    if (!this.browser) {
      this.browser = await launch({
        headless: this.options.headless,
        args: ['--no-sandbox', '--disable-setuid-sandbox']
      });
    }
  }

  async close(): Promise<void> {
    if (this.browser) {
      await this.browser.close();
      this.browser = null;
    }
  }

  async generate(options: GenerateOptions): Promise<GenerateResult> {
    const startTime = Date.now();
    const files: Array<{ path: string; size: number }> = [];

    await this.init();

    const spinner = ora('Analyzing project structure...').start();

    try {
      // Parse project data
      const projectData = await this.parseProjectData(options.inputDir, options.metadata);
      spinner.text = 'Loading template...';

      // Load template
      const templatePath = path.join(this.options.templateDir, `${options.template}.ejs`);
      const templateContent = await fs.readFile(templatePath, 'utf-8');

      // Render HTML
      spinner.text = 'Rendering HTML template...';
      const html = ejs.render(templateContent, projectData);

      // Write HTML file
      const htmlPath = path.join(options.inputDir, `${options.outputName}.html`);
      if (options.formats.includes('html')) {
        await fs.writeFile(htmlPath, html, 'utf-8');
        const stats = await fs.stat(htmlPath);
        files.push({ path: htmlPath, size: stats.size });
      }

      // Generate PDF using Puppeteer
      if (options.formats.includes('pdf')) {
        spinner.text = 'Generating PDF (using Chrome engine)...';
        const pdfPath = path.join(options.inputDir, `${options.outputName}.pdf`);
        await this.generatePDF(html, pdfPath, options.inputDir);
        const stats = await fs.stat(pdfPath);
        files.push({ path: pdfPath, size: stats.size });
      }

      // Generate DOCX
      if (options.formats.includes('docx')) {
        spinner.text = 'Generating DOCX...';
        const docxPath = path.join(options.inputDir, `${options.outputName}.docx`);
        await this.generateDOCX(projectData, docxPath);
        const stats = await fs.stat(docxPath);
        files.push({ path: docxPath, size: stats.size });
      }

      spinner.succeed('Report generation complete!');

      return {
        files,
        duration: Date.now() - startTime
      };
    } catch (error) {
      spinner.fail('Report generation failed');
      throw error;
    }
  }

  private async generatePDF(html: string, outputPath: string, basePath: string): Promise<void> {
    if (!this.browser) throw new Error('Browser not initialized');

    const page = await this.browser.newPage();

    try {
      // Write HTML to a temporary file
      const tempHtmlPath = path.join(basePath, 'temp.html');
      await fs.writeFile(tempHtmlPath, html, 'utf-8');

      // Load the HTML from file
      await page.goto(`file://${tempHtmlPath}`, {
        waitUntil: ['networkidle0', 'load', 'domcontentloaded'],
        timeout: 60000
      });

      // Wait for fonts to load
      await page.waitForFunction('document.fonts.ready', { timeout: 30000 }).catch(() => {});

      // Wait for Chart.js to render and convert canvases to images
      await page.addScriptTag({
        content: `
          setTimeout(() => {
            const canvases = document.querySelectorAll('canvas');
            canvases.forEach((canvas) => {
              try {
                const dataUrl = canvas.toDataURL('image/png');
                const img = document.createElement('img');
                img.src = dataUrl;
                img.style.width = canvas.style.width || canvas.width + 'px';
                img.style.height = canvas.style.height || canvas.height + 'px';
                img.style.maxWidth = '100%';
                img.style.height = 'auto';
                canvas.parentNode?.replaceChild(img, canvas);
              } catch (e) {
                console.warn('Could not convert canvas to image:', e);
              }
            });
          }, 2000);
        `
      });

      // Wait for canvas conversion to complete
      await new Promise(resolve => setTimeout(resolve, 4000));

      // Generate PDF with professional settings
      await page.pdf({
        path: outputPath,
        format: 'A4',
        printBackground: true,
        preferCSSPageSize: true,
        margin: {
          top: '20mm',
          right: '20mm',
          bottom: '20mm',
          left: '20mm'
        },
        displayHeaderFooter: true,
        headerTemplate: '<div></div>',
        footerTemplate: `
          <div style="font-size: 9px; text-align: center; width: 100%; color: #666; font-family: Georgia, serif;">
            Page <span class="pageNumber"></span> of <span class="totalPages"></span>
          </div>
        `
      });
    } finally {
      await page.close();
      // Clean up temp file
      try {
        await fs.unlink(path.join(basePath, 'temp.html'));
      } catch (e) {
        // Ignore cleanup errors
      }
    }
  }

  private async generateDOCX(data: any, outputPath: string): Promise<void> {
    const doc = new Document({
      sections: [{
        properties: {
          page: {
            margin: {
              top: 1440, // 1 inch = 1440 twips
              right: 1080,
              bottom: 1440,
              left: 1080
            }
          }
        },
        children: [
          new Paragraph({
            text: `NGS Analysis Report - ${data.project_id}`,
            heading: HeadingLevel.TITLE,
            alignment: AlignmentType.CENTER,
            spacing: { after: 400 }
          }),
          new Paragraph({
            text: `Project ID: ${data.project_id}`,
            spacing: { after: 200 }
          }),
          new Paragraph({
            text: `PI: ${data.pi_name}`,
            spacing: { after: 200 }
          }),
          new Paragraph({
            text: `Date: ${new Date().toLocaleDateString()}`,
            spacing: { after: 400 }
          }),
          // Add more sections based on data
        ]
      }]
    });

    const buffer = await Packer.toBuffer(doc);
    await fs.writeFile(outputPath, buffer);
  }

  private async parseProjectData(inputDir: string, metadata?: any): Promise<ProjectData> {
    return await parseProjectData(inputDir, metadata);
  }
}
