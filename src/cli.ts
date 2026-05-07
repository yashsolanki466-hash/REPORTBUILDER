#!/usr/bin/env node

import { Command } from 'commander';
import chalk from 'chalk';
import path from 'path';
import { fileURLToPath } from 'url';
import { ReportGenerator } from './reportGenerator.js';
import { validateProjectDirectory } from './utils.js';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const program = new Command();

program
  .name('ngs-report')
  .description('NGS Report Generator CLI - Generate professional HTML/PDF/DOCX reports')
  .version('1.0.0');

program
  .command('generate')
  .alias('g')
  .description('Generate report from project directory')
  .requiredOption('-i, --input <path>', 'Input project directory path')
  .option('-o, --output <filename>', 'Output filename (without extension)', 'report')
  .option('-f, --formats <formats>', 'Output formats: html,pdf,docx (comma-separated)', 'html,pdf,docx')
  .option('--project-id <id>', 'Project ID override')
  .option('--pi-name <name>', 'Principal Investigator name override')
  .option('--client-name <name>', 'Client name override')
  .option('--logo <path>', 'Custom logo path')
  .option('--template <name>', 'Template to use', 'report_full')
  .option('--headless <boolean>', 'Run puppeteer in headless mode', 'true')
  .action(async (options) => {
    console.log(chalk.blue.bold('\n🧬 NGS Report Generator\n'));
    
    try {
      // Validate input directory
      const inputDir = path.resolve(options.input);
      const validation = await validateProjectDirectory(inputDir);
      
      if (!validation.valid) {
        console.error(chalk.red('❌ Error:'), validation.error);
        process.exit(1);
      }
      
      console.log(chalk.gray(`Input: ${inputDir}`));
      console.log(chalk.gray(`Output: ${options.output}.{html,pdf,docx}`));
      console.log(chalk.gray(`Formats: ${options.formats}\n`));
      
      // Parse formats
      const formats = options.formats.split(',').map((f: string) => f.trim().toLowerCase());
      
      // Initialize generator
      const generator = new ReportGenerator({
        templateDir: path.join(__dirname, '..', 'templates'),
        headless: options.headless === 'true'
      });
      
      // Generate report
      const result = await generator.generate({
        inputDir,
        outputName: options.output,
        formats,
        metadata: {
          projectId: options.projectId,
          piName: options.piName,
          clientName: options.clientName,
          logo: options.logo
        },
        template: options.template
      });
      
      console.log(chalk.green.bold('\n✅ Report generation complete!\n'));
      
      // Display results
      result.files.forEach(file => {
        const size = formatFileSize(file.size);
        console.log(chalk.green(`  ✓ ${path.basename(file.path)} ${chalk.gray(`(${size})`)}`));
      });
      
      console.log(chalk.gray(`\nTotal time: ${result.duration}ms`));
      
    } catch (error) {
      console.error(chalk.red('❌ Error generating report:'), error);
      process.exit(1);
    }
  });

program
  .command('batch')
  .alias('b')
  .description('Generate reports for multiple projects in a directory')
  .requiredOption('-r, --root <path>', 'Root directory containing project folders')
  .option('-f, --formats <formats>', 'Output formats', 'html,pdf,docx')
  .option('--max-depth <n>', 'Maximum directory depth to scan', '2')
  .option('--pattern <pattern>', 'Pattern to identify project folders', '*')
  .action(async (options) => {
    console.log(chalk.blue.bold('\n🧬 NGS Batch Report Generator\n'));
    
    try {
      const { batchGenerate } = await import('./batch.js');
      await batchGenerate({
        rootDir: path.resolve(options.root),
        formats: options.formats.split(',').map((f: string) => f.trim()),
        maxDepth: parseInt(options.maxDepth),
        pattern: options.pattern
      });
    } catch (error) {
      console.error(chalk.red('❌ Error in batch processing:'), error);
      process.exit(1);
    }
  });

program
  .command('watch')
  .alias('w')
  .description('Watch directory and auto-generate reports when new projects are added')
  .requiredOption('-r, --root <path>', 'Root directory to watch')
  .option('-f, --formats <formats>', 'Output formats', 'html,pdf,docx')
  .option('--interval <ms>', 'Polling interval in milliseconds', '5000')
  .action(async (options) => {
    console.log(chalk.blue.bold('\n👁️  NGS Report Watcher\n'));
    console.log(chalk.gray(`Watching: ${path.resolve(options.root)}`));
    console.log(chalk.gray(`Interval: ${options.interval}ms\n`));
    
    try {
      const { startWatcher } = await import('./watch.js');
      await startWatcher({
        rootDir: path.resolve(options.root),
        formats: options.formats.split(',').map((f: string) => f.trim()),
        interval: parseInt(options.interval)
      });
    } catch (error) {
      console.error(chalk.red('❌ Error starting watcher:'), error);
      process.exit(1);
    }
  });

function formatFileSize(bytes: number): string {
  if (bytes === 0) return '0 B';
  const k = 1024;
  const sizes = ['B', 'KB', 'MB', 'GB'];
  const i = Math.floor(Math.log(bytes) / Math.log(k));
  return parseFloat((bytes / Math.pow(k, i)).toFixed(1)) + ' ' + sizes[i];
}

program.parse();
