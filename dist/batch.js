import path from 'path';
import chalk from 'chalk';
import { ReportGenerator } from './reportGenerator.js';
import { findProjects } from './utils.js';
import cliProgress from 'cli-progress';
export async function batchGenerate(options) {
    console.log(chalk.blue(`Scanning ${options.rootDir} for projects...\n`));
    const projects = await findProjects(options.rootDir, options.maxDepth, options.pattern);
    if (projects.length === 0) {
        console.log(chalk.yellow('No NGS projects found in the specified directory.'));
        return;
    }
    console.log(chalk.green(`Found ${projects.length} project(s):`));
    projects.forEach((p, i) => {
        console.log(chalk.gray(`  ${i + 1}. ${path.basename(p)}`));
    });
    console.log();
    // Initialize generator
    const generator = new ReportGenerator({
        templateDir: path.join(process.cwd(), 'templates'),
        headless: true
    });
    // Create progress bar
    const progressBar = new cliProgress.SingleBar({
        format: '{bar} {percentage}% | {value}/{total} | {project}',
        barCompleteChar: '\u2588',
        barIncompleteChar: '\u2591',
        hideCursor: true
    });
    progressBar.start(projects.length, 0, { project: 'Initializing...' });
    let successCount = 0;
    let failCount = 0;
    for (let i = 0; i < projects.length; i++) {
        const project = projects[i];
        const projectName = path.basename(project);
        progressBar.update(i, { project: projectName });
        try {
            await generator.generate({
                inputDir: project,
                outputName: `${projectName}_report`,
                formats: options.formats,
                template: 'report_template'
            });
            successCount++;
        }
        catch (error) {
            console.error(chalk.red(`\nError processing ${projectName}:`), error);
            failCount++;
        }
        progressBar.update(i + 1);
    }
    progressBar.stop();
    await generator.close();
    console.log(chalk.green.bold(`\n✅ Batch complete: ${successCount} succeeded, ${failCount} failed\n`));
}
//# sourceMappingURL=batch.js.map