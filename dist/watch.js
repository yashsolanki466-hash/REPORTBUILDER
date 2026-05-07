import chokidar from 'chokidar';
import path from 'path';
import chalk from 'chalk';
import { ReportGenerator } from './reportGenerator.js';
import { validateProjectDirectory } from './utils.js';
const processedProjects = new Set();
export async function startWatcher(options) {
    const generator = new ReportGenerator({
        templateDir: path.join(process.cwd(), 'templates'),
        headless: true
    });
    console.log(chalk.blue('Starting watcher... Press Ctrl+C to stop\n'));
    // Initial scan
    await scanAndProcess(options.rootDir, generator, options.formats);
    // Watch for changes
    const watcher = chokidar.watch(options.rootDir, {
        ignored: /(^|[\/\\])\../, // Ignore dotfiles
        persistent: true,
        depth: 2,
        ignoreInitial: true,
        awaitWriteFinish: {
            stabilityThreshold: 2000,
            pollInterval: 100
        }
    });
    watcher
        .on('addDir', async (dirPath) => {
        console.log(chalk.gray(`📁 New directory detected: ${path.basename(dirPath)}`));
        // Check if it's a valid project
        const validation = await validateProjectDirectory(dirPath);
        if (validation.valid && !processedProjects.has(dirPath)) {
            await processProject(dirPath, generator, options.formats);
        }
    })
        .on('add', async (filePath) => {
        // If a Readme.txt is added, check if parent directory is a project
        if (path.basename(filePath).toLowerCase().includes('readme')) {
            const dirPath = path.dirname(filePath);
            if (!processedProjects.has(dirPath)) {
                const validation = await validateProjectDirectory(dirPath);
                if (validation.valid) {
                    console.log(chalk.gray(`📄 Readme file detected in: ${path.basename(dirPath)}`));
                    await processProject(dirPath, generator, options.formats);
                }
            }
        }
    })
        .on('error', (error) => {
        console.error(chalk.red('Watcher error:'), error);
    });
    // Keep process alive
    process.on('SIGINT', async () => {
        console.log(chalk.yellow('\n\nStopping watcher...'));
        await watcher.close();
        await generator.close();
        console.log(chalk.green('Watcher stopped.'));
        process.exit(0);
    });
    // Prevent process from exiting
    await new Promise(() => { });
}
async function scanAndProcess(rootDir, generator, formats) {
    console.log(chalk.blue('Scanning for existing projects...'));
    const { findProjects } = await import('./utils.js');
    const projects = await findProjects(rootDir, 2);
    for (const project of projects) {
        if (!processedProjects.has(project)) {
            await processProject(project, generator, formats);
        }
    }
}
async function processProject(projectPath, generator, formats) {
    const projectName = path.basename(projectPath);
    console.log(chalk.blue(`\n🔄 Processing: ${projectName}`));
    try {
        const result = await generator.generate({
            inputDir: projectPath,
            outputName: `${projectName}_report`,
            formats,
            template: 'report_template'
        });
        processedProjects.add(projectPath);
        console.log(chalk.green(`✅ Generated ${result.files.length} file(s):`));
        result.files.forEach(file => {
            console.log(chalk.gray(`   ${path.basename(file.path)}`));
        });
    }
    catch (error) {
        console.error(chalk.red(`❌ Failed to process ${projectName}:`), error);
    }
}
//# sourceMappingURL=watch.js.map