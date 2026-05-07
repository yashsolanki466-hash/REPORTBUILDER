import fs from 'fs-extra';
import path from 'path';
export async function validateProjectDirectory(inputDir) {
    try {
        const stats = await fs.stat(inputDir);
        if (!stats.isDirectory()) {
            return { valid: false, error: `Path is not a directory: ${inputDir}` };
        }
        // Check for expected NGS project structure
        const expectedFolders = [
            '01_Raw_Data',
            '05_differential_expression_analysis'
        ];
        const contents = await fs.readdir(inputDir);
        const hasExpectedStructure = expectedFolders.some(folder => contents.includes(folder));
        if (!hasExpectedStructure) {
            console.warn(`Warning: Directory doesn't have standard NGS structure, but will attempt to process anyway.`);
        }
        return { valid: true };
    }
    catch (error) {
        return { valid: false, error: `Cannot access directory: ${inputDir}` };
    }
}
export function formatFileSize(bytes) {
    if (bytes === 0)
        return '0 B';
    const k = 1024;
    const sizes = ['B', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(1)) + ' ' + sizes[i];
}
export async function findProjects(rootDir, maxDepth = 2, pattern = '*') {
    const projects = [];
    async function scan(dir, depth) {
        if (depth > maxDepth)
            return;
        try {
            const entries = await fs.readdir(dir, { withFileTypes: true });
            for (const entry of entries) {
                if (entry.isDirectory()) {
                    const fullPath = path.join(dir, entry.name);
                    // Check if this looks like a project directory
                    const hasProjectStructure = await checkProjectStructure(fullPath);
                    if (hasProjectStructure) {
                        projects.push(fullPath);
                    }
                    else {
                        // Recurse into subdirectories
                        await scan(fullPath, depth + 1);
                    }
                }
            }
        }
        catch (error) {
            // Skip directories we can't read
        }
    }
    await scan(rootDir, 0);
    return projects;
}
async function checkProjectStructure(dir) {
    try {
        const contents = await fs.readdir(dir);
        // Check for characteristic NGS folders
        const ngsMarkers = [
            '01_Raw_Data',
            '05_differential_expression_analysis',
            '06_Significant_DGE_GO',
            'Readme.txt'
        ];
        return ngsMarkers.some(marker => contents.includes(marker));
    }
    catch {
        return false;
    }
}
//# sourceMappingURL=utils.js.map