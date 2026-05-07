export declare function validateProjectDirectory(inputDir: string): Promise<{
    valid: boolean;
    error?: string;
}>;
export declare function formatFileSize(bytes: number): string;
export declare function findProjects(rootDir: string, maxDepth?: number, pattern?: string): Promise<string[]>;
//# sourceMappingURL=utils.d.ts.map