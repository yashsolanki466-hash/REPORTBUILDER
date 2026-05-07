interface BatchOptions {
    rootDir: string;
    formats: string[];
    maxDepth: number;
    pattern: string;
}
export declare function batchGenerate(options: BatchOptions): Promise<void>;
export {};
//# sourceMappingURL=batch.d.ts.map