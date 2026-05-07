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
export declare class ReportGenerator {
    private options;
    private browser;
    constructor(options: GeneratorOptions);
    init(): Promise<void>;
    close(): Promise<void>;
    generate(options: GenerateOptions): Promise<GenerateResult>;
    private generatePDF;
    private generateDOCX;
    private parseProjectData;
}
export {};
//# sourceMappingURL=reportGenerator.d.ts.map