import { describe, it, expect, vi } from 'vitest';
import { parseGTF } from '../utils/excelParser';

// Mock XLSX since it's used in excelParser but we might not need it for GTF
vi.mock('xlsx', () => ({
  read: vi.fn(),
  utils: {
    sheet_to_json: vi.fn(),
  },
}));

describe('parseGTF', () => {
  it('should parse a GTF file content and return statistics', async () => {
    const gtfContent = `
# This is a comment
chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2";
chr1\tHAVANA\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2";
chr1\tHAVANA\texon\t12613\t12721\t.\t+\t.\tgene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2";
chr1\tHAVANA\texon\t13220\t14409\t.\t+\t.\tgene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2";
    `.trim();

    const file = new File([gtfContent], 'test.gtf', { type: 'text/plain' });

    const stats = await parseGTF(file);

    expect(stats.name).toBe('test.gtf');
    expect(stats.count).toBe(1); // 1 transcript
    // Exon 1: 12227 - 11869 + 1 = 359
    // Exon 2: 12721 - 12613 + 1 = 109
    // Exon 3: 14409 - 13220 + 1 = 1190
    // Total: 359 + 109 + 1190 = 1658
    expect(stats.totalLen).toBe(1658);
    expect(stats.meanLen).toBe(1658);
    expect(stats.maxLen).toBe(1658);
  });

  it('should handle multiple transcripts', async () => {
      const gtfContent = `
chr1\t.\texon\t100\t200\t.\t+\t.\ttranscript_id "T1";
chr1\t.\texon\t300\t400\t.\t+\t.\ttranscript_id "T2";
      `.trim();
      const file = new File([gtfContent], 'multi.gtf', { type: 'text/plain' });
      const stats = await parseGTF(file);
      
      expect(stats.count).toBe(2);
      // T1: 101, T2: 101. Total 202. Mean 101.
      expect(stats.totalLen).toBe(202);
      expect(stats.meanLen).toBe(101);
  });
  
  it('should return zeros for empty or invalid file', async () => {
       const file = new File([''], 'empty.gtf', { type: 'text/plain' });
       const stats = await parseGTF(file);
       expect(stats.count).toBe(0);
       expect(stats.totalLen).toBe(0);
  });
});
