#!/usr/bin/env python3
"""
Minimal, dependency-free DOCX text extractor.
Unzips the docx and pulls plain text from word/document.xml.
Used to align static narrative snippets with a reference DOCX.
"""

import json
import os
import re
import sys
import zipfile
from pathlib import Path
from typing import Dict, List

def extract_text_from_docx(docx_path: str) -> str:
    """Extract all text from a DOCX file (no external dependencies)."""
    if not os.path.exists(docx_path):
        raise FileNotFoundError(f"DOCX not found: {docx_path}")
    with zipfile.ZipFile(docx_path, "r") as zf:
        with zf.open("word/document.xml") as f:
            xml = f.read().decode("utf-8", errors="ignore")
    # Strip tags; keep raw text
    text = re.sub(r"</w:t[^>]*>", "", xml)
    text = re.sub(r"<w:t[^>]*>", "", text)
    text = re.sub(r"<[^>]+>", "", text)
    return text

def split_into_sections(text: str) -> Dict[str, str]:
    """Heuristic split into sections by capitalized headings (e.g., 1. Introduction)."""
    sections = {}
    # Very simple heuristic: split on lines that look like headings (digit + . + space + Capitalized)
    parts = re.split(r"\n(?=\d+\.\s+[A-Z])", text)
    for part in parts:
        part = part.strip()
        if not part:
            continue
        # Try to capture the heading as the key
        m = re.match(r"^(\d+\.\s+[^\n]+)", part)
        if m:
            key = m.group(1).replace(" ", "_").lower()
            sections[key] = part
        else:
            # Fallback: use as generic intro
            sections.setdefault("intro", sections.get("intro", "") + "\n\n" + part)
    return sections

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: python extract_docx_text.py <path-to-docx>\n")
        sys.exit(1)
    docx_path = sys.argv[1]
    out_path = Path(docx_path).with_suffix(".extracted.json")
    raw_text = extract_text_from_docx(docx_path)
    sections = split_into_sections(raw_text)
    # Save both raw and sectioned for manual alignment
    payload = {"raw": raw_text, "sections": sections}
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    print(f"Extracted text written to: {out_path}")

if __name__ == "__main__":
    main()
