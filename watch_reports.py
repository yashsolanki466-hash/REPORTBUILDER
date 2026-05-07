import argparse
import os
import subprocess
import sys
import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

def _has_readme(project_dir: str) -> bool:
    for name in ("Readme.txt", "README.txt", "readme.txt", "metadata.json"):
        if os.path.exists(os.path.join(project_dir, name)):
            return True
    return False

def _is_project_dir(path: str) -> bool:
    if not os.path.isdir(path):
        return False
    if os.path.exists(os.path.join(path, "metadata.json")):
        return True
    if _has_readme(path):
        return True
    required = [
        "01_Raw_Data",
        "02_reference_genome_and_gff",
        "03_transcript_assembly_gtf",
        "04_transcript_sequences_fasta",
        "05_differential_expression_analysis",
    ]
    return any(os.path.isdir(os.path.join(path, d)) for d in required)

class ProjectHandler(FileSystemEventHandler):
    def __init__(self, root, gen_script, strict, pdf):
        self.root = os.path.abspath(root)
        self.gen_script = gen_script
        self.strict = strict
        self.pdf = pdf
        self.seen = set()

    def on_created(self, event):
        self._check_event(event)

    def on_modified(self, event):
        self._check_event(event)

    def _check_event(self, event):
        if event.is_directory:
            return
        
        path = event.src_path
        filename = os.path.basename(path)
        if filename.lower() in ("readme.txt", "metadata.json"):
            project_dir = os.path.dirname(path)
            if project_dir not in self.seen and _is_project_dir(project_dir):
                self._generate_report(project_dir)

    def _generate_report(self, project_dir):
        print(f"[watch] Triggering report generation for: {project_dir}")
        cmd = [sys.executable, self.gen_script, "--input", project_dir]
        if self.strict:
            cmd.append("--strict")
        if self.pdf:
            cmd.append("--pdf")

        try:
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            print(res.stdout)
            if res.returncode == 0:
                self.seen.add(project_dir)
                print(f"[watch] Successfully generated report for {project_dir}")
            else:
                print(f"[watch] Failed to generate report for {project_dir}")
        except Exception as e:
            print(f"[watch] Error running generate_report: {e}")

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Watch a root folder and generate NGS reports automatically")
    p.add_argument("--root", required=True, help="Root directory containing projects")
    p.add_argument("--strict", action="store_true", help="Pass --strict to generate_report.py")
    p.add_argument("--pdf", action="store_true", help="Generate PDF reports")
    args = p.parse_args(argv)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    gen_script = os.path.join(script_dir, "generate_report.py")

    event_handler = ProjectHandler(args.root, gen_script, args.strict, args.pdf)
    observer = Observer()
    observer.schedule(event_handler, args.root, recursive=True)
    
    print(f"[watch] Starting observer on {args.root}...")
    observer.start()
    
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
