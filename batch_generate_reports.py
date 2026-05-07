import argparse
import os
import subprocess
import sys


def _is_project_dir(path: str) -> bool:
    if not os.path.isdir(path):
        return False
    if os.path.exists(os.path.join(path, "metadata.json")):
        return True
    for name in ("Readme.txt", "README.txt", "readme.txt"):
        if os.path.exists(os.path.join(path, name)):
            return True
    required = [
        "01_Raw_Data",
        "02_reference_genome_and_gff",
        "03_transcript_assembly_gtf",
        "04_transcript_sequences_fasta",
        "05_differential_expression_analysis",
    ]
    return any(os.path.isdir(os.path.join(path, d)) for d in required)


def _iter_project_dirs(root: str, max_depth: int):
    root = os.path.abspath(root)
    if _is_project_dir(root):
        yield root
        return

    root_depth = root.rstrip(os.sep).count(os.sep)
    for dirpath, dirnames, _ in os.walk(root):
        depth = dirpath.rstrip(os.sep).count(os.sep) - root_depth
        if depth > max_depth:
            dirnames[:] = []
            continue
        if dirpath == root:
            continue
        if _is_project_dir(dirpath):
            yield dirpath
            dirnames[:] = []


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Batch-generate NGS HTML reports for multiple project folders")
    p.add_argument("--root", required=True, help="Root directory containing one or more project folders")
    p.add_argument("--max-depth", type=int, default=2, help="Max directory depth to scan under root")
    p.add_argument("--strict", action="store_true", help="Pass --strict to generate_report.py")
    p.add_argument("--reference-organism", default=None, help="Override reference organism")
    p.add_argument("--logo-path", default=None, help="Override logo path")
    p.add_argument("--dry-run", action="store_true", help="Only print what would be generated")
    args = p.parse_args(argv)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    gen_script = os.path.join(script_dir, "generate_report.py")

    project_dirs = list(_iter_project_dirs(args.root, args.max_depth))
    if not project_dirs:
        print("No project directories found.")
        return 2

    failures = 0
    for proj in sorted(project_dirs):
        cmd = [sys.executable, gen_script, "--input", proj]
        if args.strict:
            cmd.append("--strict")
        if args.reference_organism is not None:
            cmd += ["--reference-organism", args.reference_organism]
        if args.logo_path is not None:
            cmd += ["--logo-path", args.logo_path]

        print(f"Generating report: {proj}")
        if args.dry_run:
            print(" ".join(cmd))
            continue

        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if res.returncode != 0:
            failures += 1
            print(res.stdout)

    if failures:
        print(f"Completed with failures: {failures}")
        return 1

    print("Completed successfully")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
