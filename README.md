# NGS Report Automation

This tool automates the generation of the UNIGENOME Final Analysis Report HTML by parsing project directories.

## Prerequisites

- Python 3.8+
- The project directory must follow the standard structure (`01_Raw_Data`, `02_reference_genome_and_gff`, etc.)

## Setup

Since the system Python environment is managed, it is recommended to use a virtual environment.

1. **Create a virtual environment:**

   ```bash
   python3 -m venv .venv
   ```
2. **Install dependencies:**

   ```bash
   .venv/bin/pip install pandas jinja2 openpyxl
   ```

## Usage

### Generate report

Run the script using the python executable from the virtual environment:

```bash
.venv/bin/python generate_report.py --input <project_dir> --output <output_html>
```

### Per-project `metadata.json`

To avoid hardcoding project details, you can place a `metadata.json` file in the root of each project folder.

Example:

```json
{
  "project_id": "240540",
  "pi_name": "Dr Madhurjit Singh Rathore",
  "application": "WTA",
  "genome": "pearl_millet_v1.1",
  "client_name": "Client Name",
  "client_org": "Client Organization",
  "reference_organism": "Solanum lycopersicum (GCF_036512215.1)",
  "logo_path": "assets/logo.png"
}
```

**Notes**

- If present, `metadata.json` is used as the primary source for `project_id`, PI, application, etc.
- CLI arguments still override metadata (e.g. `--logo-path`, `--client-name`).
- If you don’t pass `--output` (or leave it as the default `report.html`), the report is written to the project root as:
  `<project_dir>/<project_id>_report.html`

### Batch generation

Generate reports for multiple project folders under a single root directory:

```bash
.venv/bin/python batch_generate_reports.py --root /path/to/results_root
```

Optional flags:

- `--max-depth 3` to scan deeper
- `--strict` to fail on missing required inputs
- `--dry-run` to preview the commands without running

### Watcher (polling)

Continuously monitor a root directory and generate a report once a project folder contains `Readme.txt`:

```bash
.venv/bin/python watch_reports.py --root /path/to/results_root --poll-seconds 30
```

This is a polling-based watcher (no OS-specific dependencies). It triggers report generation once per project folder.

### Arguments

- `--input`: (Required) Path to the root folder of the project (e.g., `Deliverables`).
- `--output`: (Optional) Path for the generated HTML file. Defaults to `report.html`.

## Testing

To verify the setup, you can generate dummy data and run the report against it:

1. **Generate dummy data:**

   ```bash
   .venv/bin/python create_dummy_data.py
   ```
   This creates a folder `Deliverables_Dummy`.
2. **Run report generation:**

   ```bash
   .venv/bin/python generate_report.py --input Deliverables_Dummy --output dummy_report.html
   ```
3. Open `dummy_report.html` in your browser to check the results.
