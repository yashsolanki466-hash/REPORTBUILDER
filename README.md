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

Run the script using the python executable from the virtual environment:

```bash
.venv/bin/python generate_report.py --input /path/to/project_directory --output final_report.html
```

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
