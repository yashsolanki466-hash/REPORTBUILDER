# NGS Report Generator CLI

A high-performance command-line tool for generating professional NGS analysis reports in HTML, PDF, and DOCX formats. Built with Node.js and Puppeteer for perfect PDF rendering using the Chrome engine.

## Features

- **Perfect PDF Quality**: Uses Chrome browser engine via Puppeteer for pixel-perfect PDF generation
- **Multiple Output Formats**: HTML, PDF, and DOCX in a single command
- **Batch Processing**: Generate reports for multiple projects at once
- **File Watching**: Auto-generate reports when new projects are added
- **Fast & Lightweight**: No GUI overhead, pure command-line efficiency
- **Professional Templates**: Clean, scientific report styling

## Installation

```bash
npm install
npm run build
```

## Usage

### Generate Single Report

```bash
node dist/cli.js generate -i "C:/Users/pooja/Documents/260018" -o report_name
```

Options:
- `-i, --input <path>`: Project directory path (required)
- `-o, --output <name>`: Output filename without extension (default: "report")
- `-f, --formats <list>`: Comma-separated formats: html,pdf,docx (default: all)
- `--project-id <id>`: Override project ID
- `--pi-name <name>`: Override PI name
- `--client-name <name>`: Override client name

### Batch Generate Reports

```bash
node dist/cli.js batch -r "C:/Users/pooja/Documents/projects_root"
```

Options:
- `-r, --root <path>`: Root directory containing project folders (required)
- `-f, --formats <list>`: Output formats (default: html,pdf,docx)
- `--max-depth <n>`: Maximum directory depth to scan (default: 2)

### Watch Mode (Auto-Generation)

```bash
node dist/cli.js watch -r "C:/Users/pooja/Documents/projects_root"
```

Automatically generates reports when new project folders are detected.

Options:
- `-r, --root <path>`: Directory to watch (required)
- `--interval <ms>`: Polling interval in milliseconds (default: 5000)

## Examples

### Generate HTML and PDF only:
```bash
node dist/cli.js generate -i "C:/Users/pooja/Documents/260018" -o my_report -f html,pdf
```

### Batch process all projects:
```bash
node dist/cli.js batch -r "C:/Users/pooja/Documents" --max-depth 2
```

### Watch for new projects:
```bash
node dist/cli.js watch -r "C:/Users/pooja/Documents/deliverables" -f pdf
```

## Why This Beats Python/WeasyPrint

1. **Better PDF Quality**: Chrome engine renders CSS perfectly vs WeasyPrint's limitations
2. **No Font Issues**: Uses system fonts correctly
3. **JavaScript Support**: Can execute JS in templates for dynamic charts
4. **Modern CSS**: Full support for flexbox, grid, CSS variables
5. **Faster**: No Python startup overhead, optimized Node.js performance

## Project Structure Detection

The CLI automatically detects NGS project folders containing:
- `01_Raw_Data/`
- `05_differential_expression_analysis/`
- `Readme.txt`

## Output Files

Reports are saved in the input project directory:
- `{project_id}_report.html`
- `{project_id}_report.pdf`
- `{project_id}_report.docx`

## Templates

Templates are EJS files in the `templates/` folder:
- `report_template.ejs`: Main scientific report template

To customize, edit the template and re-run the CLI.

## Requirements

- Node.js 18+
- Chrome/Chromium (auto-installed with Puppeteer)

## License

MIT
