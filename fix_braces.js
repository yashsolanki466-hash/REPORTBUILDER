import fs from 'fs';

let content = fs.readFileSync('templates/report_template.ejs', 'utf-8');

// Fix remaining {{ }} patterns that weren't converted
content = content.replace(/\{\{/g, '<%= ');
content = content.replace(/\}\}/g, ' %>');

// Fix .get() method calls to bracket notation
content = content.replace(/\.get\(([^,]+)\)/g, '[$1]');
content = content.replace(/\.get\(([^,]+),\s*{\}\)/g, '[$1] || {}');
content = content.replace(/\.get\(([^,]+),\s*\[\]\)/g, '[$1] || []');
content = content.replace(/\.get\(([^,]+),\s*['"]([^'"]*)['"]\)/g, '[$1] || "$2"');

// Fix .get() with closing parenthesis issues
content = content.replace(/\[([^\]]+)\)\]/g, '[$1]');

fs.writeFileSync('templates/report_template.ejs', content);
console.log('Fixed braces and .get() calls');
