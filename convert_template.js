import fs from 'fs';

let content = fs.readFileSync('templates/report_python.ejs', 'utf-8');

// Convert Jinja2 to EJS
content = content.replace(/\{\{([^}]+)\}\}/g, '<%= $1 %>');
content = content.replace(/\{%\s*if\s+([^%]+)\s*%\}/g, '<% if ($1) { %>');
content = content.replace(/\{%\s*else\s*%\}/g, '<% } else { %>');
content = content.replace(/\{%\s*endif\s*%\}/g, '<% } %>');
content = content.replace(/\{%\s*for\s+(\w+)\s+in\s+(\w+)\s*%\}/g, '<% $2.forEach(function($1) { %>');
content = content.replace(/\{%\s*endfor\s*%\}/g, '<% }); %>');
content = content.replace(/\{%\s*set\s+(\w+)\s*=\s*([^%]+)\s*%\}/g, '<% var $1 = $2; %>');
content = content.replace(/\{%\s*set\s+(\w+)\s*=\s*([^%]+)\s*if\s+([^%]+)\s*%\}/g, '<% var $1 = $3 ? $2 : []; %>');
content = content.replace(/\|\s*join\(['"]([^'"]+)['"]\)/g, ".join('$1')");
content = content.replace(/\.get\(([^,]+),\s*{}\)/g, '[$1] || {}');
content = content.replace(/\.get\(([^,]+),\s*\[\]\)/g, '[$1] || []');
content = content.replace(/\.get\(([^,]+),\s*['"]([^'"]*)['"]\)/g, '[$1] || "$2"');
content = content.replace(/\.get\(([^,]+)\)/g, '[$1]');
content = content.replace(/loop\.index/g, 'index + 1');
content = content.replace(/\|\s*length/g, '.length');
content = content.replace(/ or /g, ' || ');

fs.writeFileSync('templates/report_template.ejs', content);
console.log('Template converted successfully');
