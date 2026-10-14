/**
 * Download plumbing for PNG, SVG and CSV exports. Works from file:// with no server.
 */

function clickDownload(href, filename) {
  const link = document.createElement('a');
  link.download = filename;
  link.href = href;
  link.click();
}

export function downloadBlob(content, type, filename) {
  const url = URL.createObjectURL(new Blob([content], { type }));
  clickDownload(url, filename);
  URL.revokeObjectURL(url);
}

export function downloadPNG(canvas, filename) {
  clickDownload(canvas.toDataURL('image/png'), filename);
}

export function downloadSVG(svg, filename) {
  downloadBlob(svg, 'image/svg+xml;charset=utf-8', filename);
}

// Rows are arrays of already-formatted cells; quote free text with csvQuote
export function downloadCSV(headers, rows, filename) {
  const csvContent = [headers, ...rows].map(row => row.join(',')).join('\n');
  downloadBlob(csvContent, 'text/csv;charset=utf-8', filename);
}

export function csvQuote(text) {
  return `"${(text || '').replace(/"/g, '""')}"`;
}
