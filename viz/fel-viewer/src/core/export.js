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

// Text for an SVG text node or attribute
export function escapeSvgText(text) {
  return String(text).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
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

// An SVG figure as a PNG at the screen's pixel density, for figures drawn as SVG markup
export function downloadSVGAsPNG(svg, width, height, filename) {
  const scale = window.devicePixelRatio || 1;
  const url = URL.createObjectURL(new Blob([svg], { type: 'image/svg+xml;charset=utf-8' }));
  const image = new Image();
  image.onload = () => {
    const canvas = document.createElement('canvas');
    canvas.width = width * scale;
    canvas.height = height * scale;
    const ctx = canvas.getContext('2d');
    ctx.scale(scale, scale);
    ctx.drawImage(image, 0, 0, width, height);
    URL.revokeObjectURL(url);
    downloadPNG(canvas, filename);
  };
  image.src = url;
}
