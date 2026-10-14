// UnsupportedNotice: what the shell shows in place of a report when it recognizes
// a result's method but has no viewer for it, rather than misreading the file as FEL

const NOTICE_ESCAPES = { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' };

function escapeNoticeText(text) {
  return String(text).replace(/[&<>"']/g, c => NOTICE_ESCAPES[c]);
}

// Markup for the notice. method is the detected id ('busted', ...) and supported
// lists the ids that do have a viewer.
export function unsupportedNotice(method, supported) {
  const name = escapeNoticeText(method.toUpperCase());
  const available = supported.map(id => escapeNoticeText(id.toUpperCase())).join(', ');
  return `
    <div class="hyphy-unsupported-notice" role="status">
      <div class="unsupported-notice-tag">Unsupported Result</div>
      <div class="unsupported-notice-title">No viewer for ${name} results yet</div>
      <div class="unsupported-notice-body">This file holds ${name} results. The viewer can display: ${available}.</div>
    </div>
  `;
}
