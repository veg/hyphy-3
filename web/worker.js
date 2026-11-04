// Web Worker for HyPhy-3 WALRUS ("ONE PING ONLY") WebAssembly Engine
let walrusModulePromise = null;

function stripAnsi(str) {
  return String(str || '').replace(/\x1B\[[0-9;]*[a-zA-Z]/g, '');
}

function initModule() {
  if (!walrusModulePromise) {
    importScripts('hyphy3_walrus.js');
    walrusModulePromise = createWalrusModule({
      print: function(text) {
        const clean = stripAnsi(text).trim();
        if (clean) {
          self.postMessage({ type: 'log', stream: 'stdout', line: clean, text: clean });
        }
      },
      printErr: function(text) {
        const clean = stripAnsi(text).trim();
        if (clean) {
          self.postMessage({ type: 'log', stream: 'stderr', line: clean, text: clean });
        }
      }
    });
  }
  return walrusModulePromise;
}

// Pre-warm the WebAssembly module as soon as the worker spawns
initModule()
  .then(() => {
    self.postMessage({ type: 'ready' });
  })
  .catch((err) => {
    const errMsg = 'WASM initialization failed: ' + (err && err.message ? err.message : String(err));
    self.postMessage({ type: 'error', error: errMsg, message: errMsg });
  });

self.onmessage = async function(e) {
  const data = e.data || {};
  if (data.type === 'ping') {
    const t0 = performance.now();
    try {
      const mod = await initModule();
      const seq = data.seqContent || data.sequence || '';
      const label = data.label || 'RED_OCTOBER_TARGET';
      const code = data.geneticCode || data.code || 'Universal';

      const banner = '>>> CONNING TOWER: Transmitting acoustic codon ping on target [' + label + '] (' + code + ', ' + seq.length + ' bytes)...';
      self.postMessage({
        type: 'log',
        stream: 'stdout',
        line: banner,
        text: banner
      });

      const html = mod.ccall(
        'walrus_ping_run',
        'string',
        ['string', 'string', 'string'],
        [seq, label, code]
      );
      const err = mod.ccall('walrus_get_last_error', 'string', [], []);
      if (!html || err) {
        throw new Error(err || 'Empty HTML report returned from WALRUS WASM engine.');
      }
      const json = mod.ccall('walrus_get_last_json', 'string', [], []);
      const elapsedMs = performance.now() - t0;

      self.postMessage({
        type: 'result',
        html: html,
        json: json,
        label: label,
        elapsedMs: elapsedMs
      });
    } catch (err) {
      const errMsg = err && err.message ? err.message : String(err);
      self.postMessage({
        type: 'error',
        error: errMsg,
        message: errMsg
      });
    }
  } else if (data.type === 'render_json') {
    const t0 = performance.now();
    try {
      const mod = await initModule();
      const jsonStr = data.jsonContent || data.json || '';
      const label = data.label || 'WALRUS_JSON_DOSSIER';

      const banner = '>>> CONNING TOWER: Re-rendering interactive WALRUS HTML dossier from JSON [' + label + ']...';
      self.postMessage({
        type: 'log',
        stream: 'stdout',
        line: banner,
        text: banner
      });

      const html = mod.ccall(
        'walrus_render_html_from_json',
        'string',
        ['string'],
        [jsonStr]
      );
      const err = mod.ccall('walrus_get_last_error', 'string', [], []);
      if (!html || err) {
        throw new Error(err || 'Failed to render HTML from JSON payload.');
      }
      const elapsedMs = performance.now() - t0;

      self.postMessage({
        type: 'result',
        html: html,
        json: jsonStr,
        label: label,
        elapsedMs: elapsedMs
      });
    } catch (err) {
      const errMsg = err && err.message ? err.message : String(err);
      self.postMessage({
        type: 'error',
        error: errMsg,
        message: errMsg
      });
    }
  }
};
