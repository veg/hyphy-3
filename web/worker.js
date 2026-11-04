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
          self.postMessage({ type: 'log', stream: 'stdout', text: clean });
        }
      },
      printErr: function(text) {
        const clean = stripAnsi(text).trim();
        if (clean) {
          self.postMessage({ type: 'log', stream: 'stderr', text: clean });
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
    self.postMessage({ type: 'error', message: 'WASM initialization failed: ' + (err && err.message ? err.message : String(err)) });
  });

self.onmessage = async function(e) {
  const data = e.data || {};
  if (data.type === 'ping') {
    const t0 = performance.now();
    try {
      const mod = await initModule();
      const seq = data.sequence || '';
      const label = data.label || 'RED_OCTOBER_TARGET';
      const code = data.code || 'Universal';

      self.postMessage({
        type: 'log',
        stream: 'stdout',
        text: '>>> CONNING TOWER: Transmitting acoustic codon ping on target [' + label + '] (' + code + ')...'
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
      self.postMessage({
        type: 'error',
        message: err && err.message ? err.message : String(err)
      });
    }
  }
};
