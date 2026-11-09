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
  if (data.type === 'stage0') {
    const t0 = performance.now();
    try {
      const mod = await initModule();
      const seq = data.seqContent || data.sequence || '';
      const label = data.label || 'STAGE0_INPUT';
      const code = data.geneticCode || data.code || 'Universal';
      const minCodons = (typeof data.minSenseCodons === 'number') ? data.minSenseCodons : 15;
      const minZ = (typeof data.minHomologyZ === 'number') ? data.minHomologyZ : 3.5;
      const stripGaps = data.stripGapsRealign ? 1 : 0;

      const jsonStr = mod.ccall(
        'walrus_stage0_run',
        'string',
        ['string', 'string', 'string', 'number', 'number', 'number'],
        [seq, label, code, minCodons, minZ, stripGaps]
      );
      const err = mod.ccall('walrus_get_last_error', 'string', [], []);
      if (!jsonStr || err) {
        throw new Error(err || 'Empty Stage 0 JSON returned from WALRUS WASM engine.');
      }
      const parsed = JSON.parse(jsonStr);
      const elapsedMs = performance.now() - t0;

      self.postMessage({
        type: 'stage0_result',
        stage0: parsed,
        json: jsonStr,
        label: label,
        elapsedMs: elapsedMs
      });
    } catch (err) {
      const errMsg = err && err.message ? err.message : String(err);
      self.postMessage({
        type: 'error',
        stage: 0,
        error: errMsg,
        message: errMsg
      });
    }
  } else if (data.type === 'stage1') {
    const t0 = performance.now();
    try {
      const mod = await initModule();
      const seq = data.seqContent || data.sequence || '';
      const label = data.label || 'STAGE1_COHORT';
      const code = data.geneticCode || data.code || 'Universal';
      const refineWindows = (data.refineIndelWindows === false) ? 0 : 1;
      const forceDeNovo = data.forceDeNovoRealign ? 1 : 0;

      const jsonStr = mod.ccall(
        'walrus_stage1_run',
        'string',
        ['string', 'string', 'string', 'number', 'number'],
        [seq, label, code, refineWindows, forceDeNovo]
      );
      const err = mod.ccall('walrus_get_last_error', 'string', [], []);
      if (!jsonStr || err) {
        throw new Error(err || 'Empty Stage 1 JSON returned from WALRUS WASM engine.');
      }
      const parsed = JSON.parse(jsonStr);
      const elapsedMs = performance.now() - t0;

      self.postMessage({
        type: 'stage1_result',
        stage1: parsed,
        json: jsonStr,
        label: label,
        elapsedMs: elapsedMs
      });
    } catch (err) {
      const errMsg = err && err.message ? err.message : String(err);
      self.postMessage({
        type: 'error',
        stage: 1,
        error: errMsg,
        message: errMsg
      });
    }
  } else if (data.type === 'stage2') {
    const t0 = performance.now();
    try {
      const mod = await initModule();
      const seq = data.seqContent || data.sequence || '';
      const treeNewick = data.treeNewick || '';
      const label = data.label || 'STAGE2_COHORT';
      const code = data.geneticCode || data.code || 'Universal';
      const enableSpr = (data.enableCoupledSpr === false) ? 0 : 1;
      const enableRecomb = (data.enableRecombScan === false) ? 0 : 1;

      const jsonStr = mod.ccall(
        'walrus_stage2_run',
        'string',
        ['string', 'string', 'string', 'string', 'number', 'number'],
        [seq, treeNewick, label, code, enableSpr, enableRecomb]
      );
      const err = mod.ccall('walrus_get_last_error', 'string', [], []);
      if (!jsonStr || err) {
        throw new Error(err || 'Empty Stage 2 JSON returned from WALRUS WASM engine.');
      }
      const parsed = JSON.parse(jsonStr);
      const elapsedMs = performance.now() - t0;

      self.postMessage({
        type: 'stage2_result',
        stage2: parsed,
        json: jsonStr,
        label: label,
        elapsedMs: elapsedMs
      });
    } catch (err) {
      const errMsg = err && err.message ? err.message : String(err);
      self.postMessage({
        type: 'error',
        stage: 2,
        error: errMsg,
        message: errMsg
      });
    }
  } else if (data.type === 'stage3') {
    const t0 = performance.now();
    try {
      const mod = await initModule();
      const seq = data.seqContent || data.sequence || '';
      const treeNewick = data.treeNewick || '';
      const label = data.label || 'STAGE3_COHORT';
      const code = data.geneticCode || data.code || 'Universal';

      const jsonStr = mod.ccall(
        'walrus_stage3_run',
        'string',
        ['string', 'string', 'string', 'string'],
        [seq, treeNewick, label, code]
      );
      const err = mod.ccall('walrus_get_last_error', 'string', [], []);
      if (!jsonStr || err) {
        throw new Error(err || 'Empty Stage 3 JSON returned from WALRUS WASM engine.');
      }
      const parsed = JSON.parse(jsonStr);
      const fullJson = mod.ccall('walrus_get_last_json', 'string', [], []);
      const fullHtml = mod.ccall('walrus_get_last_html', 'string', [], []);
      const elapsedMs = performance.now() - t0;

      self.postMessage({
        type: 'stage3_result',
        stage3: parsed,
        json: fullJson || jsonStr,
        html: fullHtml || '',
        label: label,
        elapsedMs: elapsedMs
      });
    } catch (err) {
      const errMsg = err && err.message ? err.message : String(err);
      self.postMessage({
        type: 'error',
        stage: 3,
        error: errMsg,
        message: errMsg
      });
    }
  } else if (data.type === 'ping') {
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
        ['string', 'string'],
        [jsonStr, label]
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
