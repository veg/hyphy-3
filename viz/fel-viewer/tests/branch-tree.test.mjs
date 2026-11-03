import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseNewick, lengthFrom } from '../src/core/newick.js';
import { BranchTree, branchTreeLayout } from '../src/components/branch-tree.js';
import { FakeElement } from './fake-dom.mjs';
import { nextSignificantSite } from '../src/components/overview-brush.js';

const v3Newick = '((Pig{G2}:0,Cow{G2}:0,Horse:0)Node3{G2}:0,Rat{G3}:0,Mouse:0)END;';
const v3Nodes = {
  Pig: { branch_lengths: { full_adaptive: 0.2, baseline_mg94: 0.1 } },
  Cow: { branch_lengths: { full_adaptive: 0.4 } },
  Horse: { branch_lengths: { full_adaptive: 0.3 } },
  Node3: { branch_lengths: { full_adaptive: 0.1 } },
  Rat: { branch_lengths: { full_adaptive: 0.05 } },
  Mouse: { branch_lengths: { full_adaptive: 0.15 } },
  END: { type: 'internal' }
};
const byName = layout => Object.fromEntries(layout.nodes.map(n => [n.name, n]));

test('lengthFrom reads a finite branch_lengths entry or returns null', () => {
  assert.equal(lengthFrom(v3Nodes.Pig, 'baseline_mg94'), 0.1);
  assert.equal(lengthFrom(v3Nodes.Cow, 'baseline_mg94'), null);
  assert.equal(lengthFrom(v3Nodes.END, 'full_adaptive'), null);
  assert.equal(lengthFrom(undefined, 'full_adaptive'), null);
  assert.equal(lengthFrom({ branch_lengths: { x: NaN } }, 'x'), null);
});

test('depths come from node metadata, not the zero v3 newick lengths', () => {
  const layout = branchTreeLayout(parseNewick(v3Newick), { nodeMeta: v3Nodes, lengthKey: 'full_adaptive' });
  const n = byName(layout);
  assert.equal(layout.cladogram, false);
  assert.equal(n.END.depth, 0);
  assert.equal(n.Node3.depth, 0.1);
  assert.ok(Math.abs(n.Cow.depth - 0.5) < 1e-12);
  assert.ok(Math.abs(layout.maxDepth - 0.5) < 1e-12);
  assert.equal(n.Cow.tag, 'G2');
  assert.equal(n.Cow.parent, n.Node3);
});

test('ladderizes by leaf count and centres internal rows on their children', () => {
  const layout = branchTreeLayout(parseNewick(v3Newick), { nodeMeta: v3Nodes, lengthKey: 'full_adaptive' });
  assert.deepEqual(layout.leaves.map(l => l.name), ['Rat', 'Mouse', 'Pig', 'Cow', 'Horse']);
  assert.deepEqual(layout.leaves.map(l => l.row), [0, 1, 2, 3, 4]);
  const n = byName(layout);
  assert.equal(n.Node3.row, 3);
  assert.equal(n.END.row, 1.5);

  const asWritten = branchTreeLayout(parseNewick(v3Newick), { ladderize: false });
  assert.deepEqual(asWritten.leaves.map(l => l.name), ['Pig', 'Cow', 'Horse', 'Rat', 'Mouse']);
});

test('falls back to edge counts when the tree has no lengths', () => {
  const layout = branchTreeLayout(parseNewick(v3Newick));
  const n = byName(layout);
  assert.equal(layout.cladogram, true);
  assert.equal(n.Node3.depth, 1);
  assert.equal(n.Pig.depth, 2);
  assert.equal(layout.maxDepth, 2);
});

test('newick lengths are used when no length key is given', () => {
  const layout = branchTreeLayout(parseNewick('((a:1,b:2):0.5,c:3);'));
  const n = byName(layout);
  assert.equal(layout.cladogram, false);
  assert.equal(n.b.depth, 2.5);
  assert.equal(layout.maxDepth, 3);
});

test('nextSignificantSite steps forward and back with wrap-around', () => {
  const sites = [5, 12, 40];
  assert.equal(nextSignificantSite(sites, 0, 1), 5);
  assert.equal(nextSignificantSite(sites, 5, 1), 12);
  assert.equal(nextSignificantSite(sites, 40, 1), 5);
  assert.equal(nextSignificantSite(sites, 12, -1), 5);
  assert.equal(nextSignificantSite(sites, 5, -1), 40);
  assert.equal(nextSignificantSite(sites, 0, -1), 40);
  assert.equal(nextSignificantSite([], 3, 1), null);
});

// A BranchTree over fake DOM: canvases get a context that measures 6px per character, and
// each ResizeObserver is kept so a test can fire it
function treeHarness(width) {
  const observers = [];
  globalThis.window = {};
  globalThis.ResizeObserver = class {
    constructor(callback) { this.callback = callback; observers.push(this); }
    observe(el) { this.target = el; }
    disconnect() { this.disconnected = true; }
  };
  globalThis.document = {
    createElement: () => Object.assign(new FakeElement(), {
      append() {},
      getContext: () => new Proxy({}, {
        get: (target, key) => (key === 'measureText' ? (text) => ({ width: String(text).length * 6 }) : () => {}),
        set: () => true
      }),
      getBoundingClientRect: () => ({ left: 0, top: 0 })
    })
  };
  const container = new FakeElement();
  container.clientWidth = width;
  return { container, observers };
}

const emphasised = (tree, name) => tree.scene().texts.find(t => t.text === name).bold;

test('a selected branch stays emphasised through hover and pointer leave', () => {
  const { container } = treeHarness(400);
  const selected = [];
  const tree = new BranchTree(container, { lengthKey: 'full_adaptive', onSelect: name => selected.push(name) });
  tree.setTree(v3Newick, v3Nodes);
  tree.scene();
  const pig = tree.edges.find(e => e.node.name === 'Pig');
  tree.canvas.fire('mousemove', { clientX: pig.x1 - 2, clientY: pig.y });
  tree.canvas.fire('click');
  assert.deepEqual(selected, ['Pig']);

  const rat = tree.edges.find(e => e.node.name === 'Rat');
  tree.canvas.fire('mousemove', { clientX: rat.x1 - 2, clientY: rat.y });
  assert.equal(emphasised(tree, 'Rat'), true);
  assert.equal(emphasised(tree, 'Pig'), true);
  tree.canvas.fire('mouseleave');
  tree.highlight(null);
  assert.equal(emphasised(tree, 'Rat'), false);
  assert.equal(emphasised(tree, 'Pig'), true);

  tree.select('Cow');
  assert.equal(emphasised(tree, 'Pig'), false);
  assert.equal(emphasised(tree, 'Cow'), true);
  tree.setTree(v3Newick, v3Nodes);
  assert.equal(tree.selected.name, 'Cow');
});

test('leaf labels come from the accessor; hit testing still uses node names', () => {
  const { container } = treeHarness(400);
  const meta = Object.assign({}, v3Nodes, { Pig: { ...v3Nodes.Pig, original_name: 'Sus scrofa' } });
  const tree = new BranchTree(container, {
    lengthKey: 'full_adaptive',
    leafLabel: (node, m) => m?.original_name ?? node.name
  });
  tree.setTree(v3Newick, meta);
  const texts = tree.scene().texts.map(t => t.text);
  assert.ok(texts.includes('Sus scrofa'));
  assert.ok(!texts.includes('Pig'));
  assert.ok(tree.svg().includes('>Sus scrofa</text>'));
  tree.select('Pig');
  assert.equal(emphasised(tree, 'Sus scrofa'), true);
});

test('a tree built while hidden draws once its container gets a width', () => {
  const { container, observers } = treeHarness(0);
  const tree = new BranchTree(container, { lengthKey: 'full_adaptive' });
  tree.setTree(v3Newick, v3Nodes);
  assert.equal(tree.width, 0);
  assert.equal(observers[0].target, container);

  container.clientWidth = 320;
  observers[0].callback();
  assert.equal(tree.width, 320);
  assert.equal(tree.canvas.style.width, '320px');
  assert.equal(tree.edges.length, 6);

  tree.destroy();
  assert.equal(observers[0].disconnected, true);
});
