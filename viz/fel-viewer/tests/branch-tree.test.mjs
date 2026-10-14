import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseNewick, lengthFrom } from '../src/core/newick.js';
import { branchTreeLayout } from '../src/components/branch-tree.js';
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
