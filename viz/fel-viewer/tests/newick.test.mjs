import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseNewick } from '../src/core/newick.js';

const leaves = node => node.children.length ? node.children.flatMap(leaves) : [node];
const find = (node, name) =>
  node.name === name ? node : node.children.map(c => find(c, name)).find(Boolean);

test('reads topology, leaf names and branch lengths', () => {
  const root = parseNewick('((a:0.1,b:2e-3):0.5,c:1);');
  assert.deepEqual(leaves(root).map(n => n.name), ['a', 'b', 'c']);
  assert.equal(find(root, 'a').branch_length, 0.1);
  assert.equal(find(root, 'b').branch_length, 0.002);
  assert.equal(root.children[0].branch_length, 0.5);
  assert.equal(root.branch_length, null);
});

test('strips {Tag} annotations into tag, before or after the length', () => {
  const root = parseNewick('((Pig{Group2}:0.46,Cow:0.59{Group2})Node3{Group2}:0.24,Rat{Group3}:0.16)');
  const pig = find(root, 'Pig');
  assert.equal(pig.tag, 'Group2');
  assert.equal(pig.branch_length, 0.46);
  assert.equal(find(root, 'Cow').tag, 'Group2');
  assert.equal(find(root, 'Cow').branch_length, 0.59);
  assert.equal(find(root, 'Node3').tag, 'Group2');
  assert.equal(find(root, 'Rat').tag, 'Group3');
  assert.equal(root.tag, null);
});

test("tags: 'keep' leaves the annotation in the name", () => {
  const root = parseNewick('(Pig{Group2}:0.46,Node5{Group4}:0.06)', { tags: 'keep' });
  assert.deepEqual(leaves(root).map(n => n.name), ['Pig{Group2}', 'Node5{Group4}']);
  assert.equal(root.children[0].tag, 'Group2');
});

test('quoted names keep delimiters, spaces and escaped quotes', () => {
  const root = parseNewick(`('Homo sapiens (1)':0.1,"a,b:c":0.2,'it''s':0.3)`);
  assert.deepEqual(leaves(root).map(n => n.name), ['Homo sapiens (1)', 'a,b:c', "it's"]);
  assert.equal(root.children[1].branch_length, 0.2);
});

test('keeps internal node labels and skips [comments]', () => {
  const root = parseNewick('((A,B)Node2[&support=90]:0.3,(C,D)Node5)Node1;');
  assert.equal(root.name, 'Node1');
  assert.deepEqual(root.children.map(n => n.name), ['Node2', 'Node5']);
  assert.equal(root.children[0].branch_length, 0.3);
  assert.equal(root.children[1].branch_length, null);
});

test("returns the v3 root named 'END' as the root", () => {
  const root = parseNewick('((Pig:0,Cow:0)Node3:0,Rat:0)END;');
  assert.equal(root.name, 'END');
  assert.equal(root.children.length, 2);
  assert.deepEqual(leaves(root).map(n => n.name), ['Pig', 'Cow', 'Rat']);
});

test('zero lengths are kept by default or replaced when configured', () => {
  const tree = '(a:0,b:0.0,c:0.2,d)';
  assert.deepEqual(parseNewick(tree).children.map(n => n.branch_length), [0, 0, 0.2, null]);
  assert.deepEqual(parseNewick(tree, { zeroLength: 'keep' }).children.map(n => n.branch_length), [0, 0, 0.2, null]);
  assert.deepEqual(parseNewick(tree, { zeroLength: 0.05 }).children.map(n => n.branch_length), [0.05, 0.05, 0.2, null]);
});

test('lengths that are not finite numbers read as missing', () => {
  const root = parseNewick('(a:abc,b:,c:Infinity,d:1e400,e:0.5x)');
  assert.deepEqual(root.children.map(n => n.branch_length), [null, null, null, null, null]);
});

test('empty input has no tree', () => {
  assert.equal(parseNewick(''), null);
  assert.equal(parseNewick(';'), null);
  assert.equal(parseNewick(undefined), null);
});
