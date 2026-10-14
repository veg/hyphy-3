/**
 * The one Newick parser every tree view goes through.
 *
 * A node is { name, branch_length, tag, children }. HyPhy {Tag} annotations are
 * stripped from the label and kept in tag (null when absent), quoted names may
 * contain any character ('' escapes a quote), internal labels are kept, and
 * [comments] are skipped. branch_length is null when the length is missing or
 * not a finite number.
 *
 * Options:
 *   zeroLength: 'keep' (default) leaves zero lengths at 0; a number replaces them.
 *   tags: 'strip' (default) removes {Tag} from the name; 'keep' leaves it in the
 *         name as written. Either way tag holds the annotation.
 *
 * The v3 root is written as an internal node named 'END'; it is returned as the
 * root like any other label, and callers that draw branches skip it.
 */

const NEWICK_LABEL_END = new Set([',', '(', ')', ';']);

export function parseNewick(text, options = {}) {
  const { zeroLength = 'keep', tags = 'strip' } = options;
  const s = String(text || '');
  const top = { children: [] };
  const stack = [top];
  let i = 0;

  while (i < s.length) {
    const ch = s[i];
    if (ch === ';') break;
    if (ch === '(') {
      const node = newickNode();
      stack[stack.length - 1].children.push(node);
      stack.push(node);
      i++;
    } else if (ch === ')') {
      const node = stack.length > 1 ? stack.pop() : newickNode();
      i = readNewickLabel(s, i + 1, node, zeroLength, tags);
    } else if (ch === ',' || /\s/.test(ch)) {
      i++;
    } else {
      const leaf = newickNode();
      stack[stack.length - 1].children.push(leaf);
      i = readNewickLabel(s, i, leaf, zeroLength, tags);
    }
  }
  return top.children[0] || null;
}

function newickNode() {
  return { name: '', branch_length: null, tag: null, children: [] };
}

// Reads the name, {tag}, [comment] and :length that follow a node, in any order,
// and returns the index of the delimiter that ends them.
function readNewickLabel(s, i, node, zeroLength, tags) {
  let name = '';
  let length = null;
  while (i < s.length && !NEWICK_LABEL_END.has(s[i])) {
    const ch = s[i];
    if (ch === "'" || ch === '"') {
      i++;
      while (i < s.length) {
        if (s[i] === ch && s[i + 1] === ch) {
          name += ch;
          i += 2;
        } else if (s[i] === ch) {
          i++;
          break;
        } else {
          name += s[i++];
        }
      }
    } else if (ch === '{' || ch === '[') {
      const close = s.indexOf(ch === '{' ? '}' : ']', i);
      const end = close < 0 ? s.length : close;
      if (ch === '{') {
        node.tag = s.slice(i + 1, end);
        if (tags === 'keep') name += s.slice(i, end + 1);
      }
      i = end + 1;
    } else if (ch === ':') {
      let j = i + 1;
      while (j < s.length && !NEWICK_LABEL_END.has(s[j]) && s[j] !== '{' && s[j] !== '[') j++;
      length = s.slice(i + 1, j).trim();
      i = j;
    } else {
      name += ch;
      i++;
    }
  }

  node.name = name.trim();
  const value = length ? Number(length) : NaN;
  if (Number.isFinite(value)) {
    node.branch_length = value === 0 && zeroLength !== 'keep' ? zeroLength : value;
  }
  return i;
}
