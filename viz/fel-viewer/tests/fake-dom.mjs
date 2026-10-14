// Just enough DOM for ResultsTable and Drawer to build and redraw: each selector maps
// to one stable child element, so what a component writes can be read back
export class FakeElement {
  constructor() {
    this.children = [];
    this.events = [];
    this.found = new Map();
    this.listeners = {};
    this.style = {};
    this.dataset = {};
    this.childNodes = [{ nodeValue: '' }];
    this.classList = { add() {}, remove() {}, toggle() {} };
    this.innerHTML = '';
    this.textContent = '';
    this.scrollTop = 0;
    this.clientHeight = 0;
    this.offsetHeight = 0;
  }
  appendChild(el) { this.children.push(el); return el; }
  prepend(el) { this.children.unshift(el); }
  before(el) { this.inserted = el; }
  remove() { this.removed = true; }
  replaceChildren() { this.children = []; }
  addEventListener(type, fn) { (this.listeners[type] ||= []).push(fn); }
  setAttribute(name, value) { (this.attributes ||= {})[name] = String(value); }
  getAttribute(name) { return this.attributes?.[name] ?? null; }
  dispatchEvent(e) { this.events.push(e); }
  fire(type, event = {}) { (this.listeners[type] || []).forEach(fn => fn(event)); }
  querySelector(selector) {
    if (!this.found.has(selector)) this.found.set(selector, new FakeElement());
    return this.found.get(selector);
  }
  querySelectorAll() { return []; }
}

export function installFakeDOM() {
  globalThis.document = { createElement: () => new FakeElement() };
  globalThis.window = { addEventListener() {}, removeEventListener() {} };
  globalThis.requestAnimationFrame = (fn) => fn();
}
