// Inspector Drawer: slides over the viewer from the right to show one row (a site,
// branch or segment) in detail. The shell owns the slide-over, backdrop, Esc key and
// the header's tag and badge; the method config supplies the body as tabs.
//
// options.tabs[] = {
//   id, label (markup), comment (placeholder text inside the empty panel),
//   availableWhen(data) -> false hides the tab for this result (default: always shown),
//   unavailable: muted line explaining a hidden tab,
//   mount(panel) -> view kept for the drawer's lifetime; its setData(data) runs on load,
//   render(panel, row, { data, threshold, key, view })
// }
// options.row(data, key) looks up the inspected row; options.header(row, ctx) returns
// { tag, tone: 'div' | 'pur' | 'neutral', badge }.

export function drawerTabButtons(tabs) {
  return tabs.map((tab, i) => `<button type="button" class="seg-btn${i === 0 ? ' active' : ''}" data-tab="${tab.id}" id="tab-btn-${tab.id}">${tab.label}</button>`)
    .join('\n          ');
}

export function drawerTabPanels(tabs) {
  return tabs.map((tab, i) => `<div id="drawer-${tab.id}-panel"${i === 0 ? '' : ' style="display: none;"'}>
          <!-- ${tab.comment} -->
        </div>`)
    .join('\n        ');
}

export class Drawer {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      threshold: 0.1,
      tabs: [],
      emptyTag: '—',
      row: () => null,
      header: () => ({ tag: '', tone: 'neutral', badge: '—' }),
      onClose: null
    }, options);

    this.threshold = this.options.threshold;
    this.data = null;
    this.currentKey = null;
    this.tabs = this.options.tabs;
    this.currentTab = this.tabs[0]?.id;
    this.isOpen = false;

    this.initDOM();
    this.initKeyboardListener();
  }

  initDOM() {
    // Backdrop
    this.backdrop = document.createElement('div');
    this.backdrop.className = 'hyphy-drawer-backdrop';
    this.backdrop.addEventListener('click', () => this.close());
    this.container.appendChild(this.backdrop);

    // Drawer container
    this.drawer = document.createElement('div');
    this.drawer.className = 'hyphy-site-drawer';
    this.drawer.innerHTML = `
      <div class="drawer-header">
        <div class="drawer-title-group">
          <span class="drawer-site-tag" id="drawer-site-tag">${this.options.emptyTag}</span>
          <span class="drawer-regime-badge" id="drawer-regime-badge">—</span>
        </div>
        <div class="segmented-control drawer-nav-tabs">
          ${drawerTabButtons(this.tabs)}
        </div>
        <button type="button" class="drawer-close-btn" id="drawer-close-btn" title="Close Inspector (Esc)">✕</button>
      </div>

      <div class="drawer-body" id="drawer-body">
        ${drawerTabPanels(this.tabs)}
      </div>
    `;

    this.container.appendChild(this.drawer);

    this.drawer.querySelector('#drawer-close-btn').addEventListener('click', () => this.close());
    this.siteTag = this.drawer.querySelector('#drawer-site-tag');
    this.regimeBadge = this.drawer.querySelector('#drawer-regime-badge');
    this.body = this.drawer.querySelector('#drawer-body');
    this.note = null;

    // Per-tab button, panel and long-lived view
    this.parts = new Map(this.tabs.map(tab => {
      const button = this.drawer.querySelector(`#tab-btn-${tab.id}`);
      const panel = this.drawer.querySelector(`#drawer-${tab.id}-panel`);
      button.addEventListener('click', () => this.setTab(tab.id));
      return [tab.id, { tab, button, panel, view: tab.mount?.(panel) }];
    }));
  }

  isAvailable(tab) {
    return !this.data || !tab.availableWhen || tab.availableWhen(this.data);
  }

  setTab(id) {
    const part = this.parts.get(id);
    if (!part || !this.isAvailable(part.tab)) return;

    this.currentTab = id;
    this.parts.forEach(({ button, panel }, tabId) => {
      button.classList.toggle('active', tabId === id);
      panel.style.display = tabId === id ? 'block' : 'none';
    });
    this.renderTab(part);
  }

  renderTab({ tab, panel, view }) {
    const row = this.currentKey != null && this.data ? this.options.row(this.data, this.currentKey) : null;
    if (!row) return;
    tab.render(panel, row, { data: this.data, threshold: this.threshold, key: this.currentKey, view });
  }

  initKeyboardListener() {
    this.onKeyDown = (e) => {
      if (e.key === 'Escape' && this.isOpen) {
        this.close();
      }
    };
    window.addEventListener('keydown', this.onKeyDown);
  }

  // Detach from the page when the viewer swaps in another method's components
  destroy() {
    window.removeEventListener('keydown', this.onKeyDown);
    this.backdrop.remove();
    this.drawer.remove();
  }

  // Hide the tabs this result has no data for, with one muted line saying why
  setData(data) {
    this.data = data;
    const notes = [];
    this.parts.forEach(({ tab, button, view }) => {
      const available = this.isAvailable(tab);
      if (!available) {
        button.style.display = 'none';
        if (tab.unavailable) notes.push(tab.unavailable);
      } else if (button.style.display) {
        button.style.removeProperty('display');
      }
      if (data) view?.setData?.(data);
    });
    this.setNote(notes.join(' '));

    const tab = this.availableTab(this.currentTab);
    if (tab && tab !== this.currentTab) this.setTab(tab);
  }

  // The requested tab when this result has it, else the first one it has
  availableTab(id) {
    const part = this.parts.get(id);
    if (part && this.isAvailable(part.tab)) return id;
    return this.tabs.find(tab => this.isAvailable(tab))?.id;
  }

  setNote(text) {
    if (!text) {
      this.note?.remove();
      this.note = null;
      return;
    }
    if (!this.note) {
      this.note = document.createElement('div');
      this.note.className = 'drawer-tab-note';
      this.body.prepend(this.note);
    }
    this.note.textContent = text;
  }

  setThreshold(p) {
    this.threshold = p;
    if (this.isOpen && this.currentKey != null) {
      this.inspectSite(this.currentKey, this.currentTab);
    }
  }

  open() {
    this.isOpen = true;
    this.backdrop.classList.add('active');
    this.drawer.classList.add('active');
  }

  close() {
    this.isOpen = false;
    this.backdrop.classList.remove('active');
    this.drawer.classList.remove('active');
    if (typeof this.options.onClose === 'function') {
      this.options.onClose();
    }
  }

  // Open on a row, keyed by the table's rowKey, at the given tab or the current one
  inspectSite(key, activeTab = null) {
    if (!this.data) return;
    const row = this.options.row(this.data, key);
    if (!row) return;

    this.currentKey = key;
    this.open();

    const { tag, tone, badge } = this.options.header(row, { data: this.data, threshold: this.threshold });
    this.siteTag.textContent = tag;
    this.regimeBadge.className = `drawer-regime-badge tag-${tone}`;
    this.regimeBadge.textContent = badge;

    this.setTab(this.availableTab(activeTab || this.currentTab));
  }
}
