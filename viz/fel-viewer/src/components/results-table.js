// Results Table: one row per site, branch or segment
// Dieter Rams minimalist data table with category pills, live search, sortable columns and CSV export
//
// A method config describes the table as data:
//   rowKey: the row property that names a row ('site', 'branch', 'segment'). Rows
//     carry it as data-<rowKey>, and selection reports it.
//   rows(data): the parsed result's rows in their natural order.
//   columns: [{ id, header, width?, align?, title?, cell(row, view), csv?, sort?,
//              availableWhen?(capabilities), hiddenByDefault?, alwaysShown? }]
//     or columns(data) returning that array, resolved on every setData.
//     cell returns the whole <td>; view is { data, threshold, statistic }.
//     csv: [{ header, value(row, view) }], the column's export fields at full precision.
//     A column without a cell is export-only. sort(row) is the value a header click
//     sorts by. availableWhen(data.capabilities) false drops the column for this result.
//     hiddenByDefault columns start hidden; the column chooser shows or hides any column
//     with a cell that is not alwaysShown, and the export follows what is shown.
//   columnChooser: draw the column chooser (default: when an available column is hiddenByDefault)
//   virtualize: true or { minRows, overscan } draws only the rows in view once the
//     filtered rows reach minRows (default 2000); off by default.
//   filters: [{ id, key, glyph?, label, count(stats, view), test(row, t, view), subtitle(n, t),
//              sort? }]
//     drawn as a segmented control; key names the button and badge ids, and stats is
//     data.getStats(t). sort ({ column, dir }) orders the filter's rows until the
//     viewer picks a column, in place of defaultSort. The first filter names the
//     method's finding, so its label is rewritten on load; the others keep the
//     markup's label.
//   blank: { filters, columns } drawn before data arrives (defaults to the method's own).
//     Every filter needs a blank slot with its id; blank filters the method lacks are
//     hidden once data loads.
//   defaultFilter, defaultSort: { column, dir: 'asc' | 'desc' }
//   search: { placeholder, match(row, q), noMatch: { title, desc } }; null draws no search box
//   emptyState(t, view): when the filter leaves no rows; view is { filter, stats, rows }.
//     Returns { icon?, title, desc?, jump?: { filter, label }, fallback?: rows }.
//     The fallback rows are drawn dimmed under 'Closest calls'.
//   title, placeholder (the subtitle before data), exportName (the CSV file stem)

import { downloadCSV } from '../core/export.js';

// The n rows with the smallest value: the 'Closest calls' when nothing passes
export function closestCalls(rows, value, n = 5) {
  return [...rows].sort((a, b) => value(a) - value(b)).slice(0, n);
}

function filterText(filter) {
  return filter.glyph ? `${filter.glyph} ${filter.label}` : filter.label;
}

function resultsFilterMarkup(filter, active) {
  return `
            <button type="button" class="seg-btn${active ? ' active' : ''}" data-filter="${filter.id}" id="btn-filter-${filter.key}">
              ${filterText(filter)} <span class="badge-count" id="count-${filter.key}">0</span>
            </button>`;
}

// The <tr> of column headers, indented to sit inside <thead> at depth indent
export function resultsHeaderMarkup(columns, indent, mark = () => '') {
  const cells = columns.filter(c => c.cell).map(c => {
    const style = [c.width && `width: ${c.width}px;`, c.align && `text-align: ${c.align};`].filter(Boolean).join(' ');
    const attrs = (style ? ` style="${style}"` : '') + (c.title ? ` title="${c.title}"` : '');
    return `\n${indent}  <th${attrs}>${c.header}${mark(c)}</th>`;
  });
  return `\n${indent}<tr>${cells.join('')}\n${indent}</tr>\n${indent.slice(2)}`;
}

export function resultsRowMarkup(row, columns, view, rowKey, selected) {
  const cells = columns.filter(c => c.cell).map(c => `\n            ${c.cell(row, view)}`);
  return `
          <tr data-${rowKey}="${row[rowKey]}" class="${selected ? 'row-selected' : ''}">${cells.join('')}
          </tr>
        `;
}

function resultsSearchMarkup(search) {
  return `
          <!-- Search Input -->
          <div class="candidates-search-box">
            <input type="text" class="candidates-search-input" id="search-input" placeholder="${search.placeholder || 'Search...'}" />
          </div>
`;
}

// The columns a result offers: the config's list, or columns(data), less those whose
// availableWhen rejects the result's capabilities
export function resolveColumns(columns, data) {
  const list = typeof columns === 'function' ? columns(data) : columns;
  return list.filter(c => !c.availableWhen || c.availableWhen(data?.capabilities || {}));
}

// The slice of rows a virtual table draws at scroll offset top in a viewport of the
// given height, with overscan rows either side and the spacer heights around it
export function virtualWindow(top, height, rowHeight, total, overscan) {
  const first = Math.max(0, Math.floor(top / rowHeight) - overscan);
  const last = Math.min(total, Math.ceil((top + height) / rowHeight) + overscan);
  return { first, last, before: first * rowHeight, after: (total - last) * rowHeight };
}

function virtualSpacer(height, span) {
  return height > 0 ? `<tr class="virtual-spacer" aria-hidden="true"><td colspan="${span}" style="height: ${height}px;"></td></tr>` : '';
}

export function columnChooserItems(columns, shown) {
  return columns.map(c => `<label class="column-chooser-item"><input type="checkbox" data-column="${c.id}"${shown(c) ? ' checked' : ''} /> ${c.header}</label>`).join('');
}

export class ResultsTable {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      title: 'Results',
      placeholder: '',
      rowKey: 'site',
      rows: () => [],
      columns: [],
      filters: [],
      blank: null,
      defaultFilter: null,
      defaultSort: null,
      search: null,
      emptyState: null,
      exportName: 'hyphy_results',
      columnChooser: null,
      virtualize: null,
      threshold: 0.1,
      statistic: null, // which statistic the threshold applies to ('p' | 'q'), for cells and filters
      filterType: null,
      onSelect: null
    }, options);

    this.rowKey = this.options.rowKey;
    this.columns = Array.isArray(this.options.columns) ? this.options.columns : [];
    this.columnShown = new Map(); // the chooser's choices by column id, kept across loads
    this.virtualize = this.options.virtualize && Object.assign({ minRows: 2000, overscan: 12 }, this.options.virtualize);
    this.virtual = false;
    this.drawn = null; // the virtual window last drawn
    this.frame = null;
    this.rowHeight = 32;
    this.filters = this.options.filters;
    this.blankFilters = this.options.blank?.filters || this.filters;
    this.defaultFilter = this.options.defaultFilter || this.filters[0].id;
    this.defaultSort = this.options.defaultSort;
    this.sort = null; // the column the viewer sorted by, over the default order
    this.threshold = this.options.threshold;
    this.statistic = this.options.statistic;
    this.filterType = this.knownFilter(this.options.filterType);
    this.data = null;
    this.rows = [];
    this.rowsByKey = new Map();
    this.bodyRows = [];
    this.selectedKey = null;
    this.searchQuery = '';

    this.initDOM();
  }

  initDOM() {
    const { title, placeholder, search } = this.options;
    const blankColumns = this.options.blank?.columns || this.columns.filter(c => !c.hiddenByDefault);
    const defaultBlank = this.blankFilters.find(f => f.id === this.defaultFilter) || this.blankFilters[0];
    const filterButtons = this.blankFilters.map(f => resultsFilterMarkup(f, f.id === this.filterType)).join('');

    this.wrapper = document.createElement('div');
    this.wrapper.className = 'hyphy-candidates-card';
    this.wrapper.innerHTML = `
      <div class="candidates-header">
        <div class="candidates-title-group">
          <span class="candidates-title">${title}</span>
          <span class="candidates-subtitle" id="candidates-subtitle">${placeholder}</span>
        </div>
        <div class="candidates-controls">
          <!-- Category Filter Pills: ${defaultBlank.label} Default -->
          <div class="segmented-control category-filter-control">${filterButtons}
          </div>
${search ? resultsSearchMarkup(search) : ''}
          <!-- CSV Export -->
          <button type="button" class="plot-tool-btn" id="btn-export-csv" title="Export ${title} to CSV">Export CSV</button>
        </div>
      </div>

      <!-- Table Body Container -->
      <div class="candidates-table-scroll" id="table-scroll">
        <table class="candidates-table" id="candidates-table">
          <thead>${resultsHeaderMarkup(blankColumns, ' '.repeat(12))}</thead>
          <tbody id="candidates-tbody">
            <!-- Dynamic rows -->
          </tbody>
        </table>
        <div class="candidates-empty-state" id="empty-state" style="display: none;">
          <!-- Dynamic empty message -->
        </div>
      </div>
    `;

    this.container.appendChild(this.wrapper);

    this.tbody = this.wrapper.querySelector('#candidates-tbody');
    this.thead = this.wrapper.querySelector('#candidates-table thead');
    this.emptyEl = this.wrapper.querySelector('#empty-state');
    this.tableScroll = this.wrapper.querySelector('#table-scroll');
    this.tableEl = this.wrapper.querySelector('#candidates-table');
    this.subtitle = this.wrapper.querySelector('#candidates-subtitle');

    this.wrapper.querySelectorAll('.category-filter-control .seg-btn').forEach(btn => {
      btn.addEventListener('click', () => this.setFilterType(btn.dataset.filter));
    });

    this.wrapper.querySelector('#search-input')?.addEventListener('input', (e) => {
      this.searchQuery = e.target.value;
      this.render();
    });

    this.thead.addEventListener('click', (e) => {
      const th = e.target.closest('th');
      if (th) this.sortBy(this.displayColumns()[th.cellIndex]);
    });

    // Rows of the table and of the empty state's closest calls
    this.tableScroll.addEventListener('click', (e) => {
      const tr = e.target.closest(`tr[data-${this.rowKey}]`);
      if (!tr) return;
      const row = this.rowsByKey.get(tr.getAttribute(`data-${this.rowKey}`));
      this.selectedKey = row[this.rowKey];
      this.highlightRows(false);

      const action = e.target.closest('button')?.dataset?.action || 'summary';
      if (typeof this.options.onSelect === 'function') {
        this.options.onSelect(row[this.rowKey], row, action);
      }
    });

    this.exportBtn = this.wrapper.querySelector('#btn-export-csv');
    this.exportBtn.addEventListener('click', () => this.exportCSV());

    // A virtual table redraws its window when the scroll moves past it
    this.tableScroll.addEventListener('scroll', () => {
      if (!this.virtual || this.frame) return;
      this.frame = requestAnimationFrame(() => {
        this.frame = null;
        this.drawBody(false);
      });
    });
  }

  setData(data) {
    this.data = data;
    this.columns = resolveColumns(this.options.columns, data);
    this.rows = this.options.rows(data);
    this.rowsByKey = new Map(this.rows.map(r => [String(r[this.rowKey]), r]));
    this.updateChooser();
    this.render();
  }

  isShown(column) {
    return this.columnShown.get(column.id) ?? !column.hiddenByDefault;
  }

  setColumnShown(id, shown) {
    this.columnShown.set(id, shown);
    this.render();
  }

  // A compact menu of the drawable columns, beside the export button
  updateChooser() {
    const choosable = this.columns.filter(c => c.cell && !c.alwaysShown);
    if (!(this.options.columnChooser ?? choosable.some(c => c.hiddenByDefault))) {
      this.chooser?.remove();
      this.chooser = null;
      return;
    }
    if (!this.chooser) {
      this.chooser = document.createElement('details');
      this.chooser.className = 'column-chooser';
      this.chooser.innerHTML = `<summary class="plot-tool-btn" title="Show or hide columns">Columns</summary><div class="column-chooser-menu"></div>`;
      this.chooser.addEventListener('change', (e) => this.setColumnShown(e.target.dataset.column, e.target.checked));
      this.chooser.addEventListener('toggle', () => this.trackChooserMenu());
      this.exportBtn.before(this.chooser);
    }
    this.chooser.querySelector('.column-chooser-menu').innerHTML = columnChooserItems(choosable, c => this.isShown(c));
  }

  // The menu is fixed so the card's overflow cannot clip it: placed under its button's
  // right edge on open, and kept there while anything scrolls
  trackChooserMenu() {
    this.placeChooserMenu ||= () => {
      const button = this.chooser.querySelector('summary').getBoundingClientRect();
      const menu = this.chooser.querySelector('.column-chooser-menu');
      menu.style.top = `${button.bottom + 4}px`;
      menu.style.right = `${document.documentElement.clientWidth - button.right}px`;
    };
    if (this.chooser.open) {
      this.placeChooserMenu();
      window.addEventListener('scroll', this.placeChooserMenu, true);
    } else {
      window.removeEventListener('scroll', this.placeChooserMenu, true);
    }
  }

  setThreshold(p) {
    this.threshold = p;
    this.render();
  }

  setStatistic(statistic) {
    this.statistic = statistic;
    this.render();
  }

  // Unknown filters (another method's, or the API's) fall back to the default
  knownFilter(filterType) {
    return this.filters.some(f => f.id === filterType) ? filterType : this.defaultFilter;
  }

  setFilterType(filterType) {
    this.filterType = this.knownFilter(filterType);
    this.wrapper.querySelectorAll('.category-filter-control .seg-btn').forEach(btn => {
      btn.classList.toggle('active', btn.dataset.filter === this.filterType);
    });
    this.render();
  }

  // The order in force: the viewer's column, else the filter's own, else the default
  currentSort() {
    return this.sort || this.activeFilter()?.sort || this.defaultSort;
  }

  // Header clicks sort ascending first, then flip
  sortBy(column) {
    if (!column?.sort) return;
    const current = this.currentSort();
    const dir = current?.column === column.id && current.dir === 'asc' ? 'desc' : 'asc';
    this.sort = { column: column.id, dir };
    this.render();
  }

  selectSite(key) {
    this.selectedKey = key;
    if (this.virtual) this.scrollToRow(key);
    this.highlightRows(true);
  }

  // Bring a virtual table's row into the drawn window
  scrollToRow(key) {
    const i = this.bodyRows.findIndex(r => String(r[this.rowKey]) === String(key));
    if (i < 0) return;
    const top = i * this.rowHeight;
    const bottom = top + this.rowHeight + this.thead.offsetHeight;
    const scroll = this.tableScroll;
    if (top < scroll.scrollTop) scroll.scrollTop = top;
    else if (bottom > scroll.scrollTop + scroll.clientHeight) scroll.scrollTop = bottom - scroll.clientHeight;
    this.drawBody(false);
  }

  highlightRows(scroll) {
    this.tableScroll.querySelectorAll(`tr[data-${this.rowKey}]`).forEach(tr => {
      const selected = tr.getAttribute(`data-${this.rowKey}`) === String(this.selectedKey);
      tr.classList.toggle('row-selected', selected);
      if (selected && scroll) tr.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    });
  }

  view() {
    return { data: this.data, threshold: this.threshold, statistic: this.statistic };
  }

  // The columns drawn: those with a cell that are not hidden
  displayColumns() {
    return this.columns.filter(c => c.cell && this.isShown(c));
  }

  // The drawn columns plus the export-only ones
  shownColumns() {
    return this.columns.filter(c => !c.cell || this.isShown(c));
  }

  activeFilter() {
    return this.filters.find(f => f.id === this.filterType);
  }

  // Badge counts follow the threshold; blank filters this method lacks are hidden
  updateFilters() {
    const view = this.view();
    const stats = this.data.getStats(this.threshold);
    for (const blank of this.blankFilters) {
      const filter = this.filters.find(f => f.id === blank.id);
      const btn = this.wrapper.querySelector(`#btn-filter-${blank.key}`);
      btn.querySelector('.badge-count').textContent = (filter || blank).count(stats, view);
      btn.style.display = filter ? '' : 'none';
    }
    const head = this.filters[0];
    this.wrapper.querySelector(`#btn-filter-${head.key}`).childNodes[0].nodeValue = `${filterText(head)} `;
  }

  // The search query, or nothing for a table without search
  searchText() {
    return this.options.search ? this.searchQuery.toLowerCase().trim() : '';
  }

  // The rows the filter and search pass, in the current sort order
  visibleRows() {
    const filter = this.activeFilter();
    const view = this.view();
    const q = this.searchText();
    const rows = this.rows.filter(r => (!q || this.options.search.match(r, q)) && filter.test(r, this.threshold, view));
    const sort = this.currentSort();
    const column = sort && this.columns.find(c => c.id === sort.column);
    if (!column) return rows;
    const sign = sort.dir === 'desc' ? -1 : 1;
    return rows.sort((a, b) => {
      const x = column.sort(a);
      const y = column.sort(b);
      return sign * ((x > y) - (x < y));
    });
  }

  // The default order is unmarked; a header the viewer sorted by shows its direction
  sortMark(column) {
    if (!this.sort || this.sort.column !== column.id) return '';
    return this.sort.dir === 'asc' ? ' ↑' : ' ↓';
  }

  rowsMarkup(rows, view) {
    const columns = this.shownColumns();
    return rows.map(r => resultsRowMarkup(r, columns, view, this.rowKey, this.selectedKey === r[this.rowKey])).join('');
  }

  render() {
    if (!this.data) return;
    this.updateFilters();
    this.thead.innerHTML = resultsHeaderMarkup(this.shownColumns(), ' '.repeat(8), (c) => this.sortMark(c));

    const filter = this.activeFilter();
    const rows = this.visibleRows();
    this.subtitle.textContent = filter.subtitle(rows.length, this.threshold);

    this.fallbackRows = [];
    if (rows.length === 0) {
      this.tableEl.style.display = 'none';
      this.emptyEl.style.display = 'block';
      this.renderEmptyState(filter);
      return;
    }

    this.tableEl.style.display = 'table';
    this.emptyEl.style.display = 'none';
    this.bodyRows = rows;
    this.virtual = Boolean(this.virtualize) && rows.length >= this.virtualize.minRows;
    this.drawBody(true);
  }

  // All rows, or for a virtual table the window in view between two spacer rows. The
  // row height is measured from the first drawn row, redrawing once if it was off.
  drawBody(force) {
    const view = this.view();
    if (!this.virtual) {
      this.tbody.innerHTML = this.rowsMarkup(this.bodyRows, view);
      return;
    }
    const scroll = this.tableScroll;
    const w = virtualWindow(scroll.scrollTop, scroll.clientHeight, this.rowHeight, this.bodyRows.length, this.virtualize.overscan);
    if (!force && this.drawn?.first === w.first && this.drawn?.last === w.last) return;
    this.drawn = w;
    const span = this.displayColumns().length;
    this.tbody.innerHTML = virtualSpacer(w.before, span)
      + this.rowsMarkup(this.bodyRows.slice(w.first, w.last), view)
      + virtualSpacer(w.after, span);

    const measured = this.tbody.querySelector(`tr[data-${this.rowKey}]`)?.offsetHeight;
    if (measured && Math.abs(measured - this.rowHeight) > 0.5) {
      this.rowHeight = measured;
      this.drawBody(true);
    }
  }

  // A search that matches nothing says so; an empty filter gets the method's empty state
  renderEmptyState(filter) {
    const stats = this.data.getStats(this.threshold);
    const state = this.searchText()
      ? Object.assign({ icon: '🔍' }, this.options.search.noMatch)
      : this.options.emptyState(this.threshold, { filter, stats, rows: this.rows });

    let html = state.icon ? `<div class="empty-icon">${state.icon}</div>` : '';
    html += `<div class="empty-title">${state.title}</div>`;
    if (state.desc) html += `<div class="empty-desc">${state.desc}</div>`;
    if (state.jump) {
      html += `<button type="button" class="btn-jump" data-filter="${state.jump.filter}">${state.jump.label}</button>`;
    }
    if (state.fallback?.length) {
      this.fallbackRows = state.fallback;
      html += `
        <div class="empty-fallback">
          <div class="empty-fallback-label">Closest calls</div>
          <table class="candidates-table candidates-fallback">
            <thead>${resultsHeaderMarkup(this.shownColumns(), ' '.repeat(14))}</thead>
            <tbody>${this.rowsMarkup(state.fallback, this.view())}</tbody>
          </table>
        </div>`;
    }
    this.emptyEl.innerHTML = html;
    this.emptyEl.querySelector('.btn-jump')?.addEventListener('click', () => this.setFilterType(state.jump.filter));
  }

  // The rows the drawer steps through: the table's, or the closest calls when it is empty
  stepRows() {
    const rows = this.visibleRows();
    return rows.length ? rows : (this.fallbackRows || []);
  }

  // CSV of the current filter, search and sort
  exportCSV() {
    if (!this.data) return;
    const view = this.view();
    const fields = this.shownColumns().flatMap(c => c.csv || []);
    const rows = this.visibleRows().map(r => fields.map(f => f.value(r, view)));
    downloadCSV(fields.map(f => f.header), rows, `${this.options.exportName}_${this.filterType}.csv`);
  }
}
