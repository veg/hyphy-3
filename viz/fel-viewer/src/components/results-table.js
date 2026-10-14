// Results Table: one row per site, branch or segment
// Dieter Rams minimalist data table with category pills, live search, sortable columns and CSV export
//
// A method config describes the table as data:
//   rowKey: the row property that names a row ('site', 'branch', 'segment'). Rows
//     carry it as data-<rowKey>, and selection reports it.
//   rows(data): the parsed result's rows in their natural order.
//   columns: [{ id, header, width?, align?, title?, cell(row, view), csv?, sort? }]
//     cell returns the whole <td>; view is { data, threshold }.
//     csv: [{ header, value(row, view) }], the column's export fields at full precision.
//     A column without a cell is export-only. sort(row) is the value a header click
//     sorts by.
//   filters: [{ id, key, glyph?, label, count(stats, view), test(row, t), subtitle(n, t) }]
//     drawn as a segmented control; key names the button and badge ids, and stats is
//     data.getStats(t). The first filter names the method's finding, so its label is
//     rewritten on load; the others keep the markup's label.
//   blank: { filters, columns } drawn before data arrives (defaults to the method's own).
//     Every filter needs a blank slot with its id; blank filters the method lacks are
//     hidden once data loads.
//   defaultFilter, defaultSort: { column, dir: 'asc' | 'desc' }
//   search: { placeholder, match(row, q), noMatch: { title, desc } }
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
      threshold: 0.1,
      filterType: null,
      onSelect: null
    }, options);

    this.rowKey = this.options.rowKey;
    this.columns = this.options.columns;
    this.filters = this.options.filters;
    this.blankFilters = this.options.blank?.filters || this.filters;
    this.defaultFilter = this.options.defaultFilter || this.filters[0].id;
    this.defaultSort = this.options.defaultSort;
    this.sort = this.defaultSort;
    this.threshold = this.options.threshold;
    this.filterType = this.knownFilter(this.options.filterType);
    this.data = null;
    this.rows = [];
    this.rowsByKey = new Map();
    this.selectedKey = null;
    this.searchQuery = '';

    this.initDOM();
  }

  initDOM() {
    const { title, placeholder, search } = this.options;
    const blankColumns = this.options.blank?.columns || this.columns;
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

          <!-- Search Input -->
          <div class="candidates-search-box">
            <input type="text" class="candidates-search-input" id="search-input" placeholder="${search?.placeholder || 'Search...'}" />
          </div>

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

    this.wrapper.querySelector('#search-input').addEventListener('input', (e) => {
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

    this.wrapper.querySelector('#btn-export-csv').addEventListener('click', () => this.exportCSV());
  }

  setData(data) {
    this.data = data;
    this.rows = this.options.rows(data);
    this.rowsByKey = new Map(this.rows.map(r => [String(r[this.rowKey]), r]));
    this.render();
  }

  setThreshold(p) {
    this.threshold = p;
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

  // Header clicks sort ascending first, then flip
  sortBy(column) {
    if (!column?.sort) return;
    const dir = this.sort?.column === column.id && this.sort.dir === 'asc' ? 'desc' : 'asc';
    this.sort = { column: column.id, dir };
    this.render();
  }

  selectSite(key) {
    this.selectedKey = key;
    this.highlightRows(true);
  }

  highlightRows(scroll) {
    this.tableScroll.querySelectorAll(`tr[data-${this.rowKey}]`).forEach(tr => {
      const selected = tr.getAttribute(`data-${this.rowKey}`) === String(this.selectedKey);
      tr.classList.toggle('row-selected', selected);
      if (selected && scroll) tr.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    });
  }

  view() {
    return { data: this.data, threshold: this.threshold };
  }

  displayColumns() {
    return this.columns.filter(c => c.cell);
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

  // The rows the filter and search pass, in the current sort order
  visibleRows() {
    const filter = this.activeFilter();
    const q = this.searchQuery.toLowerCase().trim();
    const rows = this.rows.filter(r => (!q || this.options.search.match(r, q)) && filter.test(r, this.threshold));
    const column = this.sort && this.columns.find(c => c.id === this.sort.column);
    if (!column) return rows;
    const sign = this.sort.dir === 'desc' ? -1 : 1;
    return rows.sort((a, b) => {
      const x = column.sort(a);
      const y = column.sort(b);
      return sign * ((x > y) - (x < y));
    });
  }

  // The default order is unmarked; a header the viewer sorted by shows its direction
  sortMark(column) {
    const isDefault = this.sort === this.defaultSort;
    if (isDefault || this.sort.column !== column.id) return '';
    return this.sort.dir === 'asc' ? ' ↑' : ' ↓';
  }

  rowsMarkup(rows, view) {
    return rows.map(r => resultsRowMarkup(r, this.columns, view, this.rowKey, this.selectedKey === r[this.rowKey])).join('');
  }

  render() {
    if (!this.data) return;
    this.updateFilters();
    this.thead.innerHTML = resultsHeaderMarkup(this.columns, ' '.repeat(8), (c) => this.sortMark(c));

    const filter = this.activeFilter();
    const rows = this.visibleRows();
    this.subtitle.textContent = filter.subtitle(rows.length, this.threshold);

    if (rows.length === 0) {
      this.tableEl.style.display = 'none';
      this.emptyEl.style.display = 'block';
      this.renderEmptyState(filter);
      return;
    }

    this.tableEl.style.display = 'table';
    this.emptyEl.style.display = 'none';
    this.tbody.innerHTML = this.rowsMarkup(rows, this.view());
  }

  // A search that matches nothing says so; an empty filter gets the method's empty state
  renderEmptyState(filter) {
    const stats = this.data.getStats(this.threshold);
    const state = this.searchQuery.trim()
      ? Object.assign({ icon: '🔍' }, this.options.search.noMatch)
      : this.options.emptyState(this.threshold, { filter, stats, rows: this.rows });

    let html = state.icon ? `<div class="empty-icon">${state.icon}</div>` : '';
    html += `<div class="empty-title">${state.title}</div>`;
    if (state.desc) html += `<div class="empty-desc">${state.desc}</div>`;
    if (state.jump) {
      html += `<button type="button" class="btn-jump" data-filter="${state.jump.filter}">${state.jump.label}</button>`;
    }
    if (state.fallback?.length) {
      html += `
        <div class="empty-fallback">
          <div class="empty-fallback-label">Closest calls</div>
          <table class="candidates-table candidates-fallback">
            <thead>${resultsHeaderMarkup(this.columns, ' '.repeat(14))}</thead>
            <tbody>${this.rowsMarkup(state.fallback, this.view())}</tbody>
          </table>
        </div>`;
    }
    this.emptyEl.innerHTML = html;
    this.emptyEl.querySelector('.btn-jump')?.addEventListener('click', () => this.setFilterType(state.jump.filter));
  }

  // CSV of the current filter, search and sort
  exportCSV() {
    if (!this.data) return;
    const view = this.view();
    const fields = this.columns.flatMap(c => c.csv || []);
    const rows = this.visibleRows().map(r => fields.map(f => f.value(r, view)));
    downloadCSV(fields.map(f => f.header), rows, `${this.options.exportName}_${this.filterType}.csv`);
  }
}
