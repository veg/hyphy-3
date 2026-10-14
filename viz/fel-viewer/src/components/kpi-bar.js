// KPI Summary Bar & Reactive Threshold Controller
// Dieter Rams / Braun Industrial aesthetic: clean cards, tactile presets, reactive slider
//
// A method config describes the bar as data:
//   cards: [{ id, label, value(stats, view), sub(stats, view), tone, widget?, ... }]
//     view is { data, threshold, statistic }; stats is data.getStats(threshold).
//     tone colours the value: 'accent' | 'purifying' | 'ink' | 'muted'.
//     widget(el, stats, view) draws into a slot below the sub line.
//     Optional: placeholder (the sub line before any data), comment (the HTML
//     comment above the card), cardClass, and valueId/subId to override the
//     default kpi-val-<id>/kpi-sub-<id> element ids.
//   blank: the cards drawn before data arrives, one per card (defaults to cards).
//     Only their markup is used; label and tone switch to the cards' on setData.
//   controller: { title, statistic, direction: 'le' | 'ge', min, max, step,
//     presets, seed, format(t) }; seed is the starting cutoff, and statistic the
//     one the cutoff applies to until the shell switches it (options.statistic).

const KPI_TONE_CLASSES = {
  accent: 'kpi-value kpi-text-orange',
  purifying: 'kpi-value kpi-text-blue',
  ink: 'kpi-value',
  muted: 'kpi-value kpi-text-muted'
};

function kpiCardMarkup(card, index) {
  const cardClass = card.cardClass ? `kpi-card ${card.cardClass}` : 'kpi-card';
  const widget = card.widget ? `
          <div class="kpi-widget"></div>` : '';
  return `
        <!-- ${index + 1}. ${card.comment || `${card.label} Card`} -->
        <div class="${cardClass}" id="kpi-${card.id}">
          <div class="kpi-label"></div>
          <div class="${KPI_TONE_CLASSES[card.tone || 'ink']}" id="kpi-val-${card.valueId || card.id}">—</div>
          <div class="kpi-sub" id="kpi-sub-${card.subId || card.id}">${card.placeholder || '—'}</div>${widget}
        </div>
`;
}

function thresholdControllerMarkup(controller, threshold, statistic, index) {
  const { format } = controller;
  const pill = (p) => {
    const active = Math.abs(p - threshold) < 0.002 ? ' active' : '';
    return `
            <button type="button" class="preset-pill${active}" data-p="${format(p)}">${thresholdLabel(controller, p, statistic)}</button>`;
  };
  return `
        <!-- ${index + 1}. Interactive Threshold Controller -->
        <div class="kpi-card kpi-card-control" id="kpi-controller">
          <div class="kpi-control-header">
            <span class="kpi-label">${controller.title}</span>
            <span class="kpi-thresh-display" id="kpi-thresh-display">${thresholdLabel(controller, threshold, statistic)}</span>
          </div>
          <div class="kpi-slider-row">
            <input type="range" class="kpi-slider" id="kpi-thresh-slider" min="${format(controller.min)}" max="${format(controller.max)}" step="${format(controller.step)}" value="${format(threshold)}" />
          </div>
          <div class="kpi-presets-row">${controller.presets.map(pill).join('')}
          </div>
        </div>
      `;
}

// 'p ≤ 0.05', 'PP ≥ 0.90': the cutoff as the controller and its presets show it
export function thresholdLabel(controller, t, statistic = controller.statistic) {
  const op = controller.direction === 'ge' ? '≥' : '≤';
  return `${statistic} ${op} ${controller.format(t)}`;
}

export class KpiBar {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      cards: [],
      blank: null,
      controller: null,
      threshold: null, // the controller's seed when unset
      statistic: null, // the controller's statistic when unset
      onThresholdChange: null
    }, options);

    this.cards = this.options.cards;
    this.blank = this.options.blank || this.cards;
    this.controller = this.options.controller;
    this.threshold = this.options.threshold ?? this.controller.seed;
    this.statistic = this.options.statistic ?? this.controller.statistic;
    this.data = null;
    this.initDOM();
  }

  initDOM() {
    this.el = document.createElement('div');
    this.el.className = 'hyphy-kpi-bar';
    this.el.innerHTML = `
      <div class="kpi-cards-grid">${this.blank.map(kpiCardMarkup).join('')}${
        thresholdControllerMarkup(this.controller, this.threshold, this.statistic, this.blank.length)}</div>
    `;

    this.container.appendChild(this.el);

    this.cardEls = this.blank.map(card => {
      const el = this.el.querySelector(`#kpi-${card.id}`);
      const label = el.querySelector('.kpi-label');
      label.textContent = card.label;
      return {
        label,
        value: el.querySelector('.kpi-value'),
        sub: el.querySelector('.kpi-sub'),
        widget: el.querySelector('.kpi-widget')
      };
    });

    // Bind slider & preset listeners
    this.slider = this.el.querySelector('#kpi-thresh-slider');
    this.display = this.el.querySelector('#kpi-thresh-display');
    this.presetBtns = this.el.querySelectorAll('.preset-pill');

    this.slider.addEventListener('input', (e) => {
      const val = parseFloat(e.target.value);
      this.setThreshold(val, true);
    });

    this.presetBtns.forEach(btn => {
      btn.addEventListener('click', () => {
        const val = parseFloat(btn.dataset.p);
        this.slider.value = val;
        this.setThreshold(val, true);
      });
    });
  }

  setThreshold(p, triggerCallback = false) {
    this.threshold = p;
    this.display.textContent = thresholdLabel(this.controller, p, this.statistic);

    // Update active preset button highlight
    this.presetBtns.forEach(btn => {
      const btnP = parseFloat(btn.dataset.p);
      if (Math.abs(btnP - p) < 0.002) {
        btn.classList.add('active');
      } else {
        btn.classList.remove('active');
      }
    });

    this.updateStats();

    if (triggerCallback && typeof this.options.onThresholdChange === 'function') {
      this.options.onThresholdChange(p);
    }
  }

  // Relabel the cutoff and its presets for the statistic it now applies to
  setStatistic(statistic) {
    this.statistic = statistic;
    this.display.textContent = thresholdLabel(this.controller, this.threshold, statistic);
    this.presetBtns.forEach(btn => {
      btn.textContent = thresholdLabel(this.controller, parseFloat(btn.dataset.p), statistic);
    });
    this.updateStats();
  }

  setData(data) {
    this.data = data;
    this.cards.forEach((card, i) => {
      this.cardEls[i].label.textContent = card.label;
      this.cardEls[i].value.className = KPI_TONE_CLASSES[card.tone || 'ink'];
    });
    this.updateStats();
  }

  updateStats() {
    if (!this.data) return;
    const stats = this.data.getStats(this.threshold);
    const view = { data: this.data, threshold: this.threshold, statistic: this.statistic };

    this.cards.forEach((card, i) => {
      const els = this.cardEls[i];
      els.value.textContent = card.value(stats, view);
      els.sub.textContent = card.sub(stats, view);
      if (card.widget) card.widget(els.widget, stats, view);
    });
  }
}
