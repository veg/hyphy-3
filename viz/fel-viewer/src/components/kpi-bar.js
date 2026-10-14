// KPI Summary Bar & Reactive Threshold Controller
// Dieter Rams / Braun Industrial aesthetic: clean cards, tactile presets, reactive slider
//
// A method config describes the bar as data:
//   cards: [{ id, label, value(stats, view), sub(stats, view), tone, widget?, ... }]
//     view is { data, threshold, statistic }; stats is data.getStats(threshold).
//     tone colours the value: 'accent' | 'purifying' | 'ink' | 'muted', or
//     tone(stats, view) for a card whose emphasis follows the cutoff.
//     widget(el, stats, view) draws into a slot below the sub line.
//     Optional: placeholder (the sub line before any data), comment (the HTML
//     comment above the card), cardClass, valueId/subId to override the
//     default kpi-val-<id>/kpi-sub-<id> element ids, and available(data) for a
//     card left out when the result cannot fill it.
//   blank: the cards drawn before data arrives (defaults to cards). When their ids
//     match the cards' one for one, setData relabels them in place; otherwise it
//     replaces them with the cards' own markup.
//   controller: { title, statistic, direction: 'le' | 'ge', min, max, step,
//     presets, seed, format(t), caption?, toggle?, scale?, presetLabel? }
//     seed is the starting cutoff, and statistic the one the cutoff applies to until
//     the shell switches it (options.statistic). format(t) prints a cutoff for labels.
//     caption: a muted line under the slider, a string or caption({ threshold,
//       statistic }) for one that follows the cutoff.
//     toggle: { statistics: ['p', 'q'] } adds a switch that reports the chosen
//       statistic through options.onStatisticChange.
//     scale: 'linear' (default) or 'log'. On a log scale min and max are still
//       cutoffs, step is in decades, and slider position x is the cutoff 10^x.
//     presetLabel(t): a shorter preset pill label than the full cutoff label.

const KPI_TONE_CLASSES = {
  accent: 'kpi-value kpi-text-orange',
  purifying: 'kpi-value kpi-text-blue',
  ink: 'kpi-value',
  muted: 'kpi-value kpi-text-muted'
};

// The tone class of a card's value: fixed, or chosen from the stats when stats are given
function kpiToneClass(card, stats, view) {
  const tone = typeof card.tone === 'function' ? (stats ? card.tone(stats, view) : 'ink') : card.tone;
  return KPI_TONE_CLASSES[tone || 'ink'];
}

function kpiCardMarkup(card, index) {
  const cardClass = card.cardClass ? `kpi-card ${card.cardClass}` : 'kpi-card';
  const widget = card.widget ? `
          <div class="kpi-widget"></div>` : '';
  return `
        <!-- ${index + 1}. ${card.comment || `${card.label} Card`} -->
        <div class="${cardClass}" id="kpi-${card.id}">
          <div class="kpi-label"></div>
          <div class="${kpiToneClass(card)}" id="kpi-val-${card.valueId || card.id}">—</div>
          <div class="kpi-sub" id="kpi-sub-${card.subId || card.id}">${card.placeholder || '—'}</div>${widget}
        </div>
`;
}

// Slider position of a cutoff: the cutoff itself, or its log10 on a log scale
export function sliderPosition(controller, t) {
  return controller.scale === 'log' ? Math.log10(t) : t;
}

// The cutoff at a slider position; a log scale keeps three significant digits
export function sliderThreshold(controller, x) {
  return controller.scale === 'log' ? Number((10 ** x).toPrecision(3)) : x;
}

// A cutoff as an attribute value: the label's form when it reads back as the same
// number (0.20 for 0.2), else the full number (0.999, which the label rounds to 1.00)
export function cutoffAttribute(controller, t) {
  const shown = controller.format(t);
  return parseFloat(shown) === t ? shown : String(t);
}

// The slider's min, max, step and value attributes for a cutoff t
export function sliderAttributes(controller, t) {
  if (controller.scale === 'log') {
    return {
      min: String(Math.log10(controller.min)),
      max: String(Math.log10(controller.max)),
      step: String(controller.step),
      value: String(Math.log10(t))
    };
  }
  return {
    min: cutoffAttribute(controller, controller.min),
    max: cutoffAttribute(controller, controller.max),
    step: cutoffAttribute(controller, controller.step),
    value: cutoffAttribute(controller, t)
  };
}

// Whether preset p is the cutoff t, to within the slider's resolution
export function presetActive(controller, p, t) {
  return Math.abs(sliderPosition(controller, p) - sliderPosition(controller, t)) < 0.002;
}

// Whether a value passes the cutoff in the controller's direction
export function thresholdPasses(controller, value, t) {
  return controller.direction === 'ge' ? value >= t : value <= t;
}

// The direction in words, for captions and sub lines: 'at most' or 'at least'
export function thresholdRelation(controller) {
  return controller.direction === 'ge' ? 'at least' : 'at most';
}

// The caption line's text for a cutoff, or null for a controller without one
export function thresholdCaption(controller, threshold, statistic) {
  const { caption } = controller;
  if (caption == null) return null;
  return typeof caption === 'function' ? caption({ threshold, statistic }) : caption;
}

function statisticToggleMarkup(controller, statistic) {
  if (!controller.toggle) return '';
  const button = (s) => `<button type="button" class="seg-btn${s === statistic ? ' active' : ''}" data-statistic="${s}">${s}</button>`;
  return `
            <span class="segmented-control kpi-stat-toggle">${controller.toggle.statistics.map(button).join('')}</span>`;
}

export function thresholdControllerMarkup(controller, threshold, statistic, index) {
  const pill = (p) => {
    const active = presetActive(controller, p, threshold) ? ' active' : '';
    return `
            <button type="button" class="preset-pill${active}" data-p="${cutoffAttribute(controller, p)}">${presetLabel(controller, p, statistic)}</button>`;
  };
  const slider = sliderAttributes(controller, threshold);
  const caption = controller.caption == null ? '' : `
          <div class="kpi-control-caption" id="kpi-control-caption"></div>`;
  return `
        <!-- ${index + 1}. Interactive Threshold Controller -->
        <div class="kpi-card kpi-card-control" id="kpi-controller">
          <div class="kpi-control-header">
            <span class="kpi-label">${controller.title}</span>${statisticToggleMarkup(controller, statistic)}
            <span class="kpi-thresh-display" id="kpi-thresh-display">${thresholdLabel(controller, threshold, statistic)}</span>
          </div>
          <div class="kpi-slider-row">
            <input type="range" class="kpi-slider" id="kpi-thresh-slider" min="${slider.min}" max="${slider.max}" step="${slider.step}" value="${slider.value}" />
          </div>${caption}
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

// A preset pill's label: the controller's short form, else the full cutoff label
export function presetLabel(controller, t, statistic = controller.statistic) {
  return controller.presetLabel ? controller.presetLabel(t) : thresholdLabel(controller, t, statistic);
}

// Whether the cards can take over the blank cards' markup: the same ids in order
export function cardsFitBlank(cards, blank) {
  return cards.length === blank.length && cards.every((card, i) => card.id === blank[i].id);
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
      onThresholdChange: null,
      onStatisticChange: null
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

    this.controlEl = this.el.querySelector('#kpi-controller');
    this.cardEls = this.bindCards(this.blank);

    // Bind slider & preset listeners
    this.slider = this.el.querySelector('#kpi-thresh-slider');
    this.display = this.el.querySelector('#kpi-thresh-display');
    this.caption = this.el.querySelector('#kpi-control-caption');
    this.presetBtns = this.el.querySelectorAll('.preset-pill');
    this.statisticBtns = this.el.querySelectorAll('.kpi-stat-toggle .seg-btn');
    this.updateCaption();

    this.slider.addEventListener('input', (e) => {
      const val = sliderThreshold(this.controller, parseFloat(e.target.value));
      this.setThreshold(val, true);
    });

    this.presetBtns.forEach(btn => {
      btn.addEventListener('click', () => {
        this.setThreshold(parseFloat(btn.dataset.p), true);
      });
    });

    this.statisticBtns.forEach(btn => {
      btn.addEventListener('click', () => this.setStatistic(btn.dataset.statistic, true));
    });
  }

  // The label, value, sub line and widget slot of each card drawn for list
  bindCards(list) {
    return list.map(card => {
      const root = this.el.querySelector(`#kpi-${card.id}`);
      const label = root.querySelector('.kpi-label');
      label.textContent = card.label;
      return {
        root,
        label,
        value: root.querySelector('.kpi-value'),
        sub: root.querySelector('.kpi-sub'),
        widget: root.querySelector('.kpi-widget')
      };
    });
  }

  // Swap the blank cards for the real ones: relabelled in place when the ids match,
  // else redrawn from the cards' own markup; every card with a widget gets a slot
  fillCards() {
    if (!cardsFitBlank(this.cards, this.blank)) {
      this.cardEls.forEach(els => els.root.remove());
      this.controlEl.insertAdjacentHTML('beforebegin', this.cards.map(kpiCardMarkup).join(''));
      this.blank = this.cards;
      this.cardEls = this.bindCards(this.cards);
    }
    this.cards.forEach((card, i) => {
      const els = this.cardEls[i];
      els.label.textContent = card.label;
      els.value.className = kpiToneClass(card);
      if (card.widget && !els.widget) {
        els.widget = document.createElement('div');
        els.widget.className = 'kpi-widget';
        els.root.appendChild(els.widget);
      }
    });
  }

  updateCaption() {
    if (this.caption) this.caption.textContent = thresholdCaption(this.controller, this.threshold, this.statistic);
  }

  setThreshold(p, triggerCallback = false) {
    this.threshold = p;
    this.slider.value = sliderPosition(this.controller, p);
    this.display.textContent = thresholdLabel(this.controller, p, this.statistic);
    this.updateCaption();

    // Update active preset button highlight
    this.presetBtns.forEach(btn => {
      btn.classList.toggle('active', presetActive(this.controller, parseFloat(btn.dataset.p), p));
    });

    this.updateStats();

    if (triggerCallback && typeof this.options.onThresholdChange === 'function') {
      this.options.onThresholdChange(p);
    }
  }

  // Relabel the cutoff and its presets for the statistic it now applies to
  setStatistic(statistic, triggerCallback = false) {
    this.statistic = statistic;
    this.display.textContent = thresholdLabel(this.controller, this.threshold, statistic);
    this.presetBtns.forEach(btn => {
      btn.textContent = presetLabel(this.controller, parseFloat(btn.dataset.p), statistic);
    });
    this.statisticBtns.forEach(btn => {
      btn.classList.toggle('active', btn.dataset.statistic === statistic);
    });
    this.updateCaption();
    this.updateStats();

    if (triggerCallback && typeof this.options.onStatisticChange === 'function') {
      this.options.onStatisticChange(statistic);
    }
  }

  setData(data) {
    this.data = data;
    this.cards = this.options.cards.filter(card => !card.available || card.available(data));
    this.fillCards();
    this.updateStats();
  }

  updateStats() {
    if (!this.data) return;
    const stats = this.data.getStats(this.threshold);
    const view = { data: this.data, threshold: this.threshold, statistic: this.statistic };

    this.cards.forEach((card, i) => {
      const els = this.cardEls[i];
      els.value.className = kpiToneClass(card, stats, view);
      els.value.textContent = card.value(stats, view);
      els.sub.textContent = card.sub(stats, view);
      if (card.widget) card.widget(els.widget, stats, view);
    });
  }
}
