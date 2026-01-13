import { E as y, U as Tt, T as ye, a as wt, a4 as U, a5 as I, D as L, t as R, a0 as E, M, ae as St, af as Ct, w as X, $, V as J, ag as Pt, ah as N, ai as z, F as Ft, v as Ie, x as We, l as H, d as Le, J as G, a7 as pe, R as ge, I as Ee, c as Y, B as A, n as Te, S as Z, z as ee, aj as Bt, ak as me, N as de, al as j, s as xe, m as He, q as Ve, a9 as Ye, ac as Xe, o as Rt, p as Mt, aa as Ut, ab as kt, ad as zt, am as Gt, an as At, e as C, ao as Dt, __tla as __tla_0 } from "./index-PYdGGKf0.js";
import { c as ie, a as Ot, b as It, B as $e } from "./colorToUniform-BXaCBwVl.js";
import { C as V, __tla as __tla_1 } from "./CanvasPool-C1Qbkdw7.js";
let __tla = Promise.all([
  (() => {
    try {
      return __tla_0;
    } catch {
    }
  })(),
  (() => {
    try {
      return __tla_1;
    } catch {
    }
  })()
]).then(async () => {
  class Ne {
    static init(e) {
      Object.defineProperty(this, "resizeTo", {
        set(t) {
          globalThis.removeEventListener("resize", this.queueResize), this._resizeTo = t, t && (globalThis.addEventListener("resize", this.queueResize), this.resize());
        },
        get() {
          return this._resizeTo;
        }
      }), this.queueResize = () => {
        this._resizeTo && (this._cancelResize(), this._resizeId = requestAnimationFrame(() => this.resize()));
      }, this._cancelResize = () => {
        this._resizeId && (cancelAnimationFrame(this._resizeId), this._resizeId = null);
      }, this.resize = () => {
        if (!this._resizeTo) return;
        this._cancelResize();
        let t, r;
        if (this._resizeTo === globalThis.window) t = globalThis.innerWidth, r = globalThis.innerHeight;
        else {
          const { clientWidth: i, clientHeight: n } = this._resizeTo;
          t = i, r = n;
        }
        this.renderer.resize(t, r), this.render();
      }, this._resizeId = null, this._resizeTo = null, this.resizeTo = e.resizeTo || null;
    }
    static destroy() {
      globalThis.removeEventListener("resize", this.queueResize), this._cancelResize(), this._cancelResize = null, this.queueResize = null, this.resizeTo = null, this.resize = null;
    }
  }
  Ne.extension = y.Application;
  class Ke {
    static init(e) {
      e = Object.assign({
        autoStart: true,
        sharedTicker: false
      }, e), Object.defineProperty(this, "ticker", {
        set(t) {
          this._ticker && this._ticker.remove(this.render, this), this._ticker = t, t && t.add(this.render, this, Tt.LOW);
        },
        get() {
          return this._ticker;
        }
      }), this.stop = () => {
        this._ticker.stop();
      }, this.start = () => {
        this._ticker.start();
      }, this._ticker = null, this.ticker = e.sharedTicker ? ye.shared : new ye(), e.autoStart && this.start();
    }
    static destroy() {
      if (this._ticker) {
        const e = this._ticker;
        this.ticker = null, e.destroy();
      }
    }
  }
  Ke.extension = y.Application;
  class Wt extends wt {
    constructor() {
      super(...arguments), this.chars = /* @__PURE__ */ Object.create(null), this.lineHeight = 0, this.fontFamily = "", this.fontMetrics = {
        fontSize: 0,
        ascent: 0,
        descent: 0
      }, this.baseLineOffset = 0, this.distanceField = {
        type: "none",
        range: 0
      }, this.pages = [], this.applyFillAsTint = true, this.baseMeasurementFontSize = 100, this.baseRenderedFontSize = 100;
    }
    get font() {
      return U(I, "BitmapFont.font is deprecated, please use BitmapFont.fontFamily instead."), this.fontFamily;
    }
    get pageTextures() {
      return U(I, "BitmapFont.pageTextures is deprecated, please use BitmapFont.pages instead."), this.pages;
    }
    get size() {
      return U(I, "BitmapFont.size is deprecated, please use BitmapFont.fontMetrics.fontSize instead."), this.fontMetrics.fontSize;
    }
    get distanceFieldRange() {
      return U(I, "BitmapFont.distanceFieldRange is deprecated, please use BitmapFont.distanceField.range instead."), this.distanceField.range;
    }
    get distanceFieldType() {
      return U(I, "BitmapFont.distanceFieldType is deprecated, please use BitmapFont.distanceField.type instead."), this.distanceField.type;
    }
    destroy(e = false) {
      var _a;
      this.emit("destroy", this), this.removeAllListeners();
      for (const t in this.chars) (_a = this.chars[t].texture) == null ? void 0 : _a.destroy();
      this.chars = null, e && (this.pages.forEach((t) => t.texture.destroy(true)), this.pages = null);
    }
  }
  class Lt {
    constructor(e = 0, t = 0, r = false) {
      this.first = null, this.items = /* @__PURE__ */ Object.create(null), this.last = null, this.max = e, this.resetTtl = r, this.size = 0, this.ttl = t;
    }
    clear() {
      return this.first = null, this.items = /* @__PURE__ */ Object.create(null), this.last = null, this.size = 0, this;
    }
    delete(e) {
      if (this.has(e)) {
        const t = this.items[e];
        delete this.items[e], this.size--, t.prev !== null && (t.prev.next = t.next), t.next !== null && (t.next.prev = t.prev), this.first === t && (this.first = t.next), this.last === t && (this.last = t.prev);
      }
      return this;
    }
    entries(e = this.keys()) {
      return e.map((t) => [
        t,
        this.get(t)
      ]);
    }
    evict(e = false) {
      if (e || this.size > 0) {
        const t = this.first;
        delete this.items[t.key], --this.size === 0 ? (this.first = null, this.last = null) : (this.first = t.next, this.first.prev = null);
      }
      return this;
    }
    expiresAt(e) {
      let t;
      return this.has(e) && (t = this.items[e].expiry), t;
    }
    get(e) {
      const t = this.items[e];
      if (t !== void 0) {
        if (this.ttl > 0 && t.expiry <= Date.now()) {
          this.delete(e);
          return;
        }
        return this.moveToEnd(t), t.value;
      }
    }
    has(e) {
      return e in this.items;
    }
    moveToEnd(e) {
      this.last !== e && (e.prev !== null && (e.prev.next = e.next), e.next !== null && (e.next.prev = e.prev), this.first === e && (this.first = e.next), e.prev = this.last, e.next = null, this.last !== null && (this.last.next = e), this.last = e, this.first === null && (this.first = e));
    }
    keys() {
      const e = [];
      let t = this.first;
      for (; t !== null; ) e.push(t.key), t = t.next;
      return e;
    }
    setWithEvicted(e, t, r = this.resetTtl) {
      let i = null;
      if (this.has(e)) this.set(e, t, true, r);
      else {
        this.max > 0 && this.size === this.max && (i = {
          ...this.first
        }, this.evict(true));
        let n = this.items[e] = {
          expiry: this.ttl > 0 ? Date.now() + this.ttl : this.ttl,
          key: e,
          prev: this.last,
          next: null,
          value: t
        };
        ++this.size === 1 ? this.first = n : this.last.next = n, this.last = n;
      }
      return i;
    }
    set(e, t, r = false, i = this.resetTtl) {
      let n = this.items[e];
      return r || n !== void 0 ? (n.value = t, r === false && i && (n.expiry = this.ttl > 0 ? Date.now() + this.ttl : this.ttl), this.moveToEnd(n)) : (this.max > 0 && this.size === this.max && this.evict(true), n = this.items[e] = {
        expiry: this.ttl > 0 ? Date.now() + this.ttl : this.ttl,
        key: e,
        prev: this.last,
        next: null,
        value: t
      }, ++this.size === 1 ? this.first = n : this.last.next = n, this.last = n), this;
    }
    values(e = this.keys()) {
      return e.map((t) => this.get(t));
    }
  }
  function je(o = 1e3, e = 0, t = false) {
    if (isNaN(o) || o < 0) throw new TypeError("Invalid max value");
    if (isNaN(e) || e < 0) throw new TypeError("Invalid ttl value");
    if (typeof t != "boolean") throw new TypeError("Invalid resetTtl value");
    return new Lt(o, e, t);
  }
  const Et = [
    "serif",
    "sans-serif",
    "monospace",
    "cursive",
    "fantasy",
    "system-ui"
  ];
  function te(o) {
    const e = typeof o.fontSize == "number" ? `${o.fontSize}px` : o.fontSize;
    let t = o.fontFamily;
    Array.isArray(o.fontFamily) || (t = o.fontFamily.split(","));
    for (let r = t.length - 1; r >= 0; r--) {
      let i = t[r].trim();
      !/([\"\'])[^\'\"]+\1/.test(i) && !Et.includes(i) && (i = `"${i}"`), t[r] = i;
    }
    return `${o.fontStyle} ${o.fontVariant} ${o.fontWeight} ${e} ${t.join(",")}`;
  }
  const se = {
    willReadFrequently: true
  }, k = class g {
    static get experimentalLetterSpacingSupported() {
      let e = g._experimentalLetterSpacingSupported;
      if (e === void 0) {
        const t = L.get().getCanvasRenderingContext2D().prototype;
        e = g._experimentalLetterSpacingSupported = "letterSpacing" in t || "textLetterSpacing" in t;
      }
      return e;
    }
    constructor(e, t, r, i, n, s, a, u, l) {
      this.text = e, this.style = t, this.width = r, this.height = i, this.lines = n, this.lineWidths = s, this.lineHeight = a, this.maxLineWidth = u, this.fontProperties = l;
    }
    static measureText(e = " ", t, r = g._canvas, i = t.wordWrap) {
      var _a;
      const n = `${e}-${t.styleKey}-wordWrap-${i}`;
      if (g._measurementCache.has(n)) return g._measurementCache.get(n);
      const s = te(t), a = g.measureFont(s);
      a.fontSize === 0 && (a.fontSize = t.fontSize, a.ascent = t.fontSize);
      const u = g.__context;
      u.font = s;
      const h = (i ? g._wordWrap(e, t, r) : e).split(/(?:\r\n|\r|\n)/), c = new Array(h.length);
      let d = 0;
      for (let _ = 0; _ < h.length; _++) {
        const b = g._measureText(h[_], t.letterSpacing, u);
        c[_] = b, d = Math.max(d, b);
      }
      const f = ((_a = t._stroke) == null ? void 0 : _a.width) || 0;
      let m = d + f;
      t.dropShadow && (m += t.dropShadow.distance);
      const v = t.lineHeight || a.fontSize;
      let x = Math.max(v, a.fontSize + f) + (h.length - 1) * (v + t.leading);
      t.dropShadow && (x += t.dropShadow.distance);
      const p = new g(e, t, m, x, h, c, v + t.leading, d, a);
      return g._measurementCache.set(n, p), p;
    }
    static _measureText(e, t, r) {
      let i = false;
      g.experimentalLetterSpacingSupported && (g.experimentalLetterSpacing ? (r.letterSpacing = `${t}px`, r.textLetterSpacing = `${t}px`, i = true) : (r.letterSpacing = "0px", r.textLetterSpacing = "0px"));
      const n = r.measureText(e);
      let s = n.width;
      const a = -n.actualBoundingBoxLeft;
      let l = n.actualBoundingBoxRight - a;
      if (s > 0) if (i) s -= t, l -= t;
      else {
        const h = (g.graphemeSegmenter(e).length - 1) * t;
        s += h, l += h;
      }
      return Math.max(s, l);
    }
    static _wordWrap(e, t, r = g._canvas) {
      const i = r.getContext("2d", se);
      let n = 0, s = "", a = "";
      const u = /* @__PURE__ */ Object.create(null), { letterSpacing: l, whiteSpace: h } = t, c = g._collapseSpaces(h), d = g._collapseNewlines(h);
      let f = !c;
      const m = t.wordWrapWidth + l, v = g._tokenize(e);
      for (let x = 0; x < v.length; x++) {
        let p = v[x];
        if (g._isNewline(p)) {
          if (!d) {
            a += g._addLine(s), f = !c, s = "", n = 0;
            continue;
          }
          p = " ";
        }
        if (c) {
          const b = g.isBreakingSpace(p), w = g.isBreakingSpace(s[s.length - 1]);
          if (b && w) continue;
        }
        const _ = g._getFromCache(p, l, u, i);
        if (_ > m) if (s !== "" && (a += g._addLine(s), s = "", n = 0), g.canBreakWords(p, t.breakWords)) {
          const b = g.wordWrapSplit(p);
          for (let w = 0; w < b.length; w++) {
            let S = b[w], F = S, P = 1;
            for (; b[w + P]; ) {
              const T = b[w + P];
              if (!g.canBreakChars(F, T, p, w, t.breakWords)) S += T;
              else break;
              F = T, P++;
            }
            w += P - 1;
            const B = g._getFromCache(S, l, u, i);
            B + n > m && (a += g._addLine(s), f = false, s = "", n = 0), s += S, n += B;
          }
        } else {
          s.length > 0 && (a += g._addLine(s), s = "", n = 0);
          const b = x === v.length - 1;
          a += g._addLine(p, !b), f = false, s = "", n = 0;
        }
        else _ + n > m && (f = false, a += g._addLine(s), s = "", n = 0), (s.length > 0 || !g.isBreakingSpace(p) || f) && (s += p, n += _);
      }
      return a += g._addLine(s, false), a;
    }
    static _addLine(e, t = true) {
      return e = g._trimRight(e), e = t ? `${e}
` : e, e;
    }
    static _getFromCache(e, t, r, i) {
      let n = r[e];
      return typeof n != "number" && (n = g._measureText(e, t, i) + t, r[e] = n), n;
    }
    static _collapseSpaces(e) {
      return e === "normal" || e === "pre-line";
    }
    static _collapseNewlines(e) {
      return e === "normal";
    }
    static _trimRight(e) {
      if (typeof e != "string") return "";
      for (let t = e.length - 1; t >= 0; t--) {
        const r = e[t];
        if (!g.isBreakingSpace(r)) break;
        e = e.slice(0, -1);
      }
      return e;
    }
    static _isNewline(e) {
      return typeof e != "string" ? false : g._newlines.includes(e.charCodeAt(0));
    }
    static isBreakingSpace(e, t) {
      return typeof e != "string" ? false : g._breakingSpaces.includes(e.charCodeAt(0));
    }
    static _tokenize(e) {
      const t = [];
      let r = "";
      if (typeof e != "string") return t;
      for (let i = 0; i < e.length; i++) {
        const n = e[i], s = e[i + 1];
        if (g.isBreakingSpace(n, s) || g._isNewline(n)) {
          r !== "" && (t.push(r), r = ""), n === "\r" && s === `
` ? (t.push(`\r
`), i++) : t.push(n);
          continue;
        }
        r += n;
      }
      return r !== "" && t.push(r), t;
    }
    static canBreakWords(e, t) {
      return t;
    }
    static canBreakChars(e, t, r, i, n) {
      return true;
    }
    static wordWrapSplit(e) {
      return g.graphemeSegmenter(e);
    }
    static measureFont(e) {
      if (g._fonts[e]) return g._fonts[e];
      const t = g._context;
      t.font = e;
      const r = t.measureText(g.METRICS_STRING + g.BASELINE_SYMBOL), i = {
        ascent: r.actualBoundingBoxAscent,
        descent: r.actualBoundingBoxDescent,
        fontSize: r.actualBoundingBoxAscent + r.actualBoundingBoxDescent
      };
      return g._fonts[e] = i, i;
    }
    static clearMetrics(e = "") {
      e ? delete g._fonts[e] : g._fonts = {};
    }
    static get _canvas() {
      var _a;
      if (!g.__canvas) {
        let e;
        try {
          const t = new OffscreenCanvas(0, 0);
          if ((_a = t.getContext("2d", se)) == null ? void 0 : _a.measureText) return g.__canvas = t, t;
          e = L.get().createCanvas();
        } catch {
          e = L.get().createCanvas();
        }
        e.width = e.height = 10, g.__canvas = e;
      }
      return g.__canvas;
    }
    static get _context() {
      return g.__context || (g.__context = g._canvas.getContext("2d", se)), g.__context;
    }
  };
  k.METRICS_STRING = "|\xC9q\xC5";
  k.BASELINE_SYMBOL = "M";
  k.BASELINE_MULTIPLIER = 1.4;
  k.HEIGHT_MULTIPLIER = 2;
  k.graphemeSegmenter = (() => {
    if (typeof (Intl == null ? void 0 : Intl.Segmenter) == "function") {
      const o = new Intl.Segmenter();
      return (e) => {
        const t = o.segment(e), r = [];
        let i = 0;
        for (const n of t) r[i++] = n.segment;
        return r;
      };
    }
    return (o) => [
      ...o
    ];
  })();
  k.experimentalLetterSpacing = false;
  k._fonts = {};
  k._newlines = [
    10,
    13
  ];
  k._breakingSpaces = [
    9,
    32,
    8192,
    8193,
    8194,
    8195,
    8196,
    8197,
    8198,
    8200,
    8201,
    8202,
    8287,
    12288
  ];
  k._measurementCache = je(1e3);
  let D = k;
  const we = 1e5;
  function re(o, e, t, r = 0) {
    if (o.texture === R.WHITE && !o.fill) return E.shared.setValue(o.color).setAlpha(o.alpha ?? 1).toHexa();
    if (o.fill) {
      if (o.fill instanceof St) {
        const i = o.fill, n = e.createPattern(i.texture.source.resource, "repeat"), s = i.transform.copyTo(M.shared);
        return s.scale(i.texture.frame.width, i.texture.frame.height), n.setTransform(s), n;
      } else if (o.fill instanceof Ct) {
        const i = o.fill, n = i.type === "linear", s = i.textureSpace === "local";
        let a = 1, u = 1;
        s && t && (a = t.width + r, u = t.height + r);
        let l, h = false;
        if (n) {
          const { start: c, end: d } = i;
          l = e.createLinearGradient(c.x * a, c.y * u, d.x * a, d.y * u), h = Math.abs(d.x - c.x) < Math.abs((d.y - c.y) * 0.1);
        } else {
          const { center: c, innerRadius: d, outerCenter: f, outerRadius: m } = i;
          l = e.createRadialGradient(c.x * a, c.y * u, d * a, f.x * a, f.y * u, m * a);
        }
        if (h && s && t) {
          const c = t.lineHeight / u;
          for (let d = 0; d < t.lines.length; d++) {
            const f = (d * t.lineHeight + r / 2) / u;
            i.colorStops.forEach((m) => {
              const v = f + m.offset * c;
              l.addColorStop(Math.floor(v * we) / we, E.shared.setValue(m.color).toHex());
            });
          }
        } else i.colorStops.forEach((c) => {
          l.addColorStop(c.offset, E.shared.setValue(c.color).toHex());
        });
        return l;
      }
    } else {
      const i = e.createPattern(o.texture.source.resource, "repeat"), n = o.matrix.copyTo(M.shared);
      return n.scale(o.texture.frame.width, o.texture.frame.height), i.setTransform(n), i;
    }
    return X("FillStyle not recognised", o), "red";
  }
  const qe = class Qe extends Wt {
    constructor(e) {
      super(), this.resolution = 1, this.pages = [], this._padding = 0, this._measureCache = /* @__PURE__ */ Object.create(null), this._currentChars = [], this._currentX = 0, this._currentY = 0, this._currentMaxCharHeight = 0, this._currentPageIndex = -1, this._skipKerning = false;
      const t = {
        ...Qe.defaultOptions,
        ...e
      };
      this._textureSize = t.textureSize, this._mipmap = t.mipmap;
      const r = t.style.clone();
      t.overrideFill && (r._fill.color = 16777215, r._fill.alpha = 1, r._fill.texture = R.WHITE, r._fill.fill = null), this.applyFillAsTint = t.overrideFill;
      const i = r.fontSize;
      r.fontSize = this.baseMeasurementFontSize;
      const n = te(r);
      t.overrideSize ? r._stroke && (r._stroke.width *= this.baseRenderedFontSize / i) : r.fontSize = this.baseRenderedFontSize = i, this._style = r, this._skipKerning = t.skipKerning ?? false, this.resolution = t.resolution ?? 1, this._padding = t.padding ?? 4, t.textureStyle && (this._textureStyle = t.textureStyle instanceof $ ? t.textureStyle : new $(t.textureStyle)), this.fontMetrics = D.measureFont(n), this.lineHeight = r.lineHeight || this.fontMetrics.fontSize || r.fontSize;
    }
    ensureCharacters(e) {
      var _a, _b;
      const t = D.graphemeSegmenter(e).filter((x) => !this._currentChars.includes(x)).filter((x, p, _) => _.indexOf(x) === p);
      if (!t.length) return;
      this._currentChars = [
        ...this._currentChars,
        ...t
      ];
      let r;
      this._currentPageIndex === -1 ? r = this._nextPage() : r = this.pages[this._currentPageIndex];
      let { canvas: i, context: n } = r.canvasAndContext, s = r.texture.source;
      const a = this._style;
      let u = this._currentX, l = this._currentY, h = this._currentMaxCharHeight;
      const c = this.baseRenderedFontSize / this.baseMeasurementFontSize, d = this._padding * c;
      let f = false;
      const m = i.width / this.resolution, v = i.height / this.resolution;
      for (let x = 0; x < t.length; x++) {
        const p = t[x], _ = D.measureText(p, a, i, false);
        _.lineHeight = _.height;
        const b = _.width * c, w = Math.ceil((a.fontStyle === "italic" ? 2 : 1) * b), S = _.height * c, F = w + d * 2, P = S + d * 2;
        if (f = false, p !== `
` && p !== "\r" && p !== "	" && p !== " " && (f = true, h = Math.ceil(Math.max(P, h))), u + F > m && (l += h, h = P, u = 0, l + h > v)) {
          s.update();
          const T = this._nextPage();
          i = T.canvasAndContext.canvas, n = T.canvasAndContext.context, s = T.texture.source, u = 0, l = 0, h = 0;
        }
        const B = b / c - (((_a = a.dropShadow) == null ? void 0 : _a.distance) ?? 0) - (((_b = a._stroke) == null ? void 0 : _b.width) ?? 0);
        if (this.chars[p] = {
          id: p.codePointAt(0),
          xOffset: -this._padding,
          yOffset: -this._padding,
          xAdvance: B,
          kerning: {}
        }, f) {
          this._drawGlyph(n, _, u + d, l + d, c, a);
          const T = s.width * c, K = s.height * c, ne = new J(u / T * s.width, l / K * s.height, F / T * s.width, P / K * s.height);
          this.chars[p].texture = new R({
            source: s,
            frame: ne
          }), u += Math.ceil(F);
        }
      }
      s.update(), this._currentX = u, this._currentY = l, this._currentMaxCharHeight = h, this._skipKerning && this._applyKerning(t, n);
    }
    get pageTextures() {
      return U(I, "BitmapFont.pageTextures is deprecated, please use BitmapFont.pages instead."), this.pages;
    }
    _applyKerning(e, t) {
      const r = this._measureCache;
      for (let i = 0; i < e.length; i++) {
        const n = e[i];
        for (let s = 0; s < this._currentChars.length; s++) {
          const a = this._currentChars[s];
          let u = r[n];
          u || (u = r[n] = t.measureText(n).width);
          let l = r[a];
          l || (l = r[a] = t.measureText(a).width);
          let h = t.measureText(n + a).width, c = h - (u + l);
          c && (this.chars[n].kerning[a] = c), h = t.measureText(n + a).width, c = h - (u + l), c && (this.chars[a].kerning[n] = c);
        }
      }
    }
    _nextPage() {
      this._currentPageIndex++;
      const e = this.resolution, t = V.getOptimalCanvasAndContext(this._textureSize, this._textureSize, e);
      this._setupContext(t.context, this._style, e);
      const r = e * (this.baseRenderedFontSize / this.baseMeasurementFontSize), i = new R({
        source: new Pt({
          resource: t.canvas,
          resolution: r,
          alphaMode: "premultiply-alpha-on-upload",
          autoGenerateMipmaps: this._mipmap
        })
      });
      this._textureStyle && (i.source.style = this._textureStyle);
      const n = {
        canvasAndContext: t,
        texture: i
      };
      return this.pages[this._currentPageIndex] = n, n;
    }
    _setupContext(e, t, r) {
      t.fontSize = this.baseRenderedFontSize, e.scale(r, r), e.font = te(t), t.fontSize = this.baseMeasurementFontSize, e.textBaseline = t.textBaseline;
      const i = t._stroke, n = (i == null ? void 0 : i.width) ?? 0;
      if (i && (e.lineWidth = n, e.lineJoin = i.join, e.miterLimit = i.miterLimit, e.strokeStyle = re(i, e)), t._fill && (e.fillStyle = re(t._fill, e)), t.dropShadow) {
        const s = t.dropShadow, a = E.shared.setValue(s.color).toArray(), u = s.blur * r, l = s.distance * r;
        e.shadowColor = `rgba(${a[0] * 255},${a[1] * 255},${a[2] * 255},${s.alpha})`, e.shadowBlur = u, e.shadowOffsetX = Math.cos(s.angle) * l, e.shadowOffsetY = Math.sin(s.angle) * l;
      } else e.shadowColor = "black", e.shadowBlur = 0, e.shadowOffsetX = 0, e.shadowOffsetY = 0;
    }
    _drawGlyph(e, t, r, i, n, s) {
      var _a;
      const a = t.text, u = t.fontProperties, h = (((_a = s._stroke) == null ? void 0 : _a.width) ?? 0) * n, c = r + h / 2, d = i - h / 2, f = u.descent * n, m = t.lineHeight * n;
      let v = false;
      s.stroke && h && (v = true, e.strokeText(a, c, d + m - f));
      const { shadowBlur: x, shadowOffsetX: p, shadowOffsetY: _ } = e;
      s._fill && (v && (e.shadowBlur = 0, e.shadowOffsetX = 0, e.shadowOffsetY = 0), e.fillText(a, c, d + m - f)), v && (e.shadowBlur = x, e.shadowOffsetX = p, e.shadowOffsetY = _);
    }
    destroy() {
      super.destroy();
      for (let e = 0; e < this.pages.length; e++) {
        const { canvasAndContext: t, texture: r } = this.pages[e];
        V.returnCanvasAndContext(t), r.destroy(true);
      }
      this.pages = null;
    }
  };
  qe.defaultOptions = {
    textureSize: 512,
    style: new N(),
    mipmap: true
  };
  let Se = qe;
  function Je(o, e, t, r) {
    const i = {
      width: 0,
      height: 0,
      offsetY: 0,
      scale: e.fontSize / t.baseMeasurementFontSize,
      lines: [
        {
          width: 0,
          charPositions: [],
          spaceWidth: 0,
          spacesIndex: [],
          chars: []
        }
      ]
    };
    i.offsetY = t.baseLineOffset;
    let n = i.lines[0], s = null, a = true;
    const u = {
      width: 0,
      start: 0,
      index: 0,
      positions: [],
      chars: []
    }, l = t.baseMeasurementFontSize / e.fontSize, h = e.letterSpacing * l, c = e.wordWrapWidth * l, d = e.lineHeight ? e.lineHeight * l : t.lineHeight, f = e.wordWrap && e.breakWords, m = (p) => {
      const _ = n.width;
      for (let b = 0; b < u.index; b++) {
        const w = p.positions[b];
        n.chars.push(p.chars[b]), n.charPositions.push(w + _);
      }
      n.width += p.width, a = false, u.width = 0, u.index = 0, u.chars.length = 0;
    }, v = () => {
      let p = n.chars.length - 1;
      if (r) {
        let _ = n.chars[p];
        for (; _ === " "; ) n.width -= t.chars[_].xAdvance, _ = n.chars[--p];
      }
      i.width = Math.max(i.width, n.width), n = {
        width: 0,
        charPositions: [],
        chars: [],
        spaceWidth: 0,
        spacesIndex: []
      }, a = true, i.lines.push(n), i.height += d;
    }, x = (p) => p - h > c;
    for (let p = 0; p < o.length + 1; p++) {
      let _;
      const b = p === o.length;
      b || (_ = o[p]);
      const w = t.chars[_] || t.chars[" "];
      if (/(?:\s)/.test(_) || _ === "\r" || _ === `
` || b) {
        if (!a && e.wordWrap && x(n.width + u.width) ? (v(), m(u), b || n.charPositions.push(0)) : (u.start = n.width, m(u), b || n.charPositions.push(0)), _ === "\r" || _ === `
`) v();
        else if (!b) {
          const B = w.xAdvance + (w.kerning[s] || 0) + h;
          n.width += B, n.spaceWidth = B, n.spacesIndex.push(n.charPositions.length), n.chars.push(_);
        }
      } else {
        const P = w.kerning[s] || 0, B = w.xAdvance + P + h;
        f && x(n.width + u.width + B) && (m(u), v()), u.positions[u.index++] = u.width + P, u.chars.push(_), u.width += B;
      }
      s = _;
    }
    return v(), e.align === "center" ? Ht(i) : e.align === "right" ? Vt(i) : e.align === "justify" && Yt(i), i;
  }
  function Ht(o) {
    for (let e = 0; e < o.lines.length; e++) {
      const t = o.lines[e], r = o.width / 2 - t.width / 2;
      for (let i = 0; i < t.charPositions.length; i++) t.charPositions[i] += r;
    }
  }
  function Vt(o) {
    for (let e = 0; e < o.lines.length; e++) {
      const t = o.lines[e], r = o.width - t.width;
      for (let i = 0; i < t.charPositions.length; i++) t.charPositions[i] += r;
    }
  }
  function Yt(o) {
    const e = o.width;
    for (let t = 0; t < o.lines.length; t++) {
      const r = o.lines[t];
      let i = 0, n = r.spacesIndex[i++], s = 0;
      const a = r.spacesIndex.length, l = (e - r.width) / a;
      for (let h = 0; h < r.charPositions.length; h++) h === n && (n = r.spacesIndex[i++], s += l), r.charPositions[h] += s;
    }
  }
  function Xt(o) {
    if (o === "") return [];
    typeof o == "string" && (o = [
      o
    ]);
    const e = [];
    for (let t = 0, r = o.length; t < r; t++) {
      const i = o[t];
      if (Array.isArray(i)) {
        if (i.length !== 2) throw new Error(`[BitmapFont]: Invalid character range length, expecting 2 got ${i.length}.`);
        if (i[0].length === 0 || i[1].length === 0) throw new Error("[BitmapFont]: Invalid character delimiter.");
        const n = i[0].charCodeAt(0), s = i[1].charCodeAt(0);
        if (s < n) throw new Error("[BitmapFont]: Invalid character range.");
        for (let a = n, u = s; a <= u; a++) e.push(String.fromCharCode(a));
      } else e.push(...Array.from(i));
    }
    if (e.length === 0) throw new Error("[BitmapFont]: Empty set when resolving characters.");
    return e;
  }
  let q = 0;
  class $t {
    constructor() {
      this.ALPHA = [
        [
          "a",
          "z"
        ],
        [
          "A",
          "Z"
        ],
        " "
      ], this.NUMERIC = [
        [
          "0",
          "9"
        ]
      ], this.ALPHANUMERIC = [
        [
          "a",
          "z"
        ],
        [
          "A",
          "Z"
        ],
        [
          "0",
          "9"
        ],
        " "
      ], this.ASCII = [
        [
          " ",
          "~"
        ]
      ], this.defaultOptions = {
        chars: this.ALPHANUMERIC,
        resolution: 1,
        padding: 4,
        skipKerning: false,
        textureStyle: null
      }, this.measureCache = je(1e3);
    }
    getFont(e, t) {
      var _a;
      let r = `${t.fontFamily}-bitmap`, i = true;
      if (t._fill.fill && !t._stroke ? (r += t._fill.fill.styleKey, i = false) : (t._stroke || t.dropShadow) && (r = `${t.styleKey}-bitmap`, i = false), !z.has(r)) {
        const s = Object.create(t);
        s.lineHeight = 0;
        const a = new Se({
          style: s,
          overrideFill: i,
          overrideSize: true,
          ...this.defaultOptions
        });
        q++, q > 50 && X("BitmapText", `You have dynamically created ${q} bitmap fonts, this can be inefficient. Try pre installing your font styles using \`BitmapFont.install({name:"style1", style})\``), a.once("destroy", () => {
          q--, z.remove(r);
        }), z.set(r, a);
      }
      const n = z.get(r);
      return (_a = n.ensureCharacters) == null ? void 0 : _a.call(n, e), n;
    }
    getLayout(e, t, r = true) {
      const i = this.getFont(e, t), n = `${e}-${t.styleKey}-${r}`;
      if (this.measureCache.has(n)) return this.measureCache.get(n);
      const s = D.graphemeSegmenter(e), a = Je(s, t, i, r);
      return this.measureCache.set(n, a), a;
    }
    measureText(e, t, r = true) {
      return this.getLayout(e, t, r);
    }
    install(...e) {
      var _a, _b, _c, _d;
      let t = e[0];
      typeof t == "string" && (t = {
        name: t,
        style: e[1],
        chars: (_a = e[2]) == null ? void 0 : _a.chars,
        resolution: (_b = e[2]) == null ? void 0 : _b.resolution,
        padding: (_c = e[2]) == null ? void 0 : _c.padding,
        skipKerning: (_d = e[2]) == null ? void 0 : _d.skipKerning
      }, U(I, "BitmapFontManager.install(name, style, options) is deprecated, use BitmapFontManager.install({name, style, ...options})"));
      const r = t == null ? void 0 : t.name;
      if (!r) throw new Error("[BitmapFontManager] Property `name` is required.");
      t = {
        ...this.defaultOptions,
        ...t
      };
      const i = t.style, n = i instanceof N ? i : new N(i), s = t.dynamicFill ?? this._canUseTintForStyle(n), a = new Se({
        style: n,
        overrideFill: s,
        skipKerning: t.skipKerning,
        padding: t.padding,
        resolution: t.resolution,
        overrideSize: false,
        textureStyle: t.textureStyle
      }), u = Xt(t.chars);
      return a.ensureCharacters(u.join("")), z.set(`${r}-bitmap`, a), a.once("destroy", () => z.remove(`${r}-bitmap`)), a;
    }
    uninstall(e) {
      const t = `${e}-bitmap`, r = z.get(t);
      r && r.destroy();
    }
    _canUseTintForStyle(e) {
      return !e._stroke && (!e.dropShadow || e.dropShadow.color === 0) && !e._fill.fill && e._fill.color === 16777215;
    }
  }
  const Nt = new $t();
  var Kt = `in vec2 aPosition;
out vec2 vTextureCoord;

uniform vec4 uInputSize;
uniform vec4 uOutputFrame;
uniform vec4 uOutputTexture;

vec4 filterVertexPosition( void )
{
    vec2 position = aPosition * uOutputFrame.zw + uOutputFrame.xy;
    
    position.x = position.x * (2.0 / uOutputTexture.x) - 1.0;
    position.y = position.y * (2.0*uOutputTexture.z / uOutputTexture.y) - uOutputTexture.z;

    return vec4(position, 0.0, 1.0);
}

vec2 filterTextureCoord( void )
{
    return aPosition * (uOutputFrame.zw * uInputSize.zw);
}

void main(void)
{
    gl_Position = filterVertexPosition();
    vTextureCoord = filterTextureCoord();
}
`, jt = `in vec2 vTextureCoord;
out vec4 finalColor;
uniform sampler2D uTexture;
void main() {
    finalColor = texture(uTexture, vTextureCoord);
}
`, Ce = `struct GlobalFilterUniforms {
  uInputSize: vec4<f32>,
  uInputPixel: vec4<f32>,
  uInputClamp: vec4<f32>,
  uOutputFrame: vec4<f32>,
  uGlobalFrame: vec4<f32>,
  uOutputTexture: vec4<f32>,
};

@group(0) @binding(0) var <uniform> gfu: GlobalFilterUniforms;
@group(0) @binding(1) var uTexture: texture_2d<f32>;
@group(0) @binding(2) var uSampler: sampler;

struct VSOutput {
  @builtin(position) position: vec4<f32>,
  @location(0) uv: vec2<f32>
};

fn filterVertexPosition(aPosition: vec2<f32>) -> vec4<f32>
{
    var position = aPosition * gfu.uOutputFrame.zw + gfu.uOutputFrame.xy;

    position.x = position.x * (2.0 / gfu.uOutputTexture.x) - 1.0;
    position.y = position.y * (2.0 * gfu.uOutputTexture.z / gfu.uOutputTexture.y) - gfu.uOutputTexture.z;

    return vec4(position, 0.0, 1.0);
}

fn filterTextureCoord(aPosition: vec2<f32>) -> vec2<f32>
{
    return aPosition * (gfu.uOutputFrame.zw * gfu.uInputSize.zw);
}

@vertex
fn mainVertex(
  @location(0) aPosition: vec2<f32>,
) -> VSOutput {
  return VSOutput(
   filterVertexPosition(aPosition),
   filterTextureCoord(aPosition)
  );
}

@fragment
fn mainFragment(
  @location(0) uv: vec2<f32>,
) -> @location(0) vec4<f32> {
    return textureSample(uTexture, uSampler, uv);
}
`;
  class qt extends Ft {
    constructor() {
      const e = Ie.from({
        vertex: {
          source: Ce,
          entryPoint: "mainVertex"
        },
        fragment: {
          source: Ce,
          entryPoint: "mainFragment"
        },
        name: "passthrough-filter"
      }), t = We.from({
        vertex: Kt,
        fragment: jt,
        name: "passthrough-filter"
      });
      super({
        gpuProgram: e,
        glProgram: t
      });
    }
  }
  class Ze {
    constructor(e) {
      this._renderer = e;
    }
    push(e, t, r) {
      this._renderer.renderPipes.batch.break(r), r.add({
        renderPipeId: "filter",
        canBundle: false,
        action: "pushFilter",
        container: t,
        filterEffect: e
      });
    }
    pop(e, t, r) {
      this._renderer.renderPipes.batch.break(r), r.add({
        renderPipeId: "filter",
        action: "popFilter",
        canBundle: false
      });
    }
    execute(e) {
      e.action === "pushFilter" ? this._renderer.filter.push(e) : e.action === "popFilter" && this._renderer.filter.pop();
    }
    destroy() {
      this._renderer = null;
    }
  }
  Ze.extension = {
    type: [
      y.WebGLPipes,
      y.WebGPUPipes,
      y.CanvasPipes
    ],
    name: "filter"
  };
  const Pe = new M();
  function Qt(o, e) {
    e.clear();
    const t = e.matrix;
    for (let r = 0; r < o.length; r++) {
      const i = o[r];
      if (i.globalDisplayStatus < 7) continue;
      const n = i.renderGroup ?? i.parentRenderGroup;
      (n == null ? void 0 : n.isCachedAsTexture) ? e.matrix = Pe.copyFrom(n.textureOffsetInverseTransform).append(i.worldTransform) : (n == null ? void 0 : n._parentCacheAsTextureRenderGroup) ? e.matrix = Pe.copyFrom(n._parentCacheAsTextureRenderGroup.inverseWorldTransform).append(i.groupTransform) : e.matrix = i.worldTransform, e.addBounds(i.bounds);
    }
    return e.matrix = t, e;
  }
  const Jt = new pe({
    attributes: {
      aPosition: {
        buffer: new Float32Array([
          0,
          0,
          1,
          0,
          1,
          1,
          0,
          1
        ]),
        format: "float32x2",
        stride: 8,
        offset: 0
      }
    },
    indexBuffer: new Uint32Array([
      0,
      1,
      2,
      0,
      2,
      3
    ])
  });
  class Zt {
    constructor() {
      this.skip = false, this.inputTexture = null, this.backTexture = null, this.filters = null, this.bounds = new Ee(), this.container = null, this.blendRequired = false, this.outputRenderSurface = null, this.globalFrame = {
        x: 0,
        y: 0,
        width: 0,
        height: 0
      }, this.firstEnabledIndex = -1, this.lastEnabledIndex = -1;
    }
  }
  class et {
    constructor(e) {
      this._filterStackIndex = 0, this._filterStack = [], this._filterGlobalUniforms = new H({
        uInputSize: {
          value: new Float32Array(4),
          type: "vec4<f32>"
        },
        uInputPixel: {
          value: new Float32Array(4),
          type: "vec4<f32>"
        },
        uInputClamp: {
          value: new Float32Array(4),
          type: "vec4<f32>"
        },
        uOutputFrame: {
          value: new Float32Array(4),
          type: "vec4<f32>"
        },
        uGlobalFrame: {
          value: new Float32Array(4),
          type: "vec4<f32>"
        },
        uOutputTexture: {
          value: new Float32Array(4),
          type: "vec4<f32>"
        }
      }), this._globalFilterBindGroup = new Le({}), this.renderer = e;
    }
    get activeBackTexture() {
      var _a;
      return (_a = this._activeFilterData) == null ? void 0 : _a.backTexture;
    }
    push(e) {
      const t = this.renderer, r = e.filterEffect.filters, i = this._pushFilterData();
      i.skip = false, i.filters = r, i.container = e.container, i.outputRenderSurface = t.renderTarget.renderSurface;
      const n = t.renderTarget.renderTarget.colorTexture.source, s = n.resolution, a = n.antialias;
      if (r.every((f) => !f.enabled)) {
        i.skip = true;
        return;
      }
      const u = i.bounds;
      if (this._calculateFilterArea(e, u), this._calculateFilterBounds(i, t.renderTarget.rootViewPort, a, s, 1), i.skip) return;
      const l = this._getPreviousFilterData(), h = this._findFilterResolution(s);
      let c = 0, d = 0;
      l && (c = l.bounds.minX, d = l.bounds.minY), this._calculateGlobalFrame(i, c, d, h, n.width, n.height), this._setupFilterTextures(i, u, t, l);
    }
    generateFilteredTexture({ texture: e, filters: t }) {
      const r = this._pushFilterData();
      this._activeFilterData = r, r.skip = false, r.filters = t;
      const i = e.source, n = i.resolution, s = i.antialias;
      if (t.every((f) => !f.enabled)) return r.skip = true, e;
      const a = r.bounds;
      if (a.addRect(e.frame), this._calculateFilterBounds(r, a.rectangle, s, n, 0), r.skip) return e;
      const u = n;
      this._calculateGlobalFrame(r, 0, 0, u, i.width, i.height), r.outputRenderSurface = G.getOptimalTexture(a.width, a.height, r.resolution, r.antialias), r.backTexture = R.EMPTY, r.inputTexture = e, this.renderer.renderTarget.finishRenderPass(), this._applyFiltersToTexture(r, true);
      const d = r.outputRenderSurface;
      return d.source.alphaMode = "premultiplied-alpha", d;
    }
    pop() {
      const e = this.renderer, t = this._popFilterData();
      t.skip || (e.globalUniforms.pop(), e.renderTarget.finishRenderPass(), this._activeFilterData = t, this._applyFiltersToTexture(t, false), t.blendRequired && G.returnTexture(t.backTexture), G.returnTexture(t.inputTexture));
    }
    getBackTexture(e, t, r) {
      const i = e.colorTexture.source._resolution, n = G.getOptimalTexture(t.width, t.height, i, false);
      let s = t.minX, a = t.minY;
      r && (s -= r.minX, a -= r.minY), s = Math.floor(s * i), a = Math.floor(a * i);
      const u = Math.ceil(t.width * i), l = Math.ceil(t.height * i);
      return this.renderer.renderTarget.copyToTexture(e, n, {
        x: s,
        y: a
      }, {
        width: u,
        height: l
      }, {
        x: 0,
        y: 0
      }), n;
    }
    applyFilter(e, t, r, i) {
      const n = this.renderer, s = this._activeFilterData, u = s.outputRenderSurface === r, l = n.renderTarget.rootRenderTarget.colorTexture.source._resolution, h = this._findFilterResolution(l);
      let c = 0, d = 0;
      if (u) {
        const m = this._findPreviousFilterOffset();
        c = m.x, d = m.y;
      }
      this._updateFilterUniforms(t, r, s, c, d, h, u, i);
      const f = e.enabled ? e : this._getPassthroughFilter();
      this._setupBindGroupsAndRender(f, t, n);
    }
    calculateSpriteMatrix(e, t) {
      const r = this._activeFilterData, i = e.set(r.inputTexture._source.width, 0, 0, r.inputTexture._source.height, r.bounds.minX, r.bounds.minY), n = t.worldTransform.copyTo(M.shared), s = t.renderGroup || t.parentRenderGroup;
      return s && s.cacheToLocalTransform && n.prepend(s.cacheToLocalTransform), n.invert(), i.prepend(n), i.scale(1 / t.texture.orig.width, 1 / t.texture.orig.height), i.translate(t.anchor.x, t.anchor.y), i;
    }
    destroy() {
      var _a;
      (_a = this._passthroughFilter) == null ? void 0 : _a.destroy(true), this._passthroughFilter = null;
    }
    _getPassthroughFilter() {
      return this._passthroughFilter ?? (this._passthroughFilter = new qt()), this._passthroughFilter;
    }
    _setupBindGroupsAndRender(e, t, r) {
      if (r.renderPipes.uniformBatch) {
        const i = r.renderPipes.uniformBatch.getUboResource(this._filterGlobalUniforms);
        this._globalFilterBindGroup.setResource(i, 0);
      } else this._globalFilterBindGroup.setResource(this._filterGlobalUniforms, 0);
      this._globalFilterBindGroup.setResource(t.source, 1), this._globalFilterBindGroup.setResource(t.source.style, 2), e.groups[0] = this._globalFilterBindGroup, r.encoder.draw({
        geometry: Jt,
        shader: e,
        state: e._state,
        topology: "triangle-list"
      }), r.type === ge.WEBGL && r.renderTarget.finishRenderPass();
    }
    _setupFilterTextures(e, t, r, i) {
      if (e.backTexture = R.EMPTY, e.inputTexture = G.getOptimalTexture(t.width, t.height, e.resolution, e.antialias), e.blendRequired) {
        r.renderTarget.finishRenderPass();
        const n = r.renderTarget.getRenderTarget(e.outputRenderSurface);
        e.backTexture = this.getBackTexture(n, t, i == null ? void 0 : i.bounds);
      }
      r.renderTarget.bind(e.inputTexture, true), r.globalUniforms.push({
        offset: t
      });
    }
    _calculateGlobalFrame(e, t, r, i, n, s) {
      const a = e.globalFrame;
      a.x = t * i, a.y = r * i, a.width = n * i, a.height = s * i;
    }
    _updateFilterUniforms(e, t, r, i, n, s, a, u) {
      const l = this._filterGlobalUniforms.uniforms, h = l.uOutputFrame, c = l.uInputSize, d = l.uInputPixel, f = l.uInputClamp, m = l.uGlobalFrame, v = l.uOutputTexture;
      a ? (h[0] = r.bounds.minX - i, h[1] = r.bounds.minY - n) : (h[0] = 0, h[1] = 0), h[2] = e.frame.width, h[3] = e.frame.height, c[0] = e.source.width, c[1] = e.source.height, c[2] = 1 / c[0], c[3] = 1 / c[1], d[0] = e.source.pixelWidth, d[1] = e.source.pixelHeight, d[2] = 1 / d[0], d[3] = 1 / d[1], f[0] = 0.5 * d[2], f[1] = 0.5 * d[3], f[2] = e.frame.width * c[2] - 0.5 * d[2], f[3] = e.frame.height * c[3] - 0.5 * d[3];
      const x = this.renderer.renderTarget.rootRenderTarget.colorTexture;
      m[0] = i * s, m[1] = n * s, m[2] = x.source.width * s, m[3] = x.source.height * s, t instanceof R && (t.source.resource = null);
      const p = this.renderer.renderTarget.getRenderTarget(t);
      this.renderer.renderTarget.bind(t, !!u), t instanceof R ? (v[0] = t.frame.width, v[1] = t.frame.height) : (v[0] = p.width, v[1] = p.height), v[2] = p.isRoot ? -1 : 1, this._filterGlobalUniforms.update();
    }
    _findFilterResolution(e) {
      let t = this._filterStackIndex - 1;
      for (; t > 0 && this._filterStack[t].skip; ) --t;
      return t > 0 && this._filterStack[t].inputTexture ? this._filterStack[t].inputTexture.source._resolution : e;
    }
    _findPreviousFilterOffset() {
      let e = 0, t = 0, r = this._filterStackIndex;
      for (; r > 0; ) {
        r--;
        const i = this._filterStack[r];
        if (!i.skip) {
          e = i.bounds.minX, t = i.bounds.minY;
          break;
        }
      }
      return {
        x: e,
        y: t
      };
    }
    _calculateFilterArea(e, t) {
      if (e.renderables ? Qt(e.renderables, t) : e.filterEffect.filterArea ? (t.clear(), t.addRect(e.filterEffect.filterArea), t.applyMatrix(e.container.worldTransform)) : e.container.getFastGlobalBounds(true, t), e.container) {
        const i = (e.container.renderGroup || e.container.parentRenderGroup).cacheToLocalTransform;
        i && t.applyMatrix(i);
      }
    }
    _applyFiltersToTexture(e, t) {
      const r = e.inputTexture, i = e.bounds, n = e.filters, s = e.firstEnabledIndex, a = e.lastEnabledIndex;
      if (this._globalFilterBindGroup.setResource(r.source.style, 2), this._globalFilterBindGroup.setResource(e.backTexture.source, 3), s === a) n[s].apply(this, r, e.outputRenderSurface, t);
      else {
        let u = e.inputTexture;
        const l = G.getOptimalTexture(i.width, i.height, u.source._resolution, false);
        let h = l;
        for (let c = s; c < a; c++) {
          const d = n[c];
          if (!d.enabled) continue;
          d.apply(this, u, h, true);
          const f = u;
          u = h, h = f;
        }
        n[a].apply(this, u, e.outputRenderSurface, t), G.returnTexture(l);
      }
    }
    _calculateFilterBounds(e, t, r, i, n) {
      var _a;
      const s = this.renderer, a = e.bounds, u = e.filters;
      let l = 1 / 0, h = 0, c = true, d = false, f = false, m = true, v = -1, x = -1;
      for (let p = 0; p < u.length; p++) {
        const _ = u[p];
        if (!_.enabled) continue;
        if (v === -1 && (v = p), x = p, l = Math.min(l, _.resolution === "inherit" ? i : _.resolution), h += _.padding, _.antialias === "off" ? c = false : _.antialias === "inherit" && c && (c = r), _.clipToViewport || (m = false), !!!(_.compatibleRenderers & s.type)) {
          f = false;
          break;
        }
        if (_.blendRequired && !(((_a = s.backBuffer) == null ? void 0 : _a.useBackBuffer) ?? true)) {
          X("Blend filter requires backBuffer on WebGL renderer to be enabled. Set `useBackBuffer: true` in the renderer options."), f = false;
          break;
        }
        f = true, d || (d = _.blendRequired);
      }
      if (!f) {
        e.skip = true;
        return;
      }
      if (m && a.fitBounds(0, t.width / i, 0, t.height / i), a.scale(l).ceil().scale(1 / l).pad((h | 0) * n), !a.isPositive) {
        e.skip = true;
        return;
      }
      e.antialias = c, e.resolution = l, e.blendRequired = d, e.firstEnabledIndex = v, e.lastEnabledIndex = x;
    }
    _popFilterData() {
      return this._filterStackIndex--, this._filterStack[this._filterStackIndex];
    }
    _getPreviousFilterData() {
      let e, t = this._filterStackIndex - 1;
      for (; t > 0 && (t--, e = this._filterStack[t], !!e.skip); ) ;
      return e;
    }
    _pushFilterData() {
      let e = this._filterStack[this._filterStackIndex];
      return e || (e = this._filterStack[this._filterStackIndex] = new Zt()), this._filterStackIndex++, e;
    }
  }
  et.extension = {
    type: [
      y.WebGLSystem,
      y.WebGPUSystem
    ],
    name: "filter"
  };
  const tt = class rt extends pe {
    constructor(...e) {
      let t = e[0] ?? {};
      t instanceof Float32Array && (U(I, "use new MeshGeometry({ positions, uvs, indices }) instead"), t = {
        positions: t,
        uvs: e[1],
        indices: e[2]
      }), t = {
        ...rt.defaultOptions,
        ...t
      };
      const r = t.positions || new Float32Array([
        0,
        0,
        1,
        0,
        1,
        1,
        0,
        1
      ]);
      let i = t.uvs;
      i || (t.positions ? i = new Float32Array(r.length) : i = new Float32Array([
        0,
        0,
        1,
        0,
        1,
        1,
        0,
        1
      ]));
      const n = t.indices || new Uint32Array([
        0,
        1,
        2,
        0,
        2,
        3
      ]), s = t.shrinkBuffersToFit, a = new Y({
        data: r,
        label: "attribute-mesh-positions",
        shrinkToFit: s,
        usage: A.VERTEX | A.COPY_DST
      }), u = new Y({
        data: i,
        label: "attribute-mesh-uvs",
        shrinkToFit: s,
        usage: A.VERTEX | A.COPY_DST
      }), l = new Y({
        data: n,
        label: "index-mesh-buffer",
        shrinkToFit: s,
        usage: A.INDEX | A.COPY_DST
      });
      super({
        attributes: {
          aPosition: {
            buffer: a,
            format: "float32x2",
            stride: 8,
            offset: 0
          },
          aUV: {
            buffer: u,
            format: "float32x2",
            stride: 8,
            offset: 0
          }
        },
        indexBuffer: l,
        topology: t.topology
      }), this.batchMode = "auto";
    }
    get positions() {
      return this.attributes.aPosition.buffer.data;
    }
    set positions(e) {
      this.attributes.aPosition.buffer.data = e;
    }
    get uvs() {
      return this.attributes.aUV.buffer.data;
    }
    set uvs(e) {
      this.attributes.aUV.buffer.data = e;
    }
    get indices() {
      return this.indexBuffer.data;
    }
    set indices(e) {
      this.indexBuffer.data = e;
    }
  };
  tt.defaultOptions = {
    topology: "triangle-list",
    shrinkBuffersToFit: false
  };
  let _e = tt, W = null, O = null;
  function er(o, e) {
    W || (W = L.get().createCanvas(256, 128), O = W.getContext("2d", {
      willReadFrequently: true
    }), O.globalCompositeOperation = "copy", O.globalAlpha = 1), (W.width < o || W.height < e) && (W.width = Te(o), W.height = Te(e));
  }
  function Fe(o, e, t) {
    for (let r = 0, i = 4 * t * e; r < e; ++r, i += 4) if (o[i + 3] !== 0) return false;
    return true;
  }
  function Be(o, e, t, r, i) {
    const n = 4 * e;
    for (let s = r, a = r * n + 4 * t; s <= i; ++s, a += n) if (o[a + 3] !== 0) return false;
    return true;
  }
  function tr(...o) {
    let e = o[0];
    e.canvas || (e = {
      canvas: o[0],
      resolution: o[1]
    });
    const { canvas: t } = e, r = Math.min(e.resolution ?? 1, 1), i = e.width ?? t.width, n = e.height ?? t.height;
    let s = e.output;
    if (er(i, n), !O) throw new TypeError("Failed to get canvas 2D context");
    O.drawImage(t, 0, 0, i, n, 0, 0, i * r, n * r);
    const u = O.getImageData(0, 0, i, n).data;
    let l = 0, h = 0, c = i - 1, d = n - 1;
    for (; h < n && Fe(u, i, h); ) ++h;
    if (h === n) return J.EMPTY;
    for (; Fe(u, i, d); ) --d;
    for (; Be(u, i, l, h, d); ) ++l;
    for (; Be(u, i, c, h, d); ) --c;
    return ++c, ++d, O.globalCompositeOperation = "source-over", O.strokeRect(l, h, c - l, d - h), O.globalCompositeOperation = "copy", s ?? (s = new J()), s.set(l / r, h / r, (c - l) / r, (d - h) / r), s;
  }
  const Re = new J();
  class rr {
    getCanvasAndContext(e) {
      const { text: t, style: r, resolution: i = 1 } = e, n = r._getFinalPadding(), s = D.measureText(t || " ", r), a = Math.ceil(Math.ceil(Math.max(1, s.width) + n * 2) * i), u = Math.ceil(Math.ceil(Math.max(1, s.height) + n * 2) * i), l = V.getOptimalCanvasAndContext(a, u);
      this._renderTextToCanvas(t, r, n, i, l);
      const h = r.trim ? tr({
        canvas: l.canvas,
        width: a,
        height: u,
        resolution: 1,
        output: Re
      }) : Re.set(0, 0, a, u);
      return {
        canvasAndContext: l,
        frame: h
      };
    }
    returnCanvasAndContext(e) {
      V.returnCanvasAndContext(e);
    }
    _renderTextToCanvas(e, t, r, i, n) {
      var _a, _b, _c, _d;
      const { canvas: s, context: a } = n, u = te(t), l = D.measureText(e || " ", t), h = l.lines, c = l.lineHeight, d = l.lineWidths, f = l.maxLineWidth, m = l.fontProperties, v = s.height;
      if (a.resetTransform(), a.scale(i, i), a.textBaseline = t.textBaseline, (_a = t._stroke) == null ? void 0 : _a.width) {
        const b = t._stroke;
        a.lineWidth = b.width, a.miterLimit = b.miterLimit, a.lineJoin = b.join, a.lineCap = b.cap;
      }
      a.font = u;
      let x, p;
      const _ = t.dropShadow ? 2 : 1;
      for (let b = 0; b < _; ++b) {
        const w = t.dropShadow && b === 0, S = w ? Math.ceil(Math.max(1, v) + r * 2) : 0, F = S * i;
        if (w) {
          a.fillStyle = "black", a.strokeStyle = "black";
          const T = t.dropShadow, K = T.color, ne = T.alpha;
          a.shadowColor = E.shared.setValue(K).setAlpha(ne).toRgbaString();
          const yt = T.blur * i, be = T.distance * i;
          a.shadowBlur = yt, a.shadowOffsetX = Math.cos(T.angle) * be, a.shadowOffsetY = Math.sin(T.angle) * be + F;
        } else {
          if (a.fillStyle = t._fill ? re(t._fill, a, l, r * 2) : null, (_b = t._stroke) == null ? void 0 : _b.width) {
            const T = t._stroke.width * 0.5 + r * 2;
            a.strokeStyle = re(t._stroke, a, l, T);
          }
          a.shadowColor = "black";
        }
        let P = (c - m.fontSize) / 2;
        c - m.fontSize < 0 && (P = 0);
        const B = ((_c = t._stroke) == null ? void 0 : _c.width) ?? 0;
        for (let T = 0; T < h.length; T++) x = B / 2, p = B / 2 + T * c + m.ascent + P, t.align === "right" ? x += f - d[T] : t.align === "center" && (x += (f - d[T]) / 2), ((_d = t._stroke) == null ? void 0 : _d.width) && this._drawLetterSpacing(h[T], t, n, x + r, p + r - S, true), t._fill !== void 0 && this._drawLetterSpacing(h[T], t, n, x + r, p + r - S);
      }
    }
    _drawLetterSpacing(e, t, r, i, n, s = false) {
      const { context: a } = r, u = t.letterSpacing;
      let l = false;
      if (D.experimentalLetterSpacingSupported && (D.experimentalLetterSpacing ? (a.letterSpacing = `${u}px`, a.textLetterSpacing = `${u}px`, l = true) : (a.letterSpacing = "0px", a.textLetterSpacing = "0px")), u === 0 || l) {
        s ? a.strokeText(e, i, n) : a.fillText(e, i, n);
        return;
      }
      let h = i;
      const c = D.graphemeSegmenter(e);
      let d = a.measureText(e).width, f = 0;
      for (let m = 0; m < c.length; ++m) {
        const v = c[m];
        s ? a.strokeText(v, h, n) : a.fillText(v, h, n);
        let x = "";
        for (let p = m + 1; p < c.length; ++p) x += c[p];
        f = a.measureText(x).width, h += d - f + u, d = f;
      }
    }
  }
  const ae = new rr(), Me = "http://www.w3.org/2000/svg", Ue = "http://www.w3.org/1999/xhtml";
  class it {
    constructor() {
      this.svgRoot = document.createElementNS(Me, "svg"), this.foreignObject = document.createElementNS(Me, "foreignObject"), this.domElement = document.createElementNS(Ue, "div"), this.styleElement = document.createElementNS(Ue, "style");
      const { foreignObject: e, svgRoot: t, styleElement: r, domElement: i } = this;
      e.setAttribute("width", "10000"), e.setAttribute("height", "10000"), e.style.overflow = "hidden", t.appendChild(e), e.appendChild(r), e.appendChild(i), this.image = L.get().createImage();
    }
    destroy() {
      this.svgRoot.remove(), this.foreignObject.remove(), this.styleElement.remove(), this.domElement.remove(), this.image.src = "", this.image.remove(), this.svgRoot = null, this.foreignObject = null, this.styleElement = null, this.domElement = null, this.image = null, this.canvasAndContext = null;
    }
  }
  let ke;
  function ir(o, e, t, r) {
    r || (r = ke || (ke = new it()));
    const { domElement: i, styleElement: n, svgRoot: s } = r;
    i.innerHTML = `<style>${e.cssStyle};</style><div style='padding:0'>${o}</div>`, i.setAttribute("style", "transform-origin: top left; display: inline-block"), t && (n.textContent = t), document.body.appendChild(s);
    const a = i.getBoundingClientRect();
    s.remove();
    const u = e.padding * 2;
    return {
      width: a.width - u,
      height: a.height - u
    };
  }
  class nr {
    constructor() {
      this.batches = [], this.batched = false;
    }
    destroy() {
      this.batches.forEach((e) => {
        ee.return(e);
      }), this.batches.length = 0;
    }
  }
  class nt {
    constructor(e, t) {
      this.state = Z.for2d(), this.renderer = e, this._adaptor = t, this.renderer.runners.contextChange.add(this);
    }
    contextChange() {
      this._adaptor.contextChange(this.renderer);
    }
    validateRenderable(e) {
      const t = e.context, r = !!e._gpuData, i = this.renderer.graphicsContext.updateGpuContext(t);
      return !!(i.isBatchable || r !== i.isBatchable);
    }
    addRenderable(e, t) {
      const r = this.renderer.graphicsContext.updateGpuContext(e.context);
      e.didViewUpdate && this._rebuild(e), r.isBatchable ? this._addToBatcher(e, t) : (this.renderer.renderPipes.batch.break(t), t.add(e));
    }
    updateRenderable(e) {
      const r = this._getGpuDataForRenderable(e).batches;
      for (let i = 0; i < r.length; i++) {
        const n = r[i];
        n._batcher.updateElement(n);
      }
    }
    execute(e) {
      if (!e.isRenderable) return;
      const t = this.renderer, r = e.context;
      if (!t.graphicsContext.getGpuContext(r).batches.length) return;
      const n = r.customShader || this._adaptor.shader;
      this.state.blendMode = e.groupBlendMode;
      const s = n.resources.localUniforms.uniforms;
      s.uTransformMatrix = e.groupTransform, s.uRound = t._roundPixels | e._roundPixels, ie(e.groupColorAlpha, s.uColor, 0), this._adaptor.execute(this, e);
    }
    _rebuild(e) {
      const t = this._getGpuDataForRenderable(e), r = this.renderer.graphicsContext.updateGpuContext(e.context);
      t.destroy(), r.isBatchable && this._updateBatchesForRenderable(e, t);
    }
    _addToBatcher(e, t) {
      const r = this.renderer.renderPipes.batch, i = this._getGpuDataForRenderable(e).batches;
      for (let n = 0; n < i.length; n++) {
        const s = i[n];
        r.addToBatch(s, t);
      }
    }
    _getGpuDataForRenderable(e) {
      return e._gpuData[this.renderer.uid] || this._initGpuDataForRenderable(e);
    }
    _initGpuDataForRenderable(e) {
      const t = new nr();
      return e._gpuData[this.renderer.uid] = t, t;
    }
    _updateBatchesForRenderable(e, t) {
      const r = e.context, i = this.renderer.graphicsContext.getGpuContext(r), n = this.renderer._roundPixels | e._roundPixels;
      t.batches = i.batches.map((s) => {
        const a = ee.get(Bt);
        return s.copyTo(a), a.renderable = e, a.roundPixels = n, a;
      });
    }
    destroy() {
      this.renderer = null, this._adaptor.destroy(), this._adaptor = null, this.state = null;
    }
  }
  nt.extension = {
    type: [
      y.WebGLPipes,
      y.WebGPUPipes,
      y.CanvasPipes
    ],
    name: "graphics"
  };
  const st = class at extends _e {
    constructor(...e) {
      super({});
      let t = e[0] ?? {};
      typeof t == "number" && (U(I, "PlaneGeometry constructor changed please use { width, height, verticesX, verticesY } instead"), t = {
        width: t,
        height: e[1],
        verticesX: e[2],
        verticesY: e[3]
      }), this.build(t);
    }
    build(e) {
      e = {
        ...at.defaultOptions,
        ...e
      }, this.verticesX = this.verticesX ?? e.verticesX, this.verticesY = this.verticesY ?? e.verticesY, this.width = this.width ?? e.width, this.height = this.height ?? e.height;
      const t = this.verticesX * this.verticesY, r = [], i = [], n = [], s = this.verticesX - 1, a = this.verticesY - 1, u = this.width / s, l = this.height / a;
      for (let c = 0; c < t; c++) {
        const d = c % this.verticesX, f = c / this.verticesX | 0;
        r.push(d * u, f * l), i.push(d / s, f / a);
      }
      const h = s * a;
      for (let c = 0; c < h; c++) {
        const d = c % s, f = c / s | 0, m = f * this.verticesX + d, v = f * this.verticesX + d + 1, x = (f + 1) * this.verticesX + d, p = (f + 1) * this.verticesX + d + 1;
        n.push(m, v, x, v, p, x);
      }
      this.buffers[0].data = new Float32Array(r), this.buffers[1].data = new Float32Array(i), this.indexBuffer.data = new Uint32Array(n), this.buffers[0].update(), this.buffers[1].update(), this.indexBuffer.update();
    }
  };
  st.defaultOptions = {
    width: 100,
    height: 100,
    verticesX: 10,
    verticesY: 10
  };
  let sr = st;
  class ve {
    constructor() {
      this.batcherName = "default", this.packAsQuad = false, this.indexOffset = 0, this.attributeOffset = 0, this.roundPixels = 0, this._batcher = null, this._batch = null, this._textureMatrixUpdateId = -1, this._uvUpdateId = -1;
    }
    get blendMode() {
      return this.renderable.groupBlendMode;
    }
    get topology() {
      return this._topology || this.geometry.topology;
    }
    set topology(e) {
      this._topology = e;
    }
    reset() {
      this.renderable = null, this.texture = null, this._batcher = null, this._batch = null, this.geometry = null, this._uvUpdateId = -1, this._textureMatrixUpdateId = -1;
    }
    setTexture(e) {
      this.texture !== e && (this.texture = e, this._textureMatrixUpdateId = -1);
    }
    get uvs() {
      const t = this.geometry.getBuffer("aUV"), r = t.data;
      let i = r;
      const n = this.texture.textureMatrix;
      return n.isSimple || (i = this._transformedUvs, (this._textureMatrixUpdateId !== n._updateID || this._uvUpdateId !== t._updateID) && ((!i || i.length < r.length) && (i = this._transformedUvs = new Float32Array(r.length)), this._textureMatrixUpdateId = n._updateID, this._uvUpdateId = t._updateID, n.multiplyUvs(r, i))), i;
    }
    get positions() {
      return this.geometry.positions;
    }
    get indices() {
      return this.geometry.indices;
    }
    get color() {
      return this.renderable.groupColorAlpha;
    }
    get groupTransform() {
      return this.renderable.groupTransform;
    }
    get attributeSize() {
      return this.geometry.positions.length / 2;
    }
    get indexSize() {
      return this.geometry.indices.length;
    }
  }
  class ze {
    destroy() {
    }
  }
  class ot {
    constructor(e, t) {
      this.localUniforms = new H({
        uTransformMatrix: {
          value: new M(),
          type: "mat3x3<f32>"
        },
        uColor: {
          value: new Float32Array([
            1,
            1,
            1,
            1
          ]),
          type: "vec4<f32>"
        },
        uRound: {
          value: 0,
          type: "f32"
        }
      }), this.localUniformsBindGroup = new Le({
        0: this.localUniforms
      }), this.renderer = e, this._adaptor = t, this._adaptor.init();
    }
    validateRenderable(e) {
      const t = this._getMeshData(e), r = t.batched, i = e.batched;
      if (t.batched = i, r !== i) return true;
      if (i) {
        const n = e._geometry;
        if (n.indices.length !== t.indexSize || n.positions.length !== t.vertexSize) return t.indexSize = n.indices.length, t.vertexSize = n.positions.length, true;
        const s = this._getBatchableMesh(e);
        return s.texture.uid !== e._texture.uid && (s._textureMatrixUpdateId = -1), !s._batcher.checkAndUpdateTexture(s, e._texture);
      }
      return false;
    }
    addRenderable(e, t) {
      var _a, _b;
      const r = this.renderer.renderPipes.batch, i = this._getMeshData(e);
      if (e.didViewUpdate && (i.indexSize = (_a = e._geometry.indices) == null ? void 0 : _a.length, i.vertexSize = (_b = e._geometry.positions) == null ? void 0 : _b.length), i.batched) {
        const n = this._getBatchableMesh(e);
        n.setTexture(e._texture), n.geometry = e._geometry, r.addToBatch(n, t);
      } else r.break(t), t.add(e);
    }
    updateRenderable(e) {
      if (e.batched) {
        const t = this._getBatchableMesh(e);
        t.setTexture(e._texture), t.geometry = e._geometry, t._batcher.updateElement(t);
      }
    }
    execute(e) {
      if (!e.isRenderable) return;
      e.state.blendMode = me(e.groupBlendMode, e.texture._source);
      const t = this.localUniforms;
      t.uniforms.uTransformMatrix = e.groupTransform, t.uniforms.uRound = this.renderer._roundPixels | e._roundPixels, t.update(), ie(e.groupColorAlpha, t.uniforms.uColor, 0), this._adaptor.execute(this, e);
    }
    _getMeshData(e) {
      var t, r;
      return (t = e._gpuData)[r = this.renderer.uid] || (t[r] = new ze()), e._gpuData[this.renderer.uid].meshData || this._initMeshData(e);
    }
    _initMeshData(e) {
      return e._gpuData[this.renderer.uid].meshData = {
        batched: e.batched,
        indexSize: 0,
        vertexSize: 0
      }, e._gpuData[this.renderer.uid].meshData;
    }
    _getBatchableMesh(e) {
      var t, r;
      return (t = e._gpuData)[r = this.renderer.uid] || (t[r] = new ze()), e._gpuData[this.renderer.uid].batchableMesh || this._initBatchableMesh(e);
    }
    _initBatchableMesh(e) {
      const t = new ve();
      return t.renderable = e, t.setTexture(e._texture), t.transform = e.groupTransform, t.roundPixels = this.renderer._roundPixels | e._roundPixels, e._gpuData[this.renderer.uid].batchableMesh = t, t;
    }
    destroy() {
      this.localUniforms = null, this.localUniformsBindGroup = null, this._adaptor.destroy(), this._adaptor = null, this.renderer = null;
    }
  }
  ot.extension = {
    type: [
      y.WebGLPipes,
      y.WebGPUPipes,
      y.CanvasPipes
    ],
    name: "mesh"
  };
  class ar {
    execute(e, t) {
      const r = e.state, i = e.renderer, n = t.shader || e.defaultShader;
      n.resources.uTexture = t.texture._source, n.resources.uniforms = e.localUniforms;
      const s = i.gl, a = e.getBuffers(t);
      i.shader.bind(n), i.state.set(r), i.geometry.bind(a.geometry, n.glProgram);
      const l = a.geometry.indexBuffer.data.BYTES_PER_ELEMENT === 2 ? s.UNSIGNED_SHORT : s.UNSIGNED_INT;
      s.drawElements(s.TRIANGLES, t.particleChildren.length * 6, l, 0);
    }
  }
  class or {
    execute(e, t) {
      const r = e.renderer, i = t.shader || e.defaultShader;
      i.groups[0] = r.renderPipes.uniformBatch.getUniformBindGroup(e.localUniforms, true), i.groups[1] = r.texture.getTextureBindGroup(t.texture);
      const n = e.state, s = e.getBuffers(t);
      r.encoder.draw({
        geometry: s.geometry,
        shader: t.shader || e.defaultShader,
        state: n,
        size: t.particleChildren.length * 6
      });
    }
  }
  function Ge(o, e = null) {
    const t = o * 6;
    if (t > 65535 ? e || (e = new Uint32Array(t)) : e || (e = new Uint16Array(t)), e.length !== t) throw new Error(`Out buffer length is incorrect, got ${e.length} and expected ${t}`);
    for (let r = 0, i = 0; r < t; r += 6, i += 4) e[r + 0] = i + 0, e[r + 1] = i + 1, e[r + 2] = i + 2, e[r + 3] = i + 0, e[r + 4] = i + 2, e[r + 5] = i + 3;
    return e;
  }
  function ur(o) {
    return {
      dynamicUpdate: Ae(o, true),
      staticUpdate: Ae(o, false)
    };
  }
  function Ae(o, e) {
    const t = [];
    t.push(`

        var index = 0;

        for (let i = 0; i < ps.length; ++i)
        {
            const p = ps[i];

            `);
    let r = 0;
    for (const n in o) {
      const s = o[n];
      if (e !== s.dynamic) continue;
      t.push(`offset = index + ${r}`), t.push(s.code);
      const a = de(s.format);
      r += a.stride / 4;
    }
    t.push(`
            index += stride * 4;
        }
    `), t.unshift(`
        var stride = ${r};
    `);
    const i = t.join(`
`);
    return new Function("ps", "f32v", "u32v", i);
  }
  class lr {
    constructor(e) {
      this._size = 0, this._generateParticleUpdateCache = {};
      const t = this._size = e.size ?? 1e3, r = e.properties;
      let i = 0, n = 0;
      for (const h in r) {
        const c = r[h], d = de(c.format);
        c.dynamic ? n += d.stride : i += d.stride;
      }
      this._dynamicStride = n / 4, this._staticStride = i / 4, this.staticAttributeBuffer = new j(t * 4 * i), this.dynamicAttributeBuffer = new j(t * 4 * n), this.indexBuffer = Ge(t);
      const s = new pe();
      let a = 0, u = 0;
      this._staticBuffer = new Y({
        data: new Float32Array(1),
        label: "static-particle-buffer",
        shrinkToFit: false,
        usage: A.VERTEX | A.COPY_DST
      }), this._dynamicBuffer = new Y({
        data: new Float32Array(1),
        label: "dynamic-particle-buffer",
        shrinkToFit: false,
        usage: A.VERTEX | A.COPY_DST
      });
      for (const h in r) {
        const c = r[h], d = de(c.format);
        c.dynamic ? (s.addAttribute(c.attributeName, {
          buffer: this._dynamicBuffer,
          stride: this._dynamicStride * 4,
          offset: a * 4,
          format: c.format
        }), a += d.size) : (s.addAttribute(c.attributeName, {
          buffer: this._staticBuffer,
          stride: this._staticStride * 4,
          offset: u * 4,
          format: c.format
        }), u += d.size);
      }
      s.addIndex(this.indexBuffer);
      const l = this.getParticleUpdate(r);
      this._dynamicUpload = l.dynamicUpdate, this._staticUpload = l.staticUpdate, this.geometry = s;
    }
    getParticleUpdate(e) {
      const t = cr(e);
      return this._generateParticleUpdateCache[t] ? this._generateParticleUpdateCache[t] : (this._generateParticleUpdateCache[t] = this.generateParticleUpdate(e), this._generateParticleUpdateCache[t]);
    }
    generateParticleUpdate(e) {
      return ur(e);
    }
    update(e, t) {
      e.length > this._size && (t = true, this._size = Math.max(e.length, this._size * 1.5 | 0), this.staticAttributeBuffer = new j(this._size * this._staticStride * 4 * 4), this.dynamicAttributeBuffer = new j(this._size * this._dynamicStride * 4 * 4), this.indexBuffer = Ge(this._size), this.geometry.indexBuffer.setDataWithSize(this.indexBuffer, this.indexBuffer.byteLength, true));
      const r = this.dynamicAttributeBuffer;
      if (this._dynamicUpload(e, r.float32View, r.uint32View), this._dynamicBuffer.setDataWithSize(this.dynamicAttributeBuffer.float32View, e.length * this._dynamicStride * 4, true), t) {
        const i = this.staticAttributeBuffer;
        this._staticUpload(e, i.float32View, i.uint32View), this._staticBuffer.setDataWithSize(i.float32View, e.length * this._staticStride * 4, true);
      }
    }
    destroy() {
      this._staticBuffer.destroy(), this._dynamicBuffer.destroy(), this.geometry.destroy();
    }
  }
  function cr(o) {
    const e = [];
    for (const t in o) {
      const r = o[t];
      e.push(t, r.code, r.dynamic ? "d" : "s");
    }
    return e.join("_");
  }
  var hr = `varying vec2 vUV;
varying vec4 vColor;

uniform sampler2D uTexture;

void main(void){
    vec4 color = texture2D(uTexture, vUV) * vColor;
    gl_FragColor = color;
}`, dr = `attribute vec2 aVertex;
attribute vec2 aUV;
attribute vec4 aColor;

attribute vec2 aPosition;
attribute float aRotation;

uniform mat3 uTranslationMatrix;
uniform float uRound;
uniform vec2 uResolution;
uniform vec4 uColor;

varying vec2 vUV;
varying vec4 vColor;

vec2 roundPixels(vec2 position, vec2 targetSize)
{       
    return (floor(((position * 0.5 + 0.5) * targetSize) + 0.5) / targetSize) * 2.0 - 1.0;
}

void main(void){
    float cosRotation = cos(aRotation);
    float sinRotation = sin(aRotation);
    float x = aVertex.x * cosRotation - aVertex.y * sinRotation;
    float y = aVertex.x * sinRotation + aVertex.y * cosRotation;

    vec2 v = vec2(x, y);
    v = v + aPosition;

    gl_Position = vec4((uTranslationMatrix * vec3(v, 1.0)).xy, 0.0, 1.0);

    if(uRound == 1.0)
    {
        gl_Position.xy = roundPixels(gl_Position.xy, uResolution);
    }

    vUV = aUV;
    vColor = vec4(aColor.rgb * aColor.a, aColor.a) * uColor;
}
`, De = `
struct ParticleUniforms {
  uTranslationMatrix:mat3x3<f32>,
  uColor:vec4<f32>,
  uRound:f32,
  uResolution:vec2<f32>,
};

fn roundPixels(position: vec2<f32>, targetSize: vec2<f32>) -> vec2<f32>
{
  return (floor(((position * 0.5 + 0.5) * targetSize) + 0.5) / targetSize) * 2.0 - 1.0;
}

@group(0) @binding(0) var<uniform> uniforms: ParticleUniforms;

@group(1) @binding(0) var uTexture: texture_2d<f32>;
@group(1) @binding(1) var uSampler : sampler;

struct VSOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv : vec2<f32>,
    @location(1) color : vec4<f32>,
  };
@vertex
fn mainVertex(
  @location(0) aVertex: vec2<f32>,
  @location(1) aPosition: vec2<f32>,
  @location(2) aUV: vec2<f32>,
  @location(3) aColor: vec4<f32>,
  @location(4) aRotation: f32,
) -> VSOutput {
  
   let v = vec2(
       aVertex.x * cos(aRotation) - aVertex.y * sin(aRotation),
       aVertex.x * sin(aRotation) + aVertex.y * cos(aRotation)
   ) + aPosition;

   var position = vec4((uniforms.uTranslationMatrix * vec3(v, 1.0)).xy, 0.0, 1.0);

   if(uniforms.uRound == 1.0) {
       position = vec4(roundPixels(position.xy, uniforms.uResolution), position.zw);
   }

    let vColor = vec4(aColor.rgb * aColor.a, aColor.a) * uniforms.uColor;

  return VSOutput(
   position,
   aUV,
   vColor,
  );
}

@fragment
fn mainFragment(
  @location(0) uv: vec2<f32>,
  @location(1) color: vec4<f32>,
  @builtin(position) position: vec4<f32>,
) -> @location(0) vec4<f32> {

    var sample = textureSample(uTexture, uSampler, uv) * color;
   
    return sample;
}`;
  class fr extends xe {
    constructor() {
      const e = We.from({
        vertex: dr,
        fragment: hr
      }), t = Ie.from({
        fragment: {
          source: De,
          entryPoint: "mainFragment"
        },
        vertex: {
          source: De,
          entryPoint: "mainVertex"
        }
      });
      super({
        glProgram: e,
        gpuProgram: t,
        resources: {
          uTexture: R.WHITE.source,
          uSampler: new $({}),
          uniforms: {
            uTranslationMatrix: {
              value: new M(),
              type: "mat3x3<f32>"
            },
            uColor: {
              value: new E(16777215),
              type: "vec4<f32>"
            },
            uRound: {
              value: 1,
              type: "f32"
            },
            uResolution: {
              value: [
                0,
                0
              ],
              type: "vec2<f32>"
            }
          }
        }
      });
    }
  }
  class ut {
    constructor(e, t) {
      this.state = Z.for2d(), this.localUniforms = new H({
        uTranslationMatrix: {
          value: new M(),
          type: "mat3x3<f32>"
        },
        uColor: {
          value: new Float32Array(4),
          type: "vec4<f32>"
        },
        uRound: {
          value: 1,
          type: "f32"
        },
        uResolution: {
          value: [
            0,
            0
          ],
          type: "vec2<f32>"
        }
      }), this.renderer = e, this.adaptor = t, this.defaultShader = new fr(), this.state = Z.for2d();
    }
    validateRenderable(e) {
      return false;
    }
    addRenderable(e, t) {
      this.renderer.renderPipes.batch.break(t), t.add(e);
    }
    getBuffers(e) {
      return e._gpuData[this.renderer.uid] || this._initBuffer(e);
    }
    _initBuffer(e) {
      return e._gpuData[this.renderer.uid] = new lr({
        size: e.particleChildren.length,
        properties: e._properties
      }), e._gpuData[this.renderer.uid];
    }
    updateRenderable(e) {
    }
    execute(e) {
      const t = e.particleChildren;
      if (t.length === 0) return;
      const r = this.renderer, i = this.getBuffers(e);
      e.texture || (e.texture = t[0].texture);
      const n = this.state;
      i.update(t, e._childrenDirty), e._childrenDirty = false, n.blendMode = me(e.blendMode, e.texture._source);
      const s = this.localUniforms.uniforms, a = s.uTranslationMatrix;
      e.worldTransform.copyTo(a), a.prepend(r.globalUniforms.globalUniformData.projectionMatrix), s.uResolution = r.globalUniforms.globalUniformData.resolution, s.uRound = r._roundPixels | e._roundPixels, ie(e.groupColorAlpha, s.uColor, 0), this.adaptor.execute(this, e);
    }
    destroy() {
      this.renderer = null, this.defaultShader && (this.defaultShader.destroy(), this.defaultShader = null);
    }
  }
  class lt extends ut {
    constructor(e) {
      super(e, new ar());
    }
  }
  lt.extension = {
    type: [
      y.WebGLPipes
    ],
    name: "particle"
  };
  class ct extends ut {
    constructor(e) {
      super(e, new or());
    }
  }
  ct.extension = {
    type: [
      y.WebGPUPipes
    ],
    name: "particle"
  };
  const ht = class dt extends sr {
    constructor(e = {}) {
      e = {
        ...dt.defaultOptions,
        ...e
      }, super({
        width: e.width,
        height: e.height,
        verticesX: 4,
        verticesY: 4
      }), this.update(e);
    }
    update(e) {
      var _a, _b;
      this.width = e.width ?? this.width, this.height = e.height ?? this.height, this._originalWidth = e.originalWidth ?? this._originalWidth, this._originalHeight = e.originalHeight ?? this._originalHeight, this._leftWidth = e.leftWidth ?? this._leftWidth, this._rightWidth = e.rightWidth ?? this._rightWidth, this._topHeight = e.topHeight ?? this._topHeight, this._bottomHeight = e.bottomHeight ?? this._bottomHeight, this._anchorX = (_a = e.anchor) == null ? void 0 : _a.x, this._anchorY = (_b = e.anchor) == null ? void 0 : _b.y, this.updateUvs(), this.updatePositions();
    }
    updatePositions() {
      const e = this.positions, { width: t, height: r, _leftWidth: i, _rightWidth: n, _topHeight: s, _bottomHeight: a, _anchorX: u, _anchorY: l } = this, h = i + n, c = t > h ? 1 : t / h, d = s + a, f = r > d ? 1 : r / d, m = Math.min(c, f), v = u * t, x = l * r;
      e[0] = e[8] = e[16] = e[24] = -v, e[2] = e[10] = e[18] = e[26] = i * m - v, e[4] = e[12] = e[20] = e[28] = t - n * m - v, e[6] = e[14] = e[22] = e[30] = t - v, e[1] = e[3] = e[5] = e[7] = -x, e[9] = e[11] = e[13] = e[15] = s * m - x, e[17] = e[19] = e[21] = e[23] = r - a * m - x, e[25] = e[27] = e[29] = e[31] = r - x, this.getBuffer("aPosition").update();
    }
    updateUvs() {
      const e = this.uvs;
      e[0] = e[8] = e[16] = e[24] = 0, e[1] = e[3] = e[5] = e[7] = 0, e[6] = e[14] = e[22] = e[30] = 1, e[25] = e[27] = e[29] = e[31] = 1;
      const t = 1 / this._originalWidth, r = 1 / this._originalHeight;
      e[2] = e[10] = e[18] = e[26] = t * this._leftWidth, e[9] = e[11] = e[13] = e[15] = r * this._topHeight, e[4] = e[12] = e[20] = e[28] = 1 - t * this._rightWidth, e[17] = e[19] = e[21] = e[23] = 1 - r * this._bottomHeight, this.getBuffer("aUV").update();
    }
  };
  ht.defaultOptions = {
    width: 100,
    height: 100,
    leftWidth: 10,
    topHeight: 10,
    rightWidth: 10,
    bottomHeight: 10,
    originalWidth: 100,
    originalHeight: 100
  };
  let pr = ht;
  class gr extends ve {
    constructor() {
      super(), this.geometry = new pr();
    }
    destroy() {
      this.geometry.destroy();
    }
  }
  class ft {
    constructor(e) {
      this._renderer = e;
    }
    addRenderable(e, t) {
      const r = this._getGpuSprite(e);
      e.didViewUpdate && this._updateBatchableSprite(e, r), this._renderer.renderPipes.batch.addToBatch(r, t);
    }
    updateRenderable(e) {
      const t = this._getGpuSprite(e);
      e.didViewUpdate && this._updateBatchableSprite(e, t), t._batcher.updateElement(t);
    }
    validateRenderable(e) {
      const t = this._getGpuSprite(e);
      return !t._batcher.checkAndUpdateTexture(t, e._texture);
    }
    _updateBatchableSprite(e, t) {
      t.geometry.update(e), t.setTexture(e._texture);
    }
    _getGpuSprite(e) {
      return e._gpuData[this._renderer.uid] || this._initGPUSprite(e);
    }
    _initGPUSprite(e) {
      const t = e._gpuData[this._renderer.uid] = new gr(), r = t;
      return r.renderable = e, r.transform = e.groupTransform, r.texture = e._texture, r.roundPixels = this._renderer._roundPixels | e._roundPixels, e.didViewUpdate || this._updateBatchableSprite(e, r), t;
    }
    destroy() {
      this._renderer = null;
    }
  }
  ft.extension = {
    type: [
      y.WebGLPipes,
      y.WebGPUPipes,
      y.CanvasPipes
    ],
    name: "nineSliceSprite"
  };
  const mr = {
    name: "tiling-bit",
    vertex: {
      header: `
            struct TilingUniforms {
                uMapCoord:mat3x3<f32>,
                uClampFrame:vec4<f32>,
                uClampOffset:vec2<f32>,
                uTextureTransform:mat3x3<f32>,
                uSizeAnchor:vec4<f32>
            };

            @group(2) @binding(0) var<uniform> tilingUniforms: TilingUniforms;
            @group(2) @binding(1) var uTexture: texture_2d<f32>;
            @group(2) @binding(2) var uSampler: sampler;
        `,
      main: `
            uv = (tilingUniforms.uTextureTransform * vec3(uv, 1.0)).xy;

            position = (position - tilingUniforms.uSizeAnchor.zw) * tilingUniforms.uSizeAnchor.xy;
        `
    },
    fragment: {
      header: `
            struct TilingUniforms {
                uMapCoord:mat3x3<f32>,
                uClampFrame:vec4<f32>,
                uClampOffset:vec2<f32>,
                uTextureTransform:mat3x3<f32>,
                uSizeAnchor:vec4<f32>
            };

            @group(2) @binding(0) var<uniform> tilingUniforms: TilingUniforms;
            @group(2) @binding(1) var uTexture: texture_2d<f32>;
            @group(2) @binding(2) var uSampler: sampler;
        `,
      main: `

            var coord = vUV + ceil(tilingUniforms.uClampOffset - vUV);
            coord = (tilingUniforms.uMapCoord * vec3(coord, 1.0)).xy;
            var unclamped = coord;
            coord = clamp(coord, tilingUniforms.uClampFrame.xy, tilingUniforms.uClampFrame.zw);

            var bias = 0.;

            if(unclamped.x == coord.x && unclamped.y == coord.y)
            {
                bias = -32.;
            }

            outColor = textureSampleBias(uTexture, uSampler, coord, bias);
        `
    }
  }, xr = {
    name: "tiling-bit",
    vertex: {
      header: `
            uniform mat3 uTextureTransform;
            uniform vec4 uSizeAnchor;

        `,
      main: `
            uv = (uTextureTransform * vec3(aUV, 1.0)).xy;

            position = (position - uSizeAnchor.zw) * uSizeAnchor.xy;
        `
    },
    fragment: {
      header: `
            uniform sampler2D uTexture;
            uniform mat3 uMapCoord;
            uniform vec4 uClampFrame;
            uniform vec2 uClampOffset;
        `,
      main: `

        vec2 coord = vUV + ceil(uClampOffset - vUV);
        coord = (uMapCoord * vec3(coord, 1.0)).xy;
        vec2 unclamped = coord;
        coord = clamp(coord, uClampFrame.xy, uClampFrame.zw);

        outColor = texture(uTexture, coord, unclamped == coord ? 0.0 : -32.0);// lod-bias very negative to force lod 0

        `
    }
  };
  let oe, ue;
  class _r extends xe {
    constructor() {
      oe ?? (oe = He({
        name: "tiling-sprite-shader",
        bits: [
          Ot,
          mr,
          Ve
        ]
      })), ue ?? (ue = Ye({
        name: "tiling-sprite-shader",
        bits: [
          It,
          xr,
          Xe
        ]
      }));
      const e = new H({
        uMapCoord: {
          value: new M(),
          type: "mat3x3<f32>"
        },
        uClampFrame: {
          value: new Float32Array([
            0,
            0,
            1,
            1
          ]),
          type: "vec4<f32>"
        },
        uClampOffset: {
          value: new Float32Array([
            0,
            0
          ]),
          type: "vec2<f32>"
        },
        uTextureTransform: {
          value: new M(),
          type: "mat3x3<f32>"
        },
        uSizeAnchor: {
          value: new Float32Array([
            100,
            100,
            0.5,
            0.5
          ]),
          type: "vec4<f32>"
        }
      });
      super({
        glProgram: ue,
        gpuProgram: oe,
        resources: {
          localUniforms: new H({
            uTransformMatrix: {
              value: new M(),
              type: "mat3x3<f32>"
            },
            uColor: {
              value: new Float32Array([
                1,
                1,
                1,
                1
              ]),
              type: "vec4<f32>"
            },
            uRound: {
              value: 0,
              type: "f32"
            }
          }),
          tilingUniforms: e,
          uTexture: R.EMPTY.source,
          uSampler: R.EMPTY.source.style
        }
      });
    }
    updateUniforms(e, t, r, i, n, s) {
      const a = this.resources.tilingUniforms, u = s.width, l = s.height, h = s.textureMatrix, c = a.uniforms.uTextureTransform;
      c.set(r.a * u / e, r.b * u / t, r.c * l / e, r.d * l / t, r.tx / e, r.ty / t), c.invert(), a.uniforms.uMapCoord = h.mapCoord, a.uniforms.uClampFrame = h.uClampFrame, a.uniforms.uClampOffset = h.uClampOffset, a.uniforms.uTextureTransform = c, a.uniforms.uSizeAnchor[0] = e, a.uniforms.uSizeAnchor[1] = t, a.uniforms.uSizeAnchor[2] = i, a.uniforms.uSizeAnchor[3] = n, s && (this.resources.uTexture = s.source, this.resources.uSampler = s.source.style);
    }
  }
  class vr extends _e {
    constructor() {
      super({
        positions: new Float32Array([
          0,
          0,
          1,
          0,
          1,
          1,
          0,
          1
        ]),
        uvs: new Float32Array([
          0,
          0,
          1,
          0,
          1,
          1,
          0,
          1
        ]),
        indices: new Uint32Array([
          0,
          1,
          2,
          0,
          2,
          3
        ])
      });
    }
  }
  function br(o, e) {
    const t = o.anchor.x, r = o.anchor.y;
    e[0] = -t * o.width, e[1] = -r * o.height, e[2] = (1 - t) * o.width, e[3] = -r * o.height, e[4] = (1 - t) * o.width, e[5] = (1 - r) * o.height, e[6] = -t * o.width, e[7] = (1 - r) * o.height;
  }
  function yr(o, e, t, r) {
    let i = 0;
    const n = o.length / e, s = r.a, a = r.b, u = r.c, l = r.d, h = r.tx, c = r.ty;
    for (t *= e; i < n; ) {
      const d = o[t], f = o[t + 1];
      o[t] = s * d + u * f + h, o[t + 1] = a * d + l * f + c, t += e, i++;
    }
  }
  function Tr(o, e) {
    const t = o.texture, r = t.frame.width, i = t.frame.height;
    let n = 0, s = 0;
    o.applyAnchorToTexture && (n = o.anchor.x, s = o.anchor.y), e[0] = e[6] = -n, e[2] = e[4] = 1 - n, e[1] = e[3] = -s, e[5] = e[7] = 1 - s;
    const a = M.shared;
    a.copyFrom(o._tileTransform.matrix), a.tx /= o.width, a.ty /= o.height, a.invert(), a.scale(o.width / r, o.height / i), yr(e, 2, 0, a);
  }
  const Q = new vr();
  class wr {
    constructor() {
      this.canBatch = true, this.geometry = new _e({
        indices: Q.indices.slice(),
        positions: Q.positions.slice(),
        uvs: Q.uvs.slice()
      });
    }
    destroy() {
      var _a;
      this.geometry.destroy(), (_a = this.shader) == null ? void 0 : _a.destroy();
    }
  }
  class pt {
    constructor(e) {
      this._state = Z.default2d, this._renderer = e;
    }
    validateRenderable(e) {
      const t = this._getTilingSpriteData(e), r = t.canBatch;
      this._updateCanBatch(e);
      const i = t.canBatch;
      if (i && i === r) {
        const { batchableMesh: n } = t;
        return !n._batcher.checkAndUpdateTexture(n, e.texture);
      }
      return r !== i;
    }
    addRenderable(e, t) {
      const r = this._renderer.renderPipes.batch;
      this._updateCanBatch(e);
      const i = this._getTilingSpriteData(e), { geometry: n, canBatch: s } = i;
      if (s) {
        i.batchableMesh || (i.batchableMesh = new ve());
        const a = i.batchableMesh;
        e.didViewUpdate && (this._updateBatchableMesh(e), a.geometry = n, a.renderable = e, a.transform = e.groupTransform, a.setTexture(e._texture)), a.roundPixels = this._renderer._roundPixels | e._roundPixels, r.addToBatch(a, t);
      } else r.break(t), i.shader || (i.shader = new _r()), this.updateRenderable(e), t.add(e);
    }
    execute(e) {
      const { shader: t } = this._getTilingSpriteData(e);
      t.groups[0] = this._renderer.globalUniforms.bindGroup;
      const r = t.resources.localUniforms.uniforms;
      r.uTransformMatrix = e.groupTransform, r.uRound = this._renderer._roundPixels | e._roundPixels, ie(e.groupColorAlpha, r.uColor, 0), this._state.blendMode = me(e.groupBlendMode, e.texture._source), this._renderer.encoder.draw({
        geometry: Q,
        shader: t,
        state: this._state
      });
    }
    updateRenderable(e) {
      const t = this._getTilingSpriteData(e), { canBatch: r } = t;
      if (r) {
        const { batchableMesh: i } = t;
        e.didViewUpdate && this._updateBatchableMesh(e), i._batcher.updateElement(i);
      } else if (e.didViewUpdate) {
        const { shader: i } = t;
        i.updateUniforms(e.width, e.height, e._tileTransform.matrix, e.anchor.x, e.anchor.y, e.texture);
      }
    }
    _getTilingSpriteData(e) {
      return e._gpuData[this._renderer.uid] || this._initTilingSpriteData(e);
    }
    _initTilingSpriteData(e) {
      const t = new wr();
      return t.renderable = e, e._gpuData[this._renderer.uid] = t, t;
    }
    _updateBatchableMesh(e) {
      const t = this._getTilingSpriteData(e), { geometry: r } = t, i = e.texture.source.style;
      i.addressMode !== "repeat" && (i.addressMode = "repeat", i.update()), Tr(e, r.uvs), br(e, r.positions);
    }
    destroy() {
      this._renderer = null;
    }
    _updateCanBatch(e) {
      const t = this._getTilingSpriteData(e), r = e.texture;
      let i = true;
      return this._renderer.type === ge.WEBGL && (i = this._renderer.context.supports.nonPowOf2wrapping), t.canBatch = r.textureMatrix.isSimple && (i || r.source.isPowerOfTwo), t.canBatch;
    }
  }
  pt.extension = {
    type: [
      y.WebGLPipes,
      y.WebGPUPipes,
      y.CanvasPipes
    ],
    name: "tilingSprite"
  };
  const Sr = {
    name: "local-uniform-msdf-bit",
    vertex: {
      header: `
            struct LocalUniforms {
                uColor:vec4<f32>,
                uTransformMatrix:mat3x3<f32>,
                uDistance: f32,
                uRound:f32,
            }

            @group(2) @binding(0) var<uniform> localUniforms : LocalUniforms;
        `,
      main: `
            vColor *= localUniforms.uColor;
            modelMatrix *= localUniforms.uTransformMatrix;
        `,
      end: `
            if(localUniforms.uRound == 1)
            {
                vPosition = vec4(roundPixels(vPosition.xy, globalUniforms.uResolution), vPosition.zw);
            }
        `
    },
    fragment: {
      header: `
            struct LocalUniforms {
                uColor:vec4<f32>,
                uTransformMatrix:mat3x3<f32>,
                uDistance: f32
            }

            @group(2) @binding(0) var<uniform> localUniforms : LocalUniforms;
         `,
      main: `
            outColor = vec4<f32>(calculateMSDFAlpha(outColor, localUniforms.uColor, localUniforms.uDistance));
        `
    }
  }, Cr = {
    name: "local-uniform-msdf-bit",
    vertex: {
      header: `
            uniform mat3 uTransformMatrix;
            uniform vec4 uColor;
            uniform float uRound;
        `,
      main: `
            vColor *= uColor;
            modelMatrix *= uTransformMatrix;
        `,
      end: `
            if(uRound == 1.)
            {
                gl_Position.xy = roundPixels(gl_Position.xy, uResolution);
            }
        `
    },
    fragment: {
      header: `
            uniform float uDistance;
         `,
      main: `
            outColor = vec4(calculateMSDFAlpha(outColor, vColor, uDistance));
        `
    }
  }, Pr = {
    name: "msdf-bit",
    fragment: {
      header: `
            fn calculateMSDFAlpha(msdfColor:vec4<f32>, shapeColor:vec4<f32>, distance:f32) -> f32 {

                // MSDF
                var median = msdfColor.r + msdfColor.g + msdfColor.b -
                    min(msdfColor.r, min(msdfColor.g, msdfColor.b)) -
                    max(msdfColor.r, max(msdfColor.g, msdfColor.b));

                // SDF
                median = min(median, msdfColor.a);

                var screenPxDistance = distance * (median - 0.5);
                var alpha = clamp(screenPxDistance + 0.5, 0.0, 1.0);
                if (median < 0.01) {
                    alpha = 0.0;
                } else if (median > 0.99) {
                    alpha = 1.0;
                }

                // Gamma correction for coverage-like alpha
                var luma: f32 = dot(shapeColor.rgb, vec3<f32>(0.299, 0.587, 0.114));
                var gamma: f32 = mix(1.0, 1.0 / 2.2, luma);
                var coverage: f32 = pow(shapeColor.a * alpha, gamma);

                return coverage;

            }
        `
    }
  }, Fr = {
    name: "msdf-bit",
    fragment: {
      header: `
            float calculateMSDFAlpha(vec4 msdfColor, vec4 shapeColor, float distance) {

                // MSDF
                float median = msdfColor.r + msdfColor.g + msdfColor.b -
                                min(msdfColor.r, min(msdfColor.g, msdfColor.b)) -
                                max(msdfColor.r, max(msdfColor.g, msdfColor.b));

                // SDF
                median = min(median, msdfColor.a);

                float screenPxDistance = distance * (median - 0.5);
                float alpha = clamp(screenPxDistance + 0.5, 0.0, 1.0);

                if (median < 0.01) {
                    alpha = 0.0;
                } else if (median > 0.99) {
                    alpha = 1.0;
                }

                // Gamma correction for coverage-like alpha
                float luma = dot(shapeColor.rgb, vec3(0.299, 0.587, 0.114));
                float gamma = mix(1.0, 1.0 / 2.2, luma);
                float coverage = pow(shapeColor.a * alpha, gamma);

                return coverage;
            }
        `
    }
  };
  let le, ce;
  class Br extends xe {
    constructor(e) {
      const t = new H({
        uColor: {
          value: new Float32Array([
            1,
            1,
            1,
            1
          ]),
          type: "vec4<f32>"
        },
        uTransformMatrix: {
          value: new M(),
          type: "mat3x3<f32>"
        },
        uDistance: {
          value: 4,
          type: "f32"
        },
        uRound: {
          value: 0,
          type: "f32"
        }
      });
      le ?? (le = He({
        name: "sdf-shader",
        bits: [
          Rt,
          Mt(e),
          Sr,
          Pr,
          Ve
        ]
      })), ce ?? (ce = Ye({
        name: "sdf-shader",
        bits: [
          Ut,
          kt(e),
          Cr,
          Fr,
          Xe
        ]
      })), super({
        glProgram: ce,
        gpuProgram: le,
        resources: {
          localUniforms: t,
          batchSamplers: zt(e)
        }
      });
    }
  }
  class Rr extends Gt {
    destroy() {
      this.context.customShader && this.context.customShader.destroy(), super.destroy();
    }
  }
  class gt {
    constructor(e) {
      this._renderer = e;
    }
    validateRenderable(e) {
      const t = this._getGpuBitmapText(e);
      return this._renderer.renderPipes.graphics.validateRenderable(t);
    }
    addRenderable(e, t) {
      const r = this._getGpuBitmapText(e);
      Oe(e, r), e._didTextUpdate && (e._didTextUpdate = false, this._updateContext(e, r)), this._renderer.renderPipes.graphics.addRenderable(r, t), r.context.customShader && this._updateDistanceField(e);
    }
    updateRenderable(e) {
      const t = this._getGpuBitmapText(e);
      Oe(e, t), this._renderer.renderPipes.graphics.updateRenderable(t), t.context.customShader && this._updateDistanceField(e);
    }
    _updateContext(e, t) {
      const { context: r } = t, i = Nt.getFont(e.text, e._style);
      r.clear(), i.distanceField.type !== "none" && (r.customShader || (r.customShader = new Br(this._renderer.limits.maxBatchableTextures)));
      const n = D.graphemeSegmenter(e.text), s = e._style;
      let a = i.baseLineOffset;
      const u = Je(n, s, i, true), l = s.padding, h = u.scale;
      let c = u.width, d = u.height + u.offsetY;
      s._stroke && (c += s._stroke.width / h, d += s._stroke.width / h), r.translate(-e._anchor._x * c - l, -e._anchor._y * d - l).scale(h, h);
      const f = i.applyFillAsTint ? s._fill.color : 16777215;
      let m = i.fontMetrics.fontSize, v = i.lineHeight;
      s.lineHeight && (m = s.fontSize / h, v = s.lineHeight / h);
      let x = (v - m) / 2;
      x - i.baseLineOffset < 0 && (x = 0);
      for (let p = 0; p < u.lines.length; p++) {
        const _ = u.lines[p];
        for (let b = 0; b < _.charPositions.length; b++) {
          const w = _.chars[b], S = i.chars[w];
          if (S == null ? void 0 : S.texture) {
            const F = S.texture;
            r.texture(F, f || "black", Math.round(_.charPositions[b] + S.xOffset), Math.round(a + S.yOffset + x), F.orig.width, F.orig.height);
          }
        }
        a += v;
      }
    }
    _getGpuBitmapText(e) {
      return e._gpuData[this._renderer.uid] || this.initGpuText(e);
    }
    initGpuText(e) {
      const t = new Rr();
      return e._gpuData[this._renderer.uid] = t, this._updateContext(e, t), t;
    }
    _updateDistanceField(e) {
      const t = this._getGpuBitmapText(e).context, r = e._style.fontFamily, i = z.get(`${r}-bitmap`), { a: n, b: s, c: a, d: u } = e.groupTransform, l = Math.sqrt(n * n + s * s), h = Math.sqrt(a * a + u * u), c = (Math.abs(l) + Math.abs(h)) / 2, d = i.baseRenderedFontSize / e._style.fontSize, f = c * i.distanceField.range * (1 / d);
      t.customShader.resources.localUniforms.uniforms.uDistance = f;
    }
    destroy() {
      this._renderer = null;
    }
  }
  gt.extension = {
    type: [
      y.WebGLPipes,
      y.WebGPUPipes,
      y.CanvasPipes
    ],
    name: "bitmapText"
  };
  function Oe(o, e) {
    e.groupTransform = o.groupTransform, e.groupColorAlpha = o.groupColorAlpha, e.groupColor = o.groupColor, e.groupBlendMode = o.groupBlendMode, e.globalDisplayStatus = o.globalDisplayStatus, e.groupTransform = o.groupTransform, e.localDisplayStatus = o.localDisplayStatus, e.groupAlpha = o.groupAlpha, e._roundPixels = o._roundPixels;
  }
  class Mr extends $e {
    constructor(e) {
      super(), this.generatingTexture = false, this.currentKey = "--", this._renderer = e, e.runners.resolutionChange.add(this);
    }
    resolutionChange() {
      const e = this.renderable;
      e._autoResolution && e.onViewUpdate();
    }
    destroy() {
      const { htmlText: e } = this._renderer;
      e.getReferenceCount(this.currentKey) === null ? e.returnTexturePromise(this.texturePromise) : e.decreaseReferenceCount(this.currentKey), this._renderer.runners.resolutionChange.remove(this), this.texturePromise = null, this._renderer = null;
    }
  }
  function fe(o, e) {
    const { texture: t, bounds: r } = o, i = e._style._getFinalPadding();
    At(r, e._anchor, t);
    const n = e._anchor._x * i * 2, s = e._anchor._y * i * 2;
    r.minX -= i - n, r.minY -= i - s, r.maxX -= i - n, r.maxY -= i - s;
  }
  class mt {
    constructor(e) {
      this._renderer = e;
    }
    validateRenderable(e) {
      const t = this._getGpuText(e), r = e.styleKey;
      return t.currentKey !== r;
    }
    addRenderable(e, t) {
      const r = this._getGpuText(e);
      if (e._didTextUpdate) {
        const i = e._autoResolution ? this._renderer.resolution : e.resolution;
        (r.currentKey !== e.styleKey || e.resolution !== i) && this._updateGpuText(e).catch((n) => {
          console.error(n);
        }), e._didTextUpdate = false, fe(r, e);
      }
      this._renderer.renderPipes.batch.addToBatch(r, t);
    }
    updateRenderable(e) {
      const t = this._getGpuText(e);
      t._batcher.updateElement(t);
    }
    async _updateGpuText(e) {
      e._didTextUpdate = false;
      const t = this._getGpuText(e);
      if (t.generatingTexture) return;
      const r = t.texturePromise;
      t.texturePromise = null, t.generatingTexture = true, e._resolution = e._autoResolution ? this._renderer.resolution : e.resolution;
      let i = this._renderer.htmlText.getTexturePromise(e);
      r && (i = i.finally(() => {
        this._renderer.htmlText.decreaseReferenceCount(t.currentKey), this._renderer.htmlText.returnTexturePromise(r);
      })), t.texturePromise = i, t.currentKey = e.styleKey, t.texture = await i;
      const n = e.renderGroup || e.parentRenderGroup;
      n && (n.structureDidChange = true), t.generatingTexture = false, fe(t, e);
    }
    _getGpuText(e) {
      return e._gpuData[this._renderer.uid] || this.initGpuText(e);
    }
    initGpuText(e) {
      const t = new Mr(this._renderer);
      return t.renderable = e, t.transform = e.groupTransform, t.texture = R.EMPTY, t.bounds = {
        minX: 0,
        maxX: 1,
        minY: 0,
        maxY: 0
      }, t.roundPixels = this._renderer._roundPixels | e._roundPixels, e._resolution = e._autoResolution ? this._renderer.resolution : e.resolution, e._gpuData[this._renderer.uid] = t, t;
    }
    destroy() {
      this._renderer = null;
    }
  }
  mt.extension = {
    type: [
      y.WebGLPipes,
      y.WebGPUPipes,
      y.CanvasPipes
    ],
    name: "htmlText"
  };
  function Ur() {
    const { userAgent: o } = L.get().getNavigator();
    return /^((?!chrome|android).)*safari/i.test(o);
  }
  const kr = new Ee();
  function xt(o, e, t, r) {
    const i = kr;
    i.minX = 0, i.minY = 0, i.maxX = o.width / r | 0, i.maxY = o.height / r | 0;
    const n = G.getOptimalTexture(i.width, i.height, r, false);
    return n.source.uploadMethodId = "image", n.source.resource = o, n.source.alphaMode = "premultiply-alpha-on-upload", n.frame.width = e / r, n.frame.height = t / r, n.source.emit("update", n.source), n.updateUvs(), n;
  }
  function zr(o, e) {
    const t = e.fontFamily, r = [], i = {}, n = /font-family:([^;"\s]+)/g, s = o.match(n);
    function a(u) {
      i[u] || (r.push(u), i[u] = true);
    }
    if (Array.isArray(t)) for (let u = 0; u < t.length; u++) a(t[u]);
    else a(t);
    s && s.forEach((u) => {
      const l = u.split(":")[1].trim();
      a(l);
    });
    for (const u in e.tagStyles) {
      const l = e.tagStyles[u].fontFamily;
      a(l);
    }
    return r;
  }
  async function Gr(o) {
    const t = await (await L.get().fetch(o)).blob(), r = new FileReader();
    return await new Promise((n, s) => {
      r.onloadend = () => n(r.result), r.onerror = s, r.readAsDataURL(t);
    });
  }
  async function Ar(o, e) {
    const t = await Gr(e);
    return `@font-face {
        font-family: "${o.fontFamily}";
        font-weight: ${o.fontWeight};
        font-style: ${o.fontStyle};
        src: url('${t}');
    }`;
  }
  const he = /* @__PURE__ */ new Map();
  async function Dr(o) {
    const e = o.filter((t) => z.has(`${t}-and-url`)).map((t) => {
      if (!he.has(t)) {
        const { entries: r } = z.get(`${t}-and-url`), i = [];
        r.forEach((n) => {
          const s = n.url, u = n.faces.map((l) => ({
            weight: l.weight,
            style: l.style
          }));
          i.push(...u.map((l) => Ar({
            fontWeight: l.weight,
            fontStyle: l.style,
            fontFamily: t
          }, s)));
        }), he.set(t, Promise.all(i).then((n) => n.join(`
`)));
      }
      return he.get(t);
    });
    return (await Promise.all(e)).join(`
`);
  }
  function Or(o, e, t, r, i) {
    const { domElement: n, styleElement: s, svgRoot: a } = i;
    n.innerHTML = `<style>${e.cssStyle}</style><div style='padding:0;'>${o}</div>`, n.setAttribute("style", `transform: scale(${t});transform-origin: top left; display: inline-block`), s.textContent = r;
    const { width: u, height: l } = i.image;
    return a.setAttribute("width", u.toString()), a.setAttribute("height", l.toString()), new XMLSerializer().serializeToString(a);
  }
  function Ir(o, e) {
    const t = V.getOptimalCanvasAndContext(o.width, o.height, e), { context: r } = t;
    return r.clearRect(0, 0, o.width, o.height), r.drawImage(o, 0, 0), t;
  }
  function Wr(o, e, t) {
    return new Promise(async (r) => {
      t && await new Promise((i) => setTimeout(i, 100)), o.onload = () => {
        r();
      }, o.src = `data:image/svg+xml;charset=utf8,${encodeURIComponent(e)}`, o.crossOrigin = "anonymous";
    });
  }
  class _t {
    constructor(e) {
      this._activeTextures = {}, this._renderer = e, this._createCanvas = e.type === ge.WEBGPU;
    }
    getTexture(e) {
      return this.getTexturePromise(e);
    }
    getManagedTexture(e) {
      const t = e.styleKey;
      if (this._activeTextures[t]) return this._increaseReferenceCount(t), this._activeTextures[t].promise;
      const r = this._buildTexturePromise(e).then((i) => (this._activeTextures[t].texture = i, i));
      return this._activeTextures[t] = {
        texture: null,
        promise: r,
        usageCount: 1
      }, r;
    }
    getReferenceCount(e) {
      var _a;
      return ((_a = this._activeTextures[e]) == null ? void 0 : _a.usageCount) ?? null;
    }
    _increaseReferenceCount(e) {
      this._activeTextures[e].usageCount++;
    }
    decreaseReferenceCount(e) {
      const t = this._activeTextures[e];
      t && (t.usageCount--, t.usageCount === 0 && (t.texture ? this._cleanUp(t.texture) : t.promise.then((r) => {
        t.texture = r, this._cleanUp(t.texture);
      }).catch(() => {
        X("HTMLTextSystem: Failed to clean texture");
      }), this._activeTextures[e] = null));
    }
    getTexturePromise(e) {
      return this._buildTexturePromise(e);
    }
    async _buildTexturePromise(e) {
      const { text: t, style: r, resolution: i, textureStyle: n } = e, s = ee.get(it), a = zr(t, r), u = await Dr(a), l = ir(t, r, u, s), h = Math.ceil(Math.ceil(Math.max(1, l.width) + r.padding * 2) * i), c = Math.ceil(Math.ceil(Math.max(1, l.height) + r.padding * 2) * i), d = s.image, f = 2;
      d.width = (h | 0) + f, d.height = (c | 0) + f;
      const m = Or(t, r, i, u, s);
      await Wr(d, m, Ur() && a.length > 0);
      const v = d;
      let x;
      this._createCanvas && (x = Ir(d, i));
      const p = xt(x ? x.canvas : v, d.width - f, d.height - f, i);
      return n && (p.source.style = n), this._createCanvas && (this._renderer.texture.initSource(p.source), V.returnCanvasAndContext(x)), ee.return(s), p;
    }
    returnTexturePromise(e) {
      e.then((t) => {
        this._cleanUp(t);
      }).catch(() => {
        X("HTMLTextSystem: Failed to clean texture");
      });
    }
    _cleanUp(e) {
      G.returnTexture(e, true), e.source.resource = null, e.source.uploadMethodId = "unknown";
    }
    destroy() {
      this._renderer = null;
      for (const e in this._activeTextures) this._activeTextures[e] && this.returnTexturePromise(this._activeTextures[e].promise);
      this._activeTextures = null;
    }
  }
  _t.extension = {
    type: [
      y.WebGLSystem,
      y.WebGPUSystem,
      y.CanvasSystem
    ],
    name: "htmlText"
  };
  class Lr extends $e {
    constructor(e) {
      super(), this._renderer = e, e.runners.resolutionChange.add(this);
    }
    resolutionChange() {
      const e = this.renderable;
      e._autoResolution && e.onViewUpdate();
    }
    destroy() {
      const { canvasText: e } = this._renderer;
      e.getReferenceCount(this.currentKey) > 0 ? e.decreaseReferenceCount(this.currentKey) : this.texture && e.returnTexture(this.texture), this._renderer.runners.resolutionChange.remove(this), this._renderer = null;
    }
  }
  class vt {
    constructor(e) {
      this._renderer = e;
    }
    validateRenderable(e) {
      const t = this._getGpuText(e), r = e.styleKey;
      return t.currentKey !== r ? true : e._didTextUpdate;
    }
    addRenderable(e, t) {
      const r = this._getGpuText(e);
      if (e._didTextUpdate) {
        const i = e._autoResolution ? this._renderer.resolution : e.resolution;
        (r.currentKey !== e.styleKey || e.resolution !== i) && this._updateGpuText(e), e._didTextUpdate = false, fe(r, e);
      }
      this._renderer.renderPipes.batch.addToBatch(r, t);
    }
    updateRenderable(e) {
      const t = this._getGpuText(e);
      t._batcher.updateElement(t);
    }
    _updateGpuText(e) {
      const t = this._getGpuText(e);
      t.texture && this._renderer.canvasText.decreaseReferenceCount(t.currentKey), e._resolution = e._autoResolution ? this._renderer.resolution : e.resolution, t.texture = this._renderer.canvasText.getManagedTexture(e), t.currentKey = e.styleKey;
    }
    _getGpuText(e) {
      return e._gpuData[this._renderer.uid] || this.initGpuText(e);
    }
    initGpuText(e) {
      const t = new Lr(this._renderer);
      return t.currentKey = "--", t.renderable = e, t.transform = e.groupTransform, t.bounds = {
        minX: 0,
        maxX: 1,
        minY: 0,
        maxY: 0
      }, t.roundPixels = this._renderer._roundPixels | e._roundPixels, e._gpuData[this._renderer.uid] = t, t;
    }
    destroy() {
      this._renderer = null;
    }
  }
  vt.extension = {
    type: [
      y.WebGLPipes,
      y.WebGPUPipes,
      y.CanvasPipes
    ],
    name: "text"
  };
  class bt {
    constructor(e) {
      this._activeTextures = {}, this._renderer = e;
    }
    getTexture(e, t, r, i) {
      typeof e == "string" && (U("8.0.0", "CanvasTextSystem.getTexture: Use object TextOptions instead of separate arguments"), e = {
        text: e,
        style: r,
        resolution: t
      }), e.style instanceof N || (e.style = new N(e.style)), e.textureStyle instanceof $ || (e.textureStyle = new $(e.textureStyle)), typeof e.text != "string" && (e.text = e.text.toString());
      const { text: n, style: s, textureStyle: a } = e, u = e.resolution ?? this._renderer.resolution, { frame: l, canvasAndContext: h } = ae.getCanvasAndContext({
        text: n,
        style: s,
        resolution: u
      }), c = xt(h.canvas, l.width, l.height, u);
      if (a && (c.source.style = a), s.trim && (l.pad(s.padding), c.frame.copyFrom(l), c.frame.scale(1 / u), c.updateUvs()), s.filters) {
        const d = this._applyFilters(c, s.filters);
        return this.returnTexture(c), ae.returnCanvasAndContext(h), d;
      }
      return this._renderer.texture.initSource(c._source), ae.returnCanvasAndContext(h), c;
    }
    returnTexture(e) {
      const t = e.source;
      t.resource = null, t.uploadMethodId = "unknown", t.alphaMode = "no-premultiply-alpha", G.returnTexture(e, true);
    }
    renderTextToCanvas() {
      U("8.10.0", "CanvasTextSystem.renderTextToCanvas: no longer supported, use CanvasTextSystem.getTexture instead");
    }
    getManagedTexture(e) {
      e._resolution = e._autoResolution ? this._renderer.resolution : e.resolution;
      const t = e.styleKey;
      if (this._activeTextures[t]) return this._increaseReferenceCount(t), this._activeTextures[t].texture;
      const r = this.getTexture({
        text: e.text,
        style: e.style,
        resolution: e._resolution,
        textureStyle: e.textureStyle
      });
      return this._activeTextures[t] = {
        texture: r,
        usageCount: 1
      }, r;
    }
    decreaseReferenceCount(e) {
      const t = this._activeTextures[e];
      t.usageCount--, t.usageCount === 0 && (this.returnTexture(t.texture), this._activeTextures[e] = null);
    }
    getReferenceCount(e) {
      var _a;
      return ((_a = this._activeTextures[e]) == null ? void 0 : _a.usageCount) ?? 0;
    }
    _increaseReferenceCount(e) {
      this._activeTextures[e].usageCount++;
    }
    _applyFilters(e, t) {
      const r = this._renderer.renderTarget.renderTarget, i = this._renderer.filter.generateFilteredTexture({
        texture: e,
        filters: t
      });
      return this._renderer.renderTarget.bind(r, false), i;
    }
    destroy() {
      this._renderer = null;
      for (const e in this._activeTextures) this._activeTextures[e] && this.returnTexture(this._activeTextures[e].texture);
      this._activeTextures = null;
    }
  }
  bt.extension = {
    type: [
      y.WebGLSystem,
      y.WebGPUSystem,
      y.CanvasSystem
    ],
    name: "canvasText"
  };
  C.add(Ne);
  C.add(Ke);
  C.add(nt);
  C.add(Dt);
  C.add(ot);
  C.add(lt);
  C.add(ct);
  C.add(bt);
  C.add(vt);
  C.add(gt);
  C.add(_t);
  C.add(mt);
  C.add(pt);
  C.add(ft);
  C.add(et);
  C.add(Ze);
});
export {
  __tla
};
