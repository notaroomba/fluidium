import { D as n, n as c, G as r, __tla as __tla_0 } from "./index-PYdGGKf0.js";
let v;
let __tla = Promise.all([
  (() => {
    try {
      return __tla_0;
    } catch {
    }
  })()
]).then(async () => {
  class l {
    constructor(a) {
      this._canvasPool = /* @__PURE__ */ Object.create(null), this.canvasOptions = a || {}, this.enableFullScreen = false;
    }
    _createCanvasAndContext(a, s) {
      const t = n.get().createCanvas();
      t.width = a, t.height = s;
      const e = t.getContext("2d");
      return {
        canvas: t,
        context: e
      };
    }
    getOptimalCanvasAndContext(a, s, t = 1) {
      a = Math.ceil(a * t - 1e-6), s = Math.ceil(s * t - 1e-6), a = c(a), s = c(s);
      const e = (a << 17) + (s << 1);
      this._canvasPool[e] || (this._canvasPool[e] = []);
      let o = this._canvasPool[e].pop();
      return o || (o = this._createCanvasAndContext(a, s)), o;
    }
    returnCanvasAndContext(a) {
      const s = a.canvas, { width: t, height: e } = s, o = (t << 17) + (e << 1);
      a.context.resetTransform(), a.context.clearRect(0, 0, t, e), this._canvasPool[o].push(a);
    }
    clear() {
      this._canvasPool = {};
    }
  }
  v = new l();
  r.register(v);
});
export {
  v as C,
  __tla
};
