<script setup lang="ts">
// The home page of the docs. The animation plays the steady cycle of the
// four-step VSA of Haghpanah et al. (2013), written by docs/home_cycle_data.jl.
import { ref, computed, onMounted, onBeforeUnmount } from 'vue'
import { withBase } from 'vitepress'
import cycle from '../assets/home_cycle.json'

const link = (path: string) => withBase(path)

// ---------------------------------------------------------------- steps
type Step = 'pressurise' | 'adsorb' | 'blowdown' | 'evacuate' | 'equalise' | 'lpp' | 'idle'
const STEP_LABEL: Record<Step, string> = {
  pressurise: 'Pressurisation', adsorb: 'Adsorption', blowdown: 'Blowdown', evacuate: 'Evacuation',
  equalise: 'Equalisation', lpp: 'Light product pressurisation', idle: 'Idle',
}
function stepOf(name: string): Step {
  const n = name.toLowerCase()
  if (n.includes('blow')) return 'blowdown'
  if (n.includes('evac')) return 'evacuate'
  if (n.includes('equal')) return 'equalise'
  if (n.includes('ads')) return 'adsorb'
  if (n.includes('press')) return 'pressurise'
  return 'idle'
}

// ------------------------------------------------------------- the cycle
const T = cycle.cycle_time_s
const stages = (() => {
  let start = 0
  return cycle.stages.map((s: { name: string, duration: number }) => {
    const st = { step: stepOf(s.name), start, end: start + s.duration }
    start += s.duration
    return st
  })
})()
const stageAt = (t: number) => Math.max(0, stages.findIndex(s => t < s.end))

const nf = cycle.t.length
// Frames either side of time t, and the weight of the second
function bracket(t: number): [number, number, number] {
  const ts = cycle.t
  if (t < ts[0] || t >= ts[nf - 1]) {
    const x = t < ts[0] ? t + T : t
    return [nf - 1, 0, (x - ts[nf - 1]) / (ts[0] + T - ts[nf - 1])]
  }
  let lo = 0, hi = nf - 1
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (ts[m] <= t) lo = m; else hi = m }
  return [lo, hi, (t - ts[lo]) / (ts[hi] - ts[lo])]
}
const lerp = (a: number, b: number, w: number) => a + (b - a) * w

const qMin = Math.min(...cycle.loading.flat())
const PMAX = 1.1

// ------------------------------------------------------------- geometry
const BED = { x: 40, w: 70, top: 44, bottom: 324 }
const BED_H = BED.bottom - BED.top
const PLOT = { x0: 252, x1: 520, y0: 250, y1: 70 }
const px = (t: number) => PLOT.x0 + (t / T) * (PLOT.x1 - PLOT.x0)
const py = (p: number) => PLOT.y0 - (p / PMAX) * (PLOT.y0 - PLOT.y1)

const nz = cycle.z.length
const cells = cycle.z.map((_: number, i: number) => ({
  y: BED.bottom - (i + 1) * BED_H / nz,
  h: BED_H / nz + 0.4,
}))

// Particles, in rows of six and five, each coloured by the loading of its cell
const beads: { x: number, y: number, i: number }[] = []
const ROWS = 26
for (let r = 0; r < ROWS; r++) {
  const y = BED.bottom - 9 - r * (BED_H - 18) / (ROWS - 1)
  const z = (BED.bottom - y) / BED_H
  const i = Math.min(nz - 1, Math.floor(z * nz))
  const odd = r % 2 === 1
  for (let c = 0; c < (odd ? 5 : 6); c++) beads.push({ x: BED.x + 10 + (odd ? 5 : 0) + c * 10, y, i })
}

// The pressure set at the open valve, from the start of the cycle, where it
// carries on from the end of the evacuation
const pressurePath = [0, ...cycle.t]
  .map((t: number, k: number) => {
    const p = k ? cycle.pressure_valve_bar[k - 1] : cycle.pressure_valve_bar[nf - 1]
    return `${k ? 'L' : 'M'}${px(t).toFixed(1)} ${py(p).toFixed(1)}`
  })
  .join(' ')
const timeTicks = Array.from({ length: Math.floor(T / 25) + 1 }, (_, k) => 25 * k)
const pressureTicks = [0, 0.5, 1.0]

// ------------------------------------------------------------- animation
const SPEEDUP = 7
const time = ref(stages[1] ? (stages[1].start + stages[1].end) / 2 : 0)  // mid-adsorption at rest
const playing = ref(false)
let raf = 0
let last = 0

const frame = computed(() => {
  const [a, b, w] = bracket(time.value)
  const q = cycle.loading[a].map((v: number, i: number) => lerp(v, cycle.loading[b][i], w))
  const y = cycle.y_co2[a].map((v: number, i: number) => lerp(v, cycle.y_co2[b][i], w))
  const p = lerp(cycle.pressure_valve_bar[a], cycle.pressure_valve_bar[b], w)
  return { q, y, p }
})
const stage = computed(() => stages[stageAt(time.value)])
const step = computed(() => stage.value.step)
const beadOpacity = (i: number) => (0.06 + 0.9 * (frame.value.q[i] - qMin) / (1 - qMin)).toFixed(3)
const gasOpacity = (i: number) => (0.6 * Math.sqrt(Math.max(0, frame.value.y[i]))).toFixed(3)

const feedIn = computed(() => step.value === 'pressurise' || step.value === 'adsorb')
const topOut = computed(() => step.value === 'adsorb' || step.value === 'blowdown')
const bottomOut = computed(() => step.value === 'evacuate')

function tick(now: number) {
  if (last) time.value = (time.value + (now - last) / 1000 * SPEEDUP) % T
  last = now
  raf = requestAnimationFrame(tick)
}
function play() { playing.value = true; last = 0; raf = requestAnimationFrame(tick) }
function pause() { playing.value = false; cancelAnimationFrame(raf) }
const toggle = () => (playing.value ? pause() : play())

onMounted(() => {
  if (!window.matchMedia('(prefers-reduced-motion: reduce)').matches) play()
})
onBeforeUnmount(pause)

// ------------------------------------------------------------- content
const pct = (x: number) => (100 * x).toFixed(1)

// Kinds of capture unit. Icons are drawn in a 64 × 56 box with currentColor.
type Status = 'available' | 'development' | 'future'
const STATUS_LABEL: Record<Status, string> = { available: 'Available', development: 'In development', future: 'Future' }
const units: { title: string, status: Status, text: string, href?: string, icon: string }[] = [
  {
    title: 'Fixed-bed adsorption', status: 'available', href: '/model',
    text: 'Packed columns cycled through pressure and vacuum swing steps: PSA, VSA, equalisation, several beds, layered adsorbents.',
    icon: '<rect x="22" y="4" width="20" height="48" rx="5" fill="none" stroke="currentColor" stroke-width="1.6"/>'
      + [10, 17, 24, 31, 38, 45].map(y => [27, 32, 37].map(x => `<circle cx="${x}" cy="${y}" r="2.2" fill="currentColor"/>`).join('')).join(''),
  },
  {
    title: 'Moving-bed adsorption', status: 'development',
    text: 'Sorbent that flows through adsorption and regeneration sections and back, such as moving-bed temperature swing adsorption.',
    icon: '<rect x="14" y="4" width="16" height="48" rx="4" fill="none" stroke="currentColor" stroke-width="1.6"/>'
      + [12, 22, 32, 42].map(y => `<circle cx="22" cy="${y}" r="2.2" fill="currentColor"/><path d="M22 ${y + 3} v3" stroke="currentColor" stroke-width="1.2"/>`).join('')
      + '<path d="M30 48 h14 a6 6 0 0 0 6 -6 v-28 a6 6 0 0 0 -6 -6 h-14" fill="none" stroke="currentColor" stroke-width="1.4" stroke-dasharray="3 2"/>',
  },
  {
    title: 'Rotating-bed adsorption', status: 'future',
    text: 'Adsorbent wheels that turn through adsorption and regeneration sectors.',
    icon: '<circle cx="32" cy="28" r="22" fill="none" stroke="currentColor" stroke-width="1.6"/>'
      + '<path d="M32 28 L32 6 M32 28 L51 39 M32 28 L13 39" stroke="currentColor" stroke-width="1.4"/>'
      + '<path d="M32 6 A22 22 0 0 1 51 39 L32 28 Z" fill="currentColor" opacity=".18"/>',
  },
  {
    title: 'Membranes', status: 'future',
    text: 'Membrane modules that let CO₂ through faster than N₂, alone or in stages.',
    icon: '<rect x="6" y="10" width="52" height="36" rx="3" fill="none" stroke="currentColor" stroke-width="1.6"/>'
      + '<path d="M6 28 H58" stroke="currentColor" stroke-width="1.4" stroke-dasharray="4 3"/>'
      + '<path d="M12 20 H46" stroke="currentColor" stroke-width="1.2" marker-end="url(#mh-arrow)"/>'
      + '<path d="M26 22 V36" stroke="currentColor" stroke-width="1.2" marker-end="url(#mh-arrow)"/>',
  },
  {
    title: 'Absorber columns', status: 'future',
    text: 'Gas–liquid contactors for solvent-based capture, with packing or trays.',
    icon: '<rect x="22" y="4" width="20" height="48" rx="5" fill="none" stroke="currentColor" stroke-width="1.6"/>'
      + [14, 22, 30, 38].map((y, k) => `<path d="M${k % 2 ? 28 : 22} ${y} h14" stroke="currentColor" stroke-width="1.2"/>`).join('')
      + '<path d="M14 48 h8 M42 8 h8" stroke="currentColor" stroke-width="1.4"/>',
  },
]

type Row = { label: string, steps: [Step, number][] }
const library: { title: string, cite: string, text: string, href: string, rows: Row[] }[] = [
  {
    title: 'Four-step VSA', cite: 'Haghpanah et al. 2013 · zeolite 13X',
    text: 'One bed, dry flue gas. The reference case of these docs.',
    href: '/examples/getting_started',
    rows: [{ label: 'Bed', steps: [['pressurise', 15], ['adsorb', 15], ['blowdown', 30], ['evacuate', 40]] }],
  },
  {
    title: 'Two-bed VSA with equalisation', cite: 'Haghpanah beds · two_bed_vsa_process()',
    text: 'The beds pass gas to each other through their tops instead of sending it to the vacuum pump.',
    href: '/examples/two_bed_vsa',
    rows: [
      { label: 'A', steps: [['pressurise', 15], ['adsorb', 15], ['equalise', 10], ['evacuate', 30], ['equalise', 10]] },
      { label: 'B', steps: [['evacuate', 30], ['equalise', 10], ['pressurise', 15], ['adsorb', 15], ['equalise', 10]] },
    ],
  },
  {
    title: 'Wet flue gas, LPP cycle', cite: 'Krishnamurthy et al. 2014 · CO₂, N₂, H₂O',
    text: 'Water builds up over hundreds of cycles; Newton\'s method still finds the steady state.',
    href: '/examples/wet_flue_gas_lpp_vsa',
    rows: [{ label: 'Bed', steps: [['lpp', 20], ['adsorb', 78.4], ['blowdown', 37.9], ['evacuate', 111.6]] }],
  },
  {
    title: 'Wet flue gas, dual adsorbent', cite: 'Krishnamurthy et al. 2014 · silica gel + 13X',
    text: 'A silica gel bed takes out the water before the zeolite bed captures the CO₂.',
    href: '/examples/wet_flue_gas_dual_adsorbent_vsa',
    rows: [
      { label: 'SG', steps: [['pressurise', 20], ['adsorb', 46.2], ['blowdown', 42.21], ['evacuate', 46.02], ['idle', 69.27]] },
      { label: '13X', steps: [['lpp', 20], ['adsorb', 46.2], ['blowdown', 56.3], ['evacuate', 101.2]] },
    ],
  },
]
const legendSteps: Step[] = ['pressurise', 'adsorb', 'blowdown', 'evacuate', 'equalise', 'lpp']
const scheduleLabel = (rows: Row[]) =>
  rows.map(r => `${r.label}: ` + r.steps.map(([s, d]) => `${STEP_LABEL[s].toLowerCase()} ${d} s`).join(', ')).join('; ')

const examples = [
  { tag: 'Basics', title: 'Building a process from its parts', text: 'Your own gas, adsorbent, column and cycle, from scratch.', href: '/examples/process_from_parts' },
  { tag: 'Basics', title: 'Breakthrough curve', text: 'A single adsorption step, and what the heat of adsorption does to it.', href: '/examples/breakthrough' },
  { tag: 'Flowsheets', title: 'Two-bed VSA with pressure equalisation', text: 'Two beds working in turn, and what the equalisation is worth.', href: '/examples/two_bed_vsa' },
  { tag: 'Literature', title: 'Wet flue gas VSA', text: 'The processes of Krishnamurthy et al. (2014), against the paper.', href: '/examples/wet_flue_gas_lpp_vsa' },
  { tag: 'Design', title: 'Gradients at cyclic steady state', text: 'How each setting moves purity and recovery, checked against finite differences.', href: '/examples/flowsheet_vsa_gradients' },
  { tag: 'Design', title: 'Optimization', text: 'Gradient-based optimisation of a cyclic process.', href: '/examples/optimization' },
]
</script>

<template>
  <div class="mh">
    <!-- Hero -->
    <section class="hero">
      <div class="hero-text">
        <p class="eyebrow">Open source · Differentiable · Julia</p>
        <h1>CO₂ capture processes, from flowsheet to <em>optimum</em>.</h1>
        <p class="lede">
          Set up a capture process in a few lines, from a published case or from its parts.
          Simulate it in seconds, and let exact gradients steer purity, recovery and energy to the
          best operating point.
        </p>
        <p class="lede-note">
          Today's units are fixed beds of adsorbent, for pressure and vacuum swing adsorption.
          The framework is built for more.
        </p>
        <div class="actions">
          <a class="btn primary" :href="link('/examples/getting_started')">Get started</a>
          <a class="btn" :href="link('/installation')">Install</a>
          <a class="btn" :href="link('/examples/optimization')">Optimise a process</a>
        </div>
      </div>

      <figure class="cycle">
        <svg viewBox="0 0 540 400" role="img"
          aria-label="Animation of the packed bed of the four-step VSA of Haghpanah et al. (2013) through one steady cycle, with the pressure at the open valve">
          <defs>
            <marker id="mh-arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
              <path class="arrowhead" d="M0,0 L10,5 L0,10 z" />
            </marker>
          </defs>

          <!-- the bed -->
          <rect v-for="(c, i) in cells" :key="'g' + i" class="gas" :x="BED.x + 1" :y="c.y" :width="BED.w - 2" :height="c.h" :fill-opacity="gasOpacity(i)" />
          <g v-for="(b, k) in beads" :key="'b' + k">
            <circle class="bead-base" :cx="b.x" :cy="b.y" r="4.1" />
            <circle class="bead" :cx="b.x" :cy="b.y" r="4.1" :fill-opacity="beadOpacity(b.i)" />
          </g>
          <rect class="bed-outline" :x="BED.x" :y="BED.top" :width="BED.w" :height="BED_H" rx="12" />
          <text class="small" :x="BED.x + BED.w + 8" :y="BED.top + 4">top</text>
          <text class="small" :x="BED.x + BED.w + 8" :y="BED.bottom + 2">bottom</text>

          <!-- streams -->
          <g :style="{ opacity: topOut ? 1 : 0 }" class="stream">
            <line class="arrow" x1="75" y1="40" x2="75" y2="8" marker-end="url(#mh-arrow)" />
            <text class="label" x="88" y="24">{{ step === 'blowdown' ? 'Blowdown gas' : 'N₂-rich product' }}</text>
          </g>
          <g :style="{ opacity: feedIn ? 1 : 0 }" class="stream">
            <line class="arrow" x1="75" y1="394" x2="75" y2="330" marker-end="url(#mh-arrow)" />
            <text class="label" x="88" y="370">Feed, 15% CO₂</text>
          </g>
          <g :style="{ opacity: bottomOut ? 1 : 0 }" class="stream">
            <line class="arrow" x1="75" y1="330" x2="75" y2="394" marker-end="url(#mh-arrow)" />
            <text class="label" x="88" y="370">CO₂ to the vacuum pump</text>
          </g>

          <!-- pressure at the open valve -->
          <text class="step-name" :x="PLOT.x0" y="40">{{ STEP_LABEL[step] }}</text>
          <line class="axis" :x1="PLOT.x0" :y1="PLOT.y0" :x2="PLOT.x1" :y2="PLOT.y0" />
          <line class="axis" :x1="PLOT.x0" :y1="PLOT.y1 - 6" :x2="PLOT.x0" :y2="PLOT.y0" />
          <g v-for="p in pressureTicks" :key="'p' + p">
            <line class="grid" :x1="PLOT.x0" :y1="py(p)" :x2="PLOT.x1" :y2="py(p)" />
            <text class="small" :x="PLOT.x0 - 6" :y="py(p) + 4" text-anchor="end">{{ p.toFixed(1) }}</text>
          </g>
          <text class="small" :x="PLOT.x0" :y="PLOT.y1 - 14">pressure at the open valve [bar]</text>
          <path class="pressure" :d="pressurePath" />
          <rect v-for="(s, k) in stages" :key="'s' + k" :class="'st-' + s.step"
            :x="px(s.start) + 0.5" y="262" :width="px(s.end) - px(s.start) - 1" height="12"
            :opacity="s === stage ? 1 : 0.35" />
          <text v-for="t in timeTicks" :key="'t' + t" class="small" :x="px(t)" y="292" text-anchor="middle">{{ t }}</text>
          <text class="small" :x="PLOT.x1" y="310" text-anchor="end">time in the cycle [s]</text>
          <line class="cursor" :x1="px(time)" :y1="PLOT.y1 - 6" :x2="px(time)" y2="276" />
          <circle class="dot" :cx="px(time)" :cy="py(frame.p)" r="4.5" />
          <text class="readout" :x="PLOT.x1" y="40" text-anchor="end">{{ time.toFixed(0) }} s · {{ frame.p.toFixed(2) }} bar</text>

          <!-- key -->
          <g transform="translate(252 336)">
            <circle class="bead" cx="6" cy="6" r="4.1" fill-opacity="0.08" />
            <circle class="bead" cx="18" cy="6" r="4.1" fill-opacity="0.5" />
            <circle class="bead" cx="30" cy="6" r="4.1" fill-opacity="0.95" />
            <text class="small" x="42" y="10">CO₂ loading, low to high</text>
            <rect class="gas" x="0" y="22" width="36" height="10" fill-opacity="0.55" />
            <text class="small" x="42" y="31">CO₂ in the gas</text>
          </g>
        </svg>
        <figcaption>
          <span class="legend">
            <span v-for="s in stages" :key="'l' + s.start" class="chip" :class="{ on: s === stage }">
              <i :class="'bg-' + s.step"></i>{{ STEP_LABEL[s.step] }}
            </span>
          </span>
          <button class="toggle" type="button" @click="toggle">{{ playing ? 'Pause' : 'Play' }}</button>
          <span class="source">
            A fixed-bed unit in the four-step VSA cycle of Haghpanah et al. (2013), as Mocca
            computes it, at {{ SPEEDUP }}× speed.
          </span>
        </figcaption>
      </figure>
    </section>

    <!-- Set up, simulate, optimise -->
    <section class="band">
      <header>
        <h2>Set up. Simulate. Optimise.</h2>
        <p>The four-step VSA of Haghpanah et al. (2013), from Mocca's library to an optimised design.</p>
      </header>
      <div class="pillars">
        <div class="pillar">
          <span class="pillar-n">1</span>
          <h3>Set up</h3>
          <p>Read a process from the library, or build it from its gas, adsorbent, column, flowsheet and cycle, in Julia or a JSON file.</p>
          <pre><code><span class="k">using</span> Mocca, Mocca.Units
process = read_process(
  <span class="s">"haghpanah_2013_four_stage_vsa"</span>)
process[:P_L] = 0.05bar</code></pre>
        </div>
        <div class="pillar">
          <span class="pillar-n">2</span>
          <h3>Simulate</h3>
          <p>Run the cycle to its steady state and measure it, with plots, stream tables and spreadsheets.</p>
          <pre><code>result = cyclic_steady_state(process)
perf = performance(result)</code></pre>
          <div class="nums">
            <div class="num"><span class="val"><b>{{ pct(cycle.purity) }}</b> <small>%</small></span><span>CO₂ purity</span></div>
            <div class="num"><span class="val"><b>{{ pct(cycle.recovery) }}</b> <small>%</small></span><span>recovery</span></div>
            <div class="num"><span class="val"><b>{{ cycle.specific_energy_kwh_per_t.toFixed(0) }}</b> <small>kWh/t</small></span><span>energy</span></div>
          </div>
        </div>
        <div class="pillar">
          <span class="pillar-n">3</span>
          <h3>Optimise</h3>
          <p>Every result is differentiable. Adjoint gradients with respect to any setting drive a gradient-based optimiser.</p>
          <pre><code>add_variable!(process, :P_L,
  (<span class="s">"evacuation"</span>, :to))
optimise!(process;
  maximise = :recovery,
  at_least = (purity = 0.9,))</code></pre>
        </div>
      </div>
      <p class="aside">The first run of a session compiles for a few minutes; after that a simulation takes seconds.</p>
    </section>

    <!-- Capture units -->
    <section id="units" class="band">
      <header>
        <h2>Capture units</h2>
        <p>
          Each kind of unit is its own module on a shared core. The flowsheet, the solvers, the
          performance measures and the gradients work the same for every unit, so a new kind of
          unit gets all of them.
        </p>
      </header>
      <div class="units">
        <component :is="u.href ? 'a' : 'div'" v-for="u in units" :key="u.title" class="unit" :class="'is-' + u.status" :href="u.href ? link(u.href) : undefined">
          <svg class="unit-icon" viewBox="0 0 64 56" aria-hidden="true" v-html="u.icon"></svg>
          <div class="unit-text">
            <span class="status">{{ STATUS_LABEL[u.status] }}</span>
            <h3>{{ u.title }}</h3>
            <p>{{ u.text }}</p>
          </div>
        </component>
      </div>
      <p class="aside">
        A new unit implements a small set of functions: its ports, its balances and how it
        exchanges streams with its neighbours. See <a :href="link('/extending')">Extending Mocca</a>.
      </p>
    </section>

    <!-- Library -->
    <section class="band">
      <header>
        <h2>Processes from the literature, ready to run</h2>
        <p>The library holds fixed-bed adsorption processes so far. Read any of them with <code>read_process</code>, change a setting and run it. Bars show each bed's steps over one cycle, to scale. <a class="inline" :href="link('/processes')">How Mocca's results compare with the papers</a></p>
      </header>
      <div class="cards">
        <a v-for="p in library" :key="p.title" class="card" :href="link(p.href)">
          <h3>{{ p.title }}</h3>
          <span class="cite">{{ p.cite }}</span>
          <div class="sched" role="img" :aria-label="scheduleLabel(p.rows)">
            <div v-for="r in p.rows" :key="r.label" class="sched-row">
              <span>{{ r.label }}</span>
              <span class="sched-bar">
                <i v-for="([s, d], k) in r.steps" :key="k" :class="'bg-' + s" :style="{ flex: d }"></i>
              </span>
            </div>
          </div>
          <p>{{ p.text }}</p>
        </a>
      </div>
      <div class="legend">
        <span v-for="s in legendSteps" :key="s" class="chip"><i :class="'bg-' + s"></i>{{ STEP_LABEL[s] }}</span>
        <span class="chip"><i class="bg-idle"></i>Idle</span>
      </div>
    </section>

    <!-- What Mocca models -->
    <section class="band">
      <header><h2>What Mocca models</h2></header>
      <div class="sheet-wrap">
        <table class="sheet">
          <tbody>
            <tr class="group"><th colspan="2">For every process</th></tr>
            <tr><th>Flowsheets</th><td>Any number of units, connected to feeds, products and each other through valves, flow controllers and streams.</td></tr>
            <tr><th>Design</th><td>Adjoint gradients of purity, recovery, energy and productivity with respect to step durations, pressures, flow rates and model parameters; gradient-based optimisation with bounds and constraints.</td></tr>
            <tr><th>Operation</th><td>Cycles of steps run in sequence or side by side in several units, solved to their cyclic steady state by Newton's method, with GMRES and grid sequencing on fine grids.</td></tr>
            <tr><th>Results</th><td>Stream tables, profiles and histories as CSV files in engineering units (bara, °C, mol%, kg/h). <a :href="link('/results')">Results and plots</a></td></tr>
            <tr class="group"><th colspan="2">Fixed-bed adsorption unit</th></tr>
            <tr><th>Bed</th><td>One-dimensional, axially dispersed plug flow through a packed column, Darcy's law for the pressure drop, a linear driving force for mass transfer, and an energy balance with the column wall (or an adiabatic or isothermal bed). <a :href="link('/model')">Model equations</a></td></tr>
            <tr><th>Isotherms</th><td>Multicomponent dual-site and single-site Langmuir, with a library of published adsorbents (zeolite 13X, silica gel).</td></tr>
            <tr><th>Steps</th><td>Pressurisation (with feed or light product), adsorption, blowdown, evacuation, purge, pressure equalisation and idle.</td></tr>
            <tr class="not"><th>Not yet</th><td>Temperature swing (heated or cooled steps), the inertial (Ergun) term of the pressure drop, other isotherm models (add your own with <code>register_isotherm!</code>).</td></tr>
          </tbody>
        </table>
      </div>
    </section>

    <!-- Examples -->
    <section class="band">
      <header>
        <h2>Examples</h2>
        <p>Each example is a script you can download and run, or open as a notebook.</p>
      </header>
      <div class="examples">
        <a v-for="e in examples" :key="e.href" class="ex" :href="link(e.href)">
          <span class="tag">{{ e.tag }}</span>
          <h3>{{ e.title }}</h3>
          <p>{{ e.text }}</p>
        </a>
      </div>
    </section>

    <!-- Where next -->
    <section class="band last">
      <header><h2>Where next</h2></header>
      <div class="next">
        <div>
          <h3>Learn</h3>
          <a :href="link('/installation')">Install Julia and Mocca</a>
          <a :href="link('/vsa_cycle')">How a VSA cycle works</a>
          <a :href="link('/examples/getting_started')">Run, change, export and save a process</a>
          <a :href="link('/model')">The equations, and how performance is defined</a>
        </div>
        <div>
          <h3>Build</h3>
          <a :href="link('/examples/process_from_parts')">Your own gas, adsorbent, column and cycle</a>
          <a :href="link('/examples/two_bed_vsa')">Connect several beds</a>
          <a :href="link('/examples/optimization')">Optimise a process</a>
          <a :href="link('/extending')">Add a new kind of unit or step</a>
        </div>
        <div>
          <h3>Look up</h3>
          <a :href="link('/user_interface')">Functions</a>
          <a :href="link('/json_input')">The process file format</a>
          <a :href="link('/cycle_steps')">The steps of a cycle</a>
          <a :href="link('/processes')">The library of processes, against the papers</a>
          <a :href="link('/adsorbents')">The library of adsorbents</a>
          <a :href="link('/results')">Numbers, plots and spreadsheets</a>
        </div>
      </div>
    </section>
  </div>
</template>

<style scoped>
.mh {
  max-width: 1152px;
  margin: 0 auto;
  padding: 0 24px;
  color: var(--vp-c-text-1);
}
.mh a { color: inherit; text-decoration: none; }
.mh code { font-size: 0.9em; }

/* ---------- hero ---------- */
.hero {
  display: grid;
  grid-template-columns: minmax(0, 1fr) minmax(0, 1.1fr);
  gap: 40px;
  align-items: center;
  padding: 56px 0 44px;
  border-bottom: 1px solid var(--vp-c-divider);
}
.eyebrow {
  margin: 0 0 18px;
  font: 500 12px/1 var(--vp-font-family-mono);
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--vp-c-text-2);
}
.hero h1 {
  margin: 0;
  font-family: var(--mc-font-display);
  font-stretch: 122%;
  font-weight: 800;
  font-size: clamp(34px, 4.6vw, 54px);
  line-height: 1.04;
  letter-spacing: -0.012em;
  text-wrap: balance;
}
.hero h1 em { font-style: normal; color: var(--mc-co2); }
.lede { margin: 18px 0 0; max-width: 44ch; font-size: 17px; line-height: 1.6; color: var(--vp-c-text-2); }
.lede-note { margin: 12px 0 0; max-width: 46ch; font-size: 15px; line-height: 1.6; color: var(--vp-c-text-2); padding-left: 12px; border-left: 2px solid var(--mc-co2); }
.actions { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 28px; }
.btn {
  font: 600 14px/1 var(--vp-font-family-base);
  padding: 12px 18px;
  border-radius: 4px;
  border: 1px solid var(--vp-c-text-1);
  transition: background-color 0.15s, color 0.15s;
}
.btn:hover { background: var(--vp-c-bg-soft); }
.btn.primary { background: var(--vp-c-brand-3); border-color: var(--vp-c-brand-3); color: #fff; }
.btn.primary:hover { background: var(--vp-c-brand-2); border-color: var(--vp-c-brand-2); }

.cycle { margin: 0; display: grid; gap: 10px; min-width: 0; }
.cycle svg { width: 100%; height: auto; display: block; }
.cycle figcaption { display: flex; flex-wrap: wrap; align-items: center; gap: 8px 14px; font-size: 13px; color: var(--vp-c-text-2); }
.source { flex-basis: 100%; }
.toggle {
  margin-left: auto;
  font: 500 12px var(--vp-font-family-mono);
  color: var(--vp-c-text-2);
  border: 1px solid var(--vp-c-divider);
  border-radius: 4px;
  padding: 3px 9px;
}
.toggle:hover { color: var(--vp-c-text-1); }

.gas { fill: var(--mc-co2); }
.bead { fill: var(--mc-co2); stroke: var(--vp-c-text-2); stroke-width: 0.6; }
.bead-base { fill: var(--mc-bead); }
.bed-outline { fill: none; stroke: var(--vp-c-text-1); stroke-width: 1.6; }
.arrow { stroke: var(--vp-c-text-1); stroke-width: 1.6; fill: none; }
.arrowhead { fill: var(--vp-c-text-1); }
.stream { transition: opacity 0.3s; }
.label { fill: var(--vp-c-text-2); font: 12px var(--vp-font-family-base); }
.small { fill: var(--vp-c-text-2); font: 11px var(--vp-font-family-mono); }
.step-name { fill: var(--vp-c-text-1); font: 700 20px var(--mc-font-display); font-stretch: 116%; }
.readout { fill: var(--vp-c-text-2); font: 12px var(--vp-font-family-mono); }
.axis { stroke: var(--vp-c-text-2); stroke-width: 1; }
.grid { stroke: var(--vp-c-divider); stroke-width: 1; }
.pressure { fill: none; stroke: var(--vp-c-brand-1); stroke-width: 2; }
.cursor { stroke: var(--vp-c-text-2); stroke-width: 1; stroke-dasharray: 3 3; }
.dot { fill: var(--vp-c-brand-1); }

.st-pressurise { fill: var(--mc-step-pressurise); }
.st-adsorb { fill: var(--mc-step-adsorb); }
.st-blowdown { fill: var(--mc-step-blowdown); }
.st-evacuate { fill: var(--mc-step-evacuate); }
.st-equalise { fill: var(--mc-step-equalise); }
.st-lpp { fill: var(--mc-step-lpp); }
.bg-pressurise { background: var(--mc-step-pressurise); }
.bg-adsorb { background: var(--mc-step-adsorb); }
.bg-blowdown { background: var(--mc-step-blowdown); }
.bg-evacuate { background: var(--mc-step-evacuate); }
.bg-equalise { background: var(--mc-step-equalise); }
.bg-lpp { background: var(--mc-step-lpp); }
.bg-idle {
  background: repeating-linear-gradient(135deg, var(--vp-c-divider) 0 3px, transparent 3px 6px);
}

.legend { display: flex; flex-wrap: wrap; gap: 6px 14px; }
.chip { display: inline-flex; align-items: center; gap: 6px; font-size: 12.5px; color: var(--vp-c-text-2); }
.chip i { display: block; width: 14px; height: 8px; border-radius: 2px; }
.chip.on { color: var(--vp-c-text-1); font-weight: 600; }

/* ---------- sections ---------- */
.band { padding: 48px 0; border-bottom: 1px solid var(--vp-c-divider); display: grid; gap: 22px; }
.band.last { border-bottom: 0; padding-bottom: 72px; }
.band > header { display: grid; gap: 6px; max-width: 66ch; }
.band h2 {
  margin: 0;
  font-family: var(--mc-font-display);
  font-stretch: 116%;
  font-weight: 700;
  font-size: 24px;
  line-height: 1.2;
  text-wrap: balance;
}
.band header p { margin: 0; color: var(--vp-c-text-2); }
.band header a.inline { color: var(--vp-c-brand-1); font-weight: 500; }

.units { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 14px; }
.unit {
  display: grid;
  grid-template-rows: auto 1fr;
  gap: 12px;
  min-width: 0;
  padding: 16px;
  border: 1px solid var(--vp-c-divider);
  border-radius: 6px;
  color: var(--vp-c-text-2);
}
.unit-icon { width: 64px; height: 56px; }
.unit-text { display: grid; gap: 4px; align-content: start; }
.unit h3 { margin: 0; font: 600 15px/1.3 var(--vp-font-family-base); color: var(--vp-c-text-1); }
.unit p { margin: 0; font-size: 13.5px; line-height: 1.5; }
.status {
  justify-self: start;
  font: 500 11px/1 var(--vp-font-family-mono);
  letter-spacing: 0.04em;
  text-transform: uppercase;
  padding: 4px 6px;
  border-radius: 3px;
  border: 1px solid currentColor;
}
.unit.is-available { border-color: var(--vp-c-text-1); color: var(--vp-c-text-1); }
.unit.is-available .unit-icon { color: var(--mc-co2); }
.unit.is-available .status { background: var(--mc-step-adsorb); border-color: var(--mc-step-adsorb); color: #fff; }
.unit.is-available:hover { border-color: var(--vp-c-brand-1); }
.unit.is-development .status { color: var(--mc-step-blowdown); }
.unit.is-future { border-style: dashed; }
.unit.is-future .unit-icon { opacity: 0.55; }
.unit.is-future .status { border-style: dashed; }
.aside { margin: 0; font-size: 14px; color: var(--vp-c-text-2); }
.aside a { color: var(--vp-c-brand-1); font-weight: 500; }
.sheet tr.group th { width: auto; padding-top: 22px; font: 600 13px var(--vp-font-family-base); letter-spacing: 0; text-transform: none; color: var(--vp-c-text-1); }
.sheet tr.group:first-child th { padding-top: 4px; }

.pillars { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 16px; }
.pillar {
  display: grid;
  grid-template-rows: auto auto auto 1fr;
  gap: 8px;
  min-width: 0;
  padding: 18px;
  border: 1px solid var(--vp-c-divider);
  border-radius: 6px;
}
.pillar-n { font: 600 12px/1 var(--vp-font-family-mono); color: var(--mc-co2); }
.pillar h3 { margin: 0; font-family: var(--mc-font-display); font-stretch: 116%; font-weight: 700; font-size: 19px; }
.pillar p { margin: 0; font-size: 14px; color: var(--vp-c-text-2); }
.pillar pre { margin: 6px 0 0; padding: 12px 14px; background: var(--vp-c-bg-soft); border-radius: 4px; overflow-x: auto; align-self: start; }
.pillar pre code { font: 13px/1.65 var(--vp-font-family-mono); color: var(--vp-c-text-1); }
.pillar .k { color: var(--vp-c-brand-1); }
.pillar .s { color: var(--mc-step-adsorb); }
.nums { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 10px; margin-top: 6px; }
.num { display: grid; gap: 2px; min-width: 0; }
.num .val { display: flex; align-items: baseline; gap: 4px; flex-wrap: wrap; }
.num b { font-family: var(--mc-font-display); font-stretch: 112%; font-weight: 700; font-size: clamp(22px, 2.4vw, 28px); line-height: 1; font-variant-numeric: tabular-nums; }
.num small { font: 12px var(--vp-font-family-mono); color: var(--vp-c-text-2); }
.num span { font-size: 13px; color: var(--vp-c-text-2); }

.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 16px; }
.card {
  display: grid;
  gap: 12px;
  align-content: start;
  min-width: 0;
  padding: 18px;
  border: 1px solid var(--vp-c-divider);
  border-radius: 6px;
  transition: border-color 0.15s;
}
.card:hover { border-color: var(--vp-c-brand-1); }
.card h3 { margin: 0; font: 600 16px/1.3 var(--vp-font-family-base); }
.cite { font: 12px var(--vp-font-family-mono); color: var(--vp-c-text-2); }
.card p { margin: 0; font-size: 14px; color: var(--vp-c-text-2); }
.sched { display: grid; gap: 4px; }
.sched-row { display: grid; grid-template-columns: 30px minmax(0, 1fr); align-items: center; gap: 8px; font: 11px var(--vp-font-family-mono); color: var(--vp-c-text-2); }
.sched-bar { display: flex; gap: 1px; height: 14px; border-radius: 2px; overflow: hidden; }
.sched-bar i { display: block; height: 100%; }

.sheet-wrap { overflow-x: auto; }
.sheet { width: 100%; min-width: 520px; border-collapse: collapse; font-size: 15px; display: table; margin: 0; }
.sheet th, .sheet td { text-align: left; vertical-align: top; padding: 11px 14px 11px 0; border: 0; border-top: 1px solid var(--vp-c-divider); background: none; }
.sheet tr { background: none !important; border: 0; }
.sheet th { width: 120px; padding-top: 13px; font: 500 12px/1.5 var(--vp-font-family-mono); letter-spacing: 0.06em; text-transform: uppercase; color: var(--vp-c-text-2); }
.sheet a { color: var(--vp-c-brand-1); font-weight: 500; }
.sheet tr.not th, .sheet tr.not td { color: var(--vp-c-text-2); }

.examples { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 0 28px; }
.ex { display: grid; gap: 4px; align-content: start; padding: 14px 0; border-top: 1px solid var(--vp-c-divider); }
.ex .tag { font: 500 11px var(--vp-font-family-mono); letter-spacing: 0.06em; text-transform: uppercase; color: var(--vp-c-text-2); }
.ex h3 { margin: 0; font: 600 15px/1.35 var(--vp-font-family-base); }
.ex:hover h3 { color: var(--vp-c-brand-1); }
.ex p { margin: 0; font-size: 14px; color: var(--vp-c-text-2); }

.next { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 20px 28px; }
.next div { display: grid; gap: 6px; align-content: start; }
.next h3 { margin: 0 0 4px; font: 500 12px var(--vp-font-family-mono); letter-spacing: 0.06em; text-transform: uppercase; color: var(--vp-c-text-2); }
.next a { color: var(--vp-c-brand-1); font-weight: 500; }
.next a:hover { text-decoration: underline; }

.mh a:focus-visible, .toggle:focus-visible { outline: 2px solid var(--vp-c-brand-1); outline-offset: 2px; }

@media (max-width: 860px) {
  .hero { grid-template-columns: 1fr; padding-top: 36px; }
  .pillars { grid-template-columns: 1fr; }
}
@media (max-width: 560px) {
  .mh { padding: 0 16px; }
  .nums { grid-template-columns: repeat(3, minmax(0, 1fr)); }
  .band { padding: 36px 0; }
}
@media (prefers-reduced-motion: reduce) {
  .stream { transition: none; }
}
</style>
