```@meta
CurrentModule = Mocca
```

# Model equations

A fixed bed in Mocca is a one-dimensional, axially dispersed plug flow model
of a packed column, as in Haghpanah et al. (2013), or, with `radial_cells`,
one that also resolves the column's radius (see
[Beds meshed across their radius](@ref)). The gas is ideal, flow
follows Darcy's law, mass transfer to the sorbent is a linear driving force,
and the energy balance covers the gas, the sorbent and the column wall. The
beds are connected to the feed, products and each other through flow devices
(valves and flow controllers), whose settings change from step to step of the
cycle.

All quantities are in SI units. The symbols are listed in the
[nomenclature](#Nomenclature) at the end of this page, with the inputs they
come from.

```@raw html
<figure class="mocca-schematic">
<svg viewBox="0 0 680 420" role="img" aria-label="Schematic of a packed column: gas between particles, adsorbent particles, the column wall, heat transfer to the wall and the surroundings, and the ports at each end">
<defs>
<pattern id="cs-hatch" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="6" style="stroke:var(--vp-c-text-2);stroke-width:1.2;opacity:.55"/></pattern>
<marker id="cs-ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path class="cs-arrowhead" d="M0,0 L10,5 L0,10 z"/></marker>
<marker id="cs-hh" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path class="cs-heathead" d="M0,0 L10,5 L0,10 z"/></marker>
</defs>
<rect class="cs-wall" x="160" y="50" width="100" height="320" rx="6"/>
<rect class="cs-bed" x="172" y="50" width="76" height="320"/>
<circle class="cs-part" cx="180" cy="61" r="4.8"/>
<circle class="cs-part" cx="192" cy="61" r="4.8"/>
<circle class="cs-part" cx="204" cy="61" r="4.8"/>
<circle class="cs-part" cx="216" cy="61" r="4.8"/>
<circle class="cs-part" cx="228" cy="61" r="4.8"/>
<circle class="cs-part" cx="240" cy="61" r="4.8"/>
<circle class="cs-part" cx="186" cy="71.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="71.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="71.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="71.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="71.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="82.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="82.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="82.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="82.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="82.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="82.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="92.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="92.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="92.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="92.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="92.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="103.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="103.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="103.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="103.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="103.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="103.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="113.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="113.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="113.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="113.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="113.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="124.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="124.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="124.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="124.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="124.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="124.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="134.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="134.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="134.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="134.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="134.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="145.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="145.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="145.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="145.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="145.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="145.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="155.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="155.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="155.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="155.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="155.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="166.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="166.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="166.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="166.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="166.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="166.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="176.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="176.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="176.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="176.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="176.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="187.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="187.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="187.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="187.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="187.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="187.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="197.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="197.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="197.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="197.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="197.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="208.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="208.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="208.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="208.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="208.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="208.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="218.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="218.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="218.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="218.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="218.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="229.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="229.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="229.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="229.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="229.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="229.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="239.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="239.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="239.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="239.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="239.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="250.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="250.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="250.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="250.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="250.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="250.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="260.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="260.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="260.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="260.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="260.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="271.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="271.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="271.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="271.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="271.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="271.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="281.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="281.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="281.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="281.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="281.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="292.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="292.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="292.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="292.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="292.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="292.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="302.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="302.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="302.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="302.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="302.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="313.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="313.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="313.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="313.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="313.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="313.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="323.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="323.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="323.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="323.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="323.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="334.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="334.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="334.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="334.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="334.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="334.0" r="4.8"/>
<circle class="cs-part" cx="186" cy="344.5" r="4.8"/>
<circle class="cs-part" cx="198" cy="344.5" r="4.8"/>
<circle class="cs-part" cx="210" cy="344.5" r="4.8"/>
<circle class="cs-part" cx="222" cy="344.5" r="4.8"/>
<circle class="cs-part" cx="234" cy="344.5" r="4.8"/>
<circle class="cs-part" cx="180" cy="355.0" r="4.8"/>
<circle class="cs-part" cx="192" cy="355.0" r="4.8"/>
<circle class="cs-part" cx="204" cy="355.0" r="4.8"/>
<circle class="cs-part" cx="216" cy="355.0" r="4.8"/>
<circle class="cs-part" cx="228" cy="355.0" r="4.8"/>
<circle class="cs-part" cx="240" cy="355.0" r="4.8"/>
<line class="cs-line" x1="160" y1="50" x2="260" y2="50"/><line class="cs-line" x1="160" y1="370" x2="260" y2="370"/>
<line class="cs-line" x1="122" y1="370" x2="122" y2="58" marker-end="url(#cs-ah)"/>
<text class="cs-sym" x="112" y="374" text-anchor="end">0</text><text class="cs-sym" x="112" y="64" text-anchor="end">L</text><text class="cs-sym" x="130" y="58">z</text>
<line class="cs-line" x1="210" y1="412" x2="210" y2="378" marker-start="url(#cs-ah)" marker-end="url(#cs-ah)"/>
<text class="cs-text" x="222" y="404">bottom port, <tspan class="cs-sym">z</tspan> = 0</text>
<line class="cs-line" x1="210" y1="42" x2="210" y2="8" marker-start="url(#cs-ah)" marker-end="url(#cs-ah)"/>
<text class="cs-text" x="222" y="22">top port, <tspan class="cs-sym">z</tspan> = <tspan class="cs-sym">L</tspan></text>
<line class="cs-thin" x1="243" y1="96" x2="296" y2="96"/><circle cx="243" cy="96" r="2" class="cs-pore"/>
<text class="cs-head" x="302" y="92">Gas between the particles</text>
<text class="cs-text" x="302" y="110"><tspan class="cs-sym">P</tspan>, <tspan class="cs-sym">T</tspan>, <tspan class="cs-sym">y<tspan class="cs-sub" dy="3">i</tspan><tspan dy="-3">&#8203;</tspan></tspan>, velocity <tspan class="cs-sym">u</tspan> · bed porosity <tspan class="cs-sym">ε</tspan></text>
<line class="cs-thin" x1="240" y1="176" x2="402" y2="176"/>
<circle cx="460" cy="176" r="58" style="fill:var(--vp-c-bg);stroke:var(--vp-c-text-1);stroke-width:1.2"/>
<circle class="cs-part" cx="440" cy="152" r="20"/>
<circle class="cs-pore" cx="440.7" cy="161.0" r="1.6"/>
<circle class="cs-pore" cx="433.4" cy="159.0" r="1.6"/>
<circle class="cs-pore" cx="437.4" cy="149.4" r="1.6"/>
<circle class="cs-pore" cx="452.2" cy="153.0" r="1.6"/>
<circle class="cs-pore" cx="439.7" cy="157.6" r="1.6"/>
<circle class="cs-pore" cx="448.2" cy="151.8" r="1.6"/>
<circle class="cs-pore" cx="444.3" cy="144.9" r="1.6"/>
<circle class="cs-pore" cx="437.0" cy="148.4" r="1.6"/>
<circle class="cs-pore" cx="431.7" cy="142.6" r="1.6"/>
<circle class="cs-part" cx="484" cy="160" r="20"/>
<circle class="cs-pore" cx="473.0" cy="158.4" r="1.6"/>
<circle class="cs-pore" cx="482.2" cy="156.7" r="1.6"/>
<circle class="cs-pore" cx="484.5" cy="150.5" r="1.6"/>
<circle class="cs-pore" cx="482.9" cy="163.2" r="1.6"/>
<circle class="cs-pore" cx="489.4" cy="153.9" r="1.6"/>
<circle class="cs-pore" cx="481.5" cy="147.6" r="1.6"/>
<circle class="cs-pore" cx="481.1" cy="147.2" r="1.6"/>
<circle class="cs-pore" cx="474.7" cy="167.2" r="1.6"/>
<circle class="cs-pore" cx="471.5" cy="164.5" r="1.6"/>
<circle class="cs-part" cx="446" cy="200" r="20"/>
<circle class="cs-pore" cx="448.9" cy="197.2" r="1.6"/>
<circle class="cs-pore" cx="449.5" cy="204.1" r="1.6"/>
<circle class="cs-pore" cx="453.6" cy="198.3" r="1.6"/>
<circle class="cs-pore" cx="441.6" cy="195.5" r="1.6"/>
<circle class="cs-pore" cx="438.8" cy="199.7" r="1.6"/>
<circle class="cs-pore" cx="440.4" cy="207.6" r="1.6"/>
<circle class="cs-pore" cx="434.8" cy="193.5" r="1.6"/>
<circle class="cs-pore" cx="440.5" cy="188.0" r="1.6"/>
<circle class="cs-pore" cx="454.6" cy="189.1" r="1.6"/>
<circle class="cs-part" cx="490" cy="206" r="20"/>
<circle class="cs-pore" cx="487.7" cy="201.8" r="1.6"/>
<circle class="cs-pore" cx="498.7" cy="195.5" r="1.6"/>
<circle class="cs-pore" cx="497.6" cy="200.8" r="1.6"/>
<circle class="cs-pore" cx="488.8" cy="200.8" r="1.6"/>
<circle class="cs-pore" cx="494.6" cy="197.9" r="1.6"/>
<circle class="cs-pore" cx="489.2" cy="209.6" r="1.6"/>
<circle class="cs-pore" cx="498.4" cy="195.0" r="1.6"/>
<circle class="cs-pore" cx="500.0" cy="212.2" r="1.6"/>
<circle class="cs-pore" cx="486.1" cy="208.5" r="1.6"/>
<text class="cs-sym" x="462" y="186">ε</text>
<text class="cs-head" x="530" y="140">Adsorbent particle</text>
<text class="cs-text" x="530" y="158">diameter <tspan class="cs-sym">d<tspan class="cs-sub" dy="3">p</tspan><tspan dy="-3">&#8203;</tspan></tspan></text>
<text class="cs-text" x="530" y="175">porosity <tspan class="cs-sym">ε<tspan class="cs-sub" dy="3">p</tspan><tspan dy="-3">&#8203;</tspan></tspan></text>
<text class="cs-text" x="530" y="192">tortuosity <tspan class="cs-sym">τ</tspan></text>
<text class="cs-text" x="530" y="209">loading <tspan class="cs-sym">q<tspan class="cs-sub" dy="3">i</tspan><tspan dy="-3">&#8203;</tspan></tspan> → <tspan class="cs-sym">q<tspan class="cs-sub" dy="3">i</tspan><tspan dy="-3">&#8203;</tspan></tspan>*</text>
<line class="cs-thin" x1="256" y1="252" x2="296" y2="252"/>
<text class="cs-head" x="302" y="248">Column wall</text>
<text class="cs-text" x="302" y="266">temperature <tspan class="cs-sym">T<tspan class="cs-sub" dy="3">w</tspan><tspan dy="-3">&#8203;</tspan></tspan> · <tspan class="cs-sym">ρ<tspan class="cs-sub" dy="3">w</tspan><tspan dy="-3">&#8203;</tspan></tspan>, <tspan class="cs-sym">c<tspan class="cs-sub" dy="3">pw</tspan><tspan dy="-3">&#8203;</tspan></tspan>, <tspan class="cs-sym">K<tspan class="cs-sub" dy="3">w</tspan><tspan dy="-3">&#8203;</tspan></tspan></text>
<line class="cs-heat" x1="228" y1="306" x2="257" y2="306" marker-end="url(#cs-hh)"/>
<line class="cs-heat" x1="262" y1="306" x2="292" y2="306" marker-end="url(#cs-hh)"/>
<text class="cs-text" x="302" y="302"><tspan class="cs-sym">h<tspan class="cs-sub" dy="3">in</tspan><tspan dy="-3">&#8203;</tspan></tspan>: bed to wall</text>
<text class="cs-text" x="302" y="320"><tspan class="cs-sym">h<tspan class="cs-sub" dy="3">out</tspan><tspan dy="-3">&#8203;</tspan></tspan>: wall to surroundings at <tspan class="cs-sym">T<tspan class="cs-sub" dy="3">a</tspan><tspan dy="-3">&#8203;</tspan></tspan></text>
<circle cx="600" cy="330" r="38" style="fill:url(#cs-hatch);stroke:var(--vp-c-text-1);stroke-width:1.2"/>
<circle class="cs-bed" cx="600" cy="330" r="31" style="stroke:var(--vp-c-text-1);stroke-width:1.2"/>
<line class="cs-line" x1="600" y1="330" x2="631" y2="330" marker-end="url(#cs-ah)"/>
<line class="cs-line" x1="600" y1="330" x2="573" y2="357" marker-end="url(#cs-ah)"/>
<text class="cs-sym" x="604" y="324">r<tspan class="cs-sub" dy="3">in</tspan><tspan dy="-3">&#8203;</tspan></text>
<text class="cs-sym" x="556" y="370" text-anchor="end">r<tspan class="cs-sub" dy="3">out</tspan><tspan dy="-3">&#8203;</tspan></text>
<text class="cs-text" x="600" y="282" text-anchor="middle">cross-section</text>
</svg>
<figcaption>A fixed bed. Gas flows through the voids between the particles in either direction, entering or leaving through the ports at each end. The particles take up CO₂ towards the equilibrium loading of the isotherm, and heat passes from the bed to the wall and from the wall to the surroundings. Symbols as in the <a href="#Nomenclature">nomenclature</a>.</figcaption>
</figure>

```

## Fixed bed

### Component mass balance

For each component ``i`` of the gas,

```math
\varepsilon \frac{\partial c_i}{\partial t}
+ \frac{\partial (u\, c_i)}{\partial z}
- \frac{\partial}{\partial z}\left(\varepsilon D_L C \frac{\partial y_i}{\partial z}\right)
+ (1-\varepsilon)\frac{\partial q_i}{\partial t} = 0,
\qquad c_i = y_i C = \frac{y_i P}{RT}.
```

``q_i`` is the loading in mol per m³ of adsorbent (particles), so
``(1-\varepsilon)`` converts it to bed volume.

### Momentum: Darcy's law

The superficial velocity follows Darcy's law with the Kozeny–Carman
(Blake–Kozeny) permeability of a packed bed,

```math
u = -\frac{k}{\mu}\frac{\partial P}{\partial z},
\qquad
k = \frac{\varepsilon^3 d_p^2}{150\,(1-\varepsilon)^2}.
```

This is the viscous (laminar) term of the Ergun equation; the inertial term is
not included. The interstitial velocity is ``v = u/\varepsilon``. Give
`permeability` to [`add_bed!`](@ref) to override ``k``.

### Axial dispersion

```math
D_L = 0.7\, D_m + 0.5\, v_0\, d_p,
```

where ``D_m`` is the molecular diffusivity of the gas (`diffusivity` of
[`Gas`](@ref), that of the first component if given per component) and ``v_0`` the interstitial feed velocity of the bed's
adsorption step. ``D_L`` is evaluated once, when the process is built, and is
constant along the bed and through the cycle. Set `dispersion_velocity` or
`dispersion` in [`add_bed!`](@ref) to fix it yourself.

### Adsorption rate: linear driving force

```math
\frac{\partial q_i}{\partial t} = k_i \left(q_i^* - q_i\right),
\qquad
k_i = \frac{c_i}{q_i^*}\,\frac{15\,\varepsilon_p D_{m,i}}{\tau\, r_p^2},
\qquad r_p = \frac{d_p}{2}.
```

This is macropore diffusion control: the particle porosity ``\varepsilon_p``,
tortuosity ``\tau`` and diameter ``d_p`` of the [`Adsorbent`](@ref) set the rate.

### Isotherm: dual-site Langmuir

```math
q_i^* = \frac{q_{sb,i}\, b_i c_i}{1 + \sum_j b_j c_j}
      + \frac{q_{sd,i}\, d_i c_i}{1 + \sum_j d_j c_j},
\qquad
b_i = b_{0,i}\, e^{-\Delta U_{b,i}/RT},
\quad
d_i = d_{0,i}\, e^{-\Delta U_{d,i}/RT}.
```

The affinities are on a concentration basis (m³/mol), and the saturation
capacities in mol per m³ of adsorbent, or in mol/kg with `basis = :mass`,
which Mocca multiplies by the adsorbent density. The single-site
[`langmuir`](@ref) isotherm has ``q_{sd,i} = 0``. Check an isotherm against its
source with [`equilibrium_loading`](@ref) and [`isotherm_table`](@ref), which
report mol/kg.

The pure-component isotherms of zeolite 13X from Haghpanah et al. (2013), at
three temperatures. CO₂ is held far more strongly than N₂, and both less at
higher temperature, so the heat released by adsorption, which warms the bed,
lowers what it can hold:

```@example model
using Mocca, CairoMakie # hide
gas = Mocca.haghpanah_2013_gas() # hide
ads = Mocca.load_adsorbent("zeolite_13x_haghpanah_2013") # hide
p = range(0, 1e5, length = 101) # hide
T = [298.15, 323.15, 348.15] # hide
fig = Figure(size = (760, 300)) # hide
axes = [Axis(fig[1, j], title = name, xlabel = "partial pressure [bar]", ylabel = j == 1 ? "loading [mol/kg]" : "") for (j, name) in enumerate(["CO2", "N2"])] # hide
for (ax, name) in zip(axes, ["CO2", "N2"]) # hide
    q = Mocca.isotherm_table(ads, gas, name; pressures = p, temperatures = T) # hide
    for k in eachindex(T) # hide
        lines!(ax, p ./ 1e5, q[:, k], label = "$(round(Int, T[k] - 273.15)) °C") # hide
    end # hide
end # hide
axislegend(axes[1], position = :rb) # hide
fig # hide
```

The heat of adsorption of each component is constant, the average of its two
site energies weighted by their capacities at the reference temperature
``T_0`` (`reference_temperature`, default 298.15 K):

```math
\Delta H_i = \frac{q_{sb,i}\left(\Delta U_{b,i} - RT_0\right) + q_{sd,i}\left(\Delta U_{d,i} - RT_0\right)}{q_{sb,i} + q_{sd,i}}.
```

### Energy balance of the column

Written in terms of temperature and pressure:

```math
\begin{aligned}
&\left[(1-\varepsilon)\left(\rho_s c_{ps} + \tilde c_{pa} \sum_i q_i\right)\right]\frac{\partial T}{\partial t}
+ \varepsilon\,\frac{\tilde c_{pg}}{R}\frac{\partial P}{\partial t}
+ \frac{\tilde c_{pg}}{R}\frac{\partial (uP)}{\partial z}
- \frac{\partial}{\partial z}\left(K_z \frac{\partial T}{\partial z}\right) \\
&\qquad + (1-\varepsilon)\sum_i \left(\tilde c_{pa} T + \Delta H_i\right)\frac{\partial q_i}{\partial t}
+ \frac{2 h_{in}}{r_{in}}\left(T - T_w\right) = 0,
\end{aligned}
```

with molar heat capacities of the gas and the adsorbed phase from the
per-mass values of each component,

```math
\tilde c_{pg} = \bar M \sum_i y_i\, c_{pg,i},
\qquad
\tilde c_{pa} = \bar M \sum_i y_i\, c_{pa,i},
\qquad
\bar M = \sum_i y_i M_i.
```

``\Delta H_i < 0``, so adsorption heats the bed.

### Wall energy balance

```math
\rho_w c_{pw}\frac{\partial T_w}{\partial t}
= K_w \frac{\partial^2 T_w}{\partial z^2}
+ \frac{2 r_{in} h_{in}}{r_{out}^2 - r_{in}^2}\left(T - T_w\right)
- \frac{2 r_{out} h_{out}}{r_{out}^2 - r_{in}^2}\left(T_w - T_a\right).
```

At both ends the wall conducts heat to the surroundings at ``T_a``.

### Thermal options

The `thermal` keyword of [`add_bed!`](@ref) chooses the energy model:

| `thermal` | Energy model |
|---|---|
| `WithWall()` (default) | Column and wall energy balances, as above |
| `Adiabatic()` | Column energy balance without the wall term (``h_{in} = 0``) |
| `Isothermal(T)` | Temperature fixed at `T`, no energy balance |

## Beds meshed across their radius

With `radial_cells` in [`add_bed!`](@ref), the bed resolves the gradients
across its radius as well as along it: the variables depend on ``z`` and the
radius ``r``, and the bed is symmetric about its axis. This matters where
heat leaves or enters through the wall faster than it crosses the bed, which
the Biot number of the cross-section measures,

```math
\mathrm{Bi} = \frac{h_w\, r_{in}}{K_r}.
```

Well below 1 the temperature is nearly flat across the radius and the 1D bed
is enough; around 1 and above, the centre and the wall differ. Large columns,
beds of low conductivity and strong heat exchange at the wall raise it.

### Balances

The component mass balance gains dispersion and flow across the radius,

```math
\varepsilon \frac{\partial c_i}{\partial t}
+ \frac{\partial (u_z c_i)}{\partial z}
+ \frac{1}{r}\frac{\partial (r\, u_r c_i)}{\partial r}
- \frac{\partial}{\partial z}\left(\varepsilon D_L C \frac{\partial y_i}{\partial z}\right)
- \frac{1}{r}\frac{\partial}{\partial r}\left(r\, \varepsilon D_R C \frac{\partial y_i}{\partial r}\right)
+ (1-\varepsilon)\frac{\partial q_i}{\partial t} = 0,
```

with Darcy's law in both directions, ``u_z = -(k/\mu)\,\partial P/\partial z``
and ``u_r = -(k/\mu)\,\partial P/\partial r``. The radial dispersion is

```math
D_R = 0.7\, D_m + 0.1\, v_0\, d_p,
```

the form of the axial dispersion with a radial Péclet number ``v_0 d_p / D_R``
of 10 at high velocity instead of 2 (see
[`compute_radial_dispersion`](@ref)); give `radial_dispersion` to
[`add_bed!`](@ref) to set it yourself. The energy balance of the column gains
conduction across the radius, ``\frac{1}{r}\frac{\partial}{\partial r}\left(r
K_r \frac{\partial T}{\partial r}\right)``, and the flow of enthalpy across it,
and loses the wall term, which becomes a boundary condition.

### The wall

The wall surrounds the bed at ``r = r_{in}``, where the heat conducted to it
crosses a film at the wall,

```math
-K_r \frac{\partial T}{\partial r}\bigg|_{r_{in}} = h_w \left(T(r_{in}) - T_w\right),
```

and the wall energy balance is that of the 1D bed with ``h_w`` and the gas
temperature at the wall in place of ``h_{in}`` and ``T``. At the axis the
gradients vanish. ``h_w`` (`wall_heat_transfer` of [`Column`](@ref)) is the
coefficient of the film alone: the ``h_{in}`` of a 1D bed (`inner_heat_transfer`)
also lumps the resistance of the bed across its radius. For heat released
evenly in the bed, the 1D bed whose mean temperature is that of the bed
across its radius has

```math
\frac{1}{h_{in}} = \frac{1}{h_w} + \frac{r_{in}}{4 K_r},
```

which [`cyclic_steady_state`](@ref) uses for the 1D bed it starts a radial
bed's steady state from.

### Discretisation

Each of the `ncells` cells along the column is split into `radial_cells`
rings, uniform in radius unless `radii` gives their edges. A ring from
``r_1`` to ``r_2`` holds ``\pi (r_2^2 - r_1^2)\,\Delta z``; the face between
two rings at radius ``r`` has the area ``2\pi r\, \Delta z``, and the face
between two cells along the column the ring's cross-section. The fluxes are
two-point, as along the column, over the distances between the centres of
the rings.

- The ports are the rings at each end. A stream into the bed is shared
  between them as their cross-sections are; a stream out of it takes each
  ring's own gas, and the device sees the rings mixed, at the mixing-cup
  temperature, so that every component and the enthalpy are conserved.
- The wall exchanges heat with the outer ring only, through the film and the
  half of the ring between its centre and the wall,
  ``1/h = 1/h_w + r_{in} \ln(r_{in}/r_c)/K_r`` for the ring's centre ``r_c``
  ([`RingWallResistance`](@ref)). Without the half ring the error would be
  first order in the width of the rings; with it, the cooling of a bed by
  conduction alone matches the series solution to second order (see
  [Conduction across the radius](@ref)).
- The cells inside the outer ring have no wall: their wall temperature
  follows their gas temperature, and has no meaning.

A bed of one radial cell is the 1D bed of the sections above.

## Flow devices and boundary conditions

Each end of a bed (`:bottom` at ``z = 0``, `:top` at ``z = L``) connects to one
or more flow devices. A device has an inlet and an outlet port and carries a
molar flow ``F`` [mol/s], positive from inlet to outlet. The port attached to
the bed takes the pressure, temperature and composition of the end cell. The
other port is either another bed (a link) or a boundary at a prescribed
pressure, temperature and composition. The flow enters or leaves the bed's
end cell carrying the composition and enthalpy ``F\,\tilde c_{pg} T`` of the
upstream port.

In each step of the cycle every device is closed, a valve or a flow
controller:

| Law | Flow |
|---|---|
| Closed | ``F = 0`` |
| [`LinearValve`](@ref) | ``F = C_v\,(P_{in} - P_{out})\,\dfrac{P_{up}}{R\,T_{up}}`` |
| [`VolumetricFlow`](@ref) (flow controller) | ``F = Q\,\dfrac{P_{out}}{R\,T_{in}}`` |

Here ``C_v`` is the valve conductance [m³/(s·Pa)] and "up" is the upstream port.
A feed given as an interstitial velocity ``v`` is a flow controller with
``Q = \varepsilon v A``.

By default a valve has the conductance of half a cell at the bed end it is
attached to, ``C_v = kA/(\mu\,\Delta z/2)``, so that it adds no more resistance
than the bed itself and the bed end follows the pressure on the far side of
the valve. Give `conductance` to a step to model a smaller valve.

During pressurisation, blowdown and evacuation, the pressure on the far side
of the valve moves exponentially from where the step starts to where it ends:

```math
P(t) = P_{end} + \left(P_{start} - P_{end}\right) e^{-\lambda t},
```

with ``\lambda`` the `ramp_rate` (default 0.5 s⁻¹) of the step or the
[`Process`](@ref). A blowdown from 1 to 0.1 bar at three ramp rates; at the
default rate the pressure on the far side of the valve has made 99% of its
change after about 10 s:

```@example model
using CairoMakie # hide
t = range(0, 30, length = 301) # hide
fig = Figure(size = (560, 300)) # hide
ax = Axis(fig[1, 1], xlabel = "time in the step [s]", ylabel = "pressure [bar]") # hide
for λ in (0.1, 0.5, 2.0) # hide
    lines!(ax, t, 0.1 .+ (1.0 - 0.1) .* exp.(-λ .* t), label = "λ = $λ s⁻¹") # hide
end # hide
axislegend(ax) # hide
fig # hide
```

For the four-step VSA cycle this gives the boundary conditions of Haghpanah et
al. (2013):

| Step | Bottom (``z = 0``) | Top (``z = L``) |
|---|---|---|
| Pressurisation | feed valve, ``P`` ramps ``P_L \to P_H`` | closed |
| Adsorption | feed at velocity ``v_{feed}`` | product valve at ``P_H`` |
| Blowdown | closed | product valve, ``P`` ramps ``P_H \to P_I`` |
| Evacuation | vacuum valve, ``P`` ramps ``P_I \to P_L`` | closed |

## Performance

[`performance`](@ref) measures a process over its last simulated cycle, from
the streams through its devices. The device roles (see [`add_device!`](@ref))
say which streams count. For the key component ``k`` (the first component of
the gas unless `key_component` is given to [`Process`](@ref)), with ``N_i^{r}``
the moles of component ``i`` that passed through all devices of role ``r``
over the cycle:

```math
\text{purity} = \frac{N_k^{extract}}{\sum_i N_i^{extract}},
\qquad
\text{recovery} = \frac{N_k^{extract}}{N_k^{feed}},
\qquad
\text{productivity} = \frac{N_k^{extract}}{V_{ads}\; t_{cycle}},
```

with ``V_{ads}`` the volume of adsorbent, ``(1-\varepsilon)`` times the bed
volume, summed over the beds. Productivity is in mol/(m³ adsorbent · s);
multiply by ``M_k \cdot 86400/1000`` for tonnes per m³ per day.

The energy counts a vacuum pump on every outlet device (roles `:product`,
`:extract`, `:waste`) while it discharges below ambient pressure, and a
compressor on every inlet device (`:feed`, `:light_inlet`) while it delivers
above ambient pressure. Both compress adiabatically between the ambient pressure
``P_{amb}`` and the device's pressure ``P``:

```math
E = \int \frac{F R T}{\eta}\,\frac{\gamma}{\gamma - 1}
\left[\left(\frac{P_{high}}{P_{low}}\right)^{(\gamma-1)/\gamma} - 1\right] dt,
\qquad
\text{specific energy} = \frac{E}{N_k^{extract} M_k},
```

reported in kWh per tonne of key component. ``P_{amb}``, ``\eta`` and
``\gamma`` are the `ambient_pressure` (101325 Pa), `pump_efficiency` (0.72) and
`heat_capacity_ratio` (1.4) of the [`Process`](@ref). Haghpanah et al. count
the vacuum pumps only; [`vacuum_pump_energy`](@ref) gives the energy of one
device, to report it their way.

## Steady state of the cycle

The cycle maps the state of the beds at its start, ``x``, to the state at its
end, ``\Phi(x)``. Performance is measured at the cyclic steady state,
``x = \Phi(x)``. How Mocca finds it, and its options, are on the
[Cyclic steady state](@ref) page.

## Numerical method

Each bed is divided into `ncells` equal cells (30 by default), and each of
those into `radial_cells` rings for a bed meshed across its radius. Fluxes
between cells are two-point, with upwinded advection, and time stepping is fully
implicit (backward Euler), solved by Newton's method with exact derivatives
from automatic differentiation. Each step of the cycle starts with time steps
of `first_dt` that double up to `max_dt` (see [`TimeStepping`](@ref)), which
carries the solver through the abrupt start of a step. Gradients for
optimisation come from the adjoint of the same discrete equations.

## Nomenclature

| Symbol | Meaning | Unit | Set by |
|---|---|---|---|
| ``A`` | Column cross-section, ``\pi r_{in}^2`` | m² | [`Column`](@ref) `radius` |
| ``b_0, d_0`` | Affinity at infinite temperature, sites b and d | m³/mol | isotherm `b0`, `d0` |
| ``c_i,\ C`` | Gas concentration of component ``i``, total | mol/m³ | state |
| ``c_{pa,i}`` | Heat capacity of adsorbed component ``i`` | J/(kg·K) | [`Adsorbent`](@ref) `adsorbed_heat_capacity` |
| ``c_{pg,i}`` | Heat capacity of gas component ``i`` | J/(kg·K) | [`Gas`](@ref) `heat_capacities` |
| ``c_{ps}`` | Heat capacity of the adsorbent | J/(kg·K) | [`Adsorbent`](@ref) `heat_capacity` |
| ``c_{pw}`` | Heat capacity of the wall | J/(kg·K) | [`Wall`](@ref) `heat_capacity` |
| ``C_v`` | Valve conductance | m³/(s·Pa) | step `conductance` |
| ``d_p`` | Particle diameter | m | [`Adsorbent`](@ref) `particle_diameter` |
| ``D_L`` | Axial dispersion coefficient | m²/s | computed, or [`add_bed!`](@ref) `dispersion` |
| ``D_R`` | Radial dispersion coefficient | m²/s | computed, or [`add_bed!`](@ref) `radial_dispersion` |
| ``D_m`` | Molecular diffusivity | m²/s | [`Gas`](@ref) `diffusivity` |
| ``F`` | Molar flow through a device | mol/s | state |
| ``h_{in}, h_{out}`` | Heat transfer coefficient, bed to wall and wall to ambient | W/(m²·K) | [`Column`](@ref) `inner_heat_transfer`, `outer_heat_transfer` |
| ``h_w`` | Heat transfer coefficient of the film at the wall, for a bed meshed across its radius | W/(m²·K) | [`Column`](@ref) `wall_heat_transfer` |
| ``k`` | Bed permeability | m² | computed, or [`add_bed!`](@ref) `permeability` |
| ``K_z`` | Axial thermal conductivity of the bed | W/(m·K) | [`Column`](@ref) `axial_conductivity` |
| ``K_r`` | Radial thermal conductivity of the bed | W/(m·K) | [`Column`](@ref) `radial_conductivity` |
| ``K_w`` | Thermal conductivity of the wall | W/(m·K) | [`Wall`](@ref) `conductivity` |
| ``L`` | Column length | m | [`Column`](@ref) `length` |
| ``M_i`` | Molar mass | kg/mol | [`Gas`](@ref) `molar_masses` |
| ``P`` | Pressure | Pa | state |
| ``q_i,\ q_i^*`` | Loading and equilibrium loading | mol/m³ adsorbent | state |
| ``q_{sb}, q_{sd}`` | Saturation capacity, sites b and d | mol/m³ adsorbent | isotherm `qsb`, `qsd` |
| ``r_{in}, r_{out}`` | Inner and outer radius of the column wall | m | [`Column`](@ref) `radius`, `outer_radius` |
| ``r`` | Radius, from the axis of the column | m | |
| ``R`` | Gas constant, 8.314 | J/(mol·K) | |
| ``T,\ T_w,\ T_a`` | Gas and sorbent, wall and ambient temperature | K | state; [`Column`](@ref) `ambient_temperature` |
| ``u,\ v`` | Superficial and interstitial velocity | m/s | state; step `velocity` |
| ``u_z,\ u_r`` | Superficial velocity along the column and across its radius | m/s | state |
| ``y_i`` | Gas mole fraction | – | state |
| ``\Delta U_b, \Delta U_d`` | Internal energy of adsorption, sites b and d | J/mol | isotherm `ΔUb`, `ΔUd` |
| ``\varepsilon`` | Bed porosity (void fraction between particles) | – | [`Column`](@ref) `porosity` |
| ``\varepsilon_p`` | Particle porosity | – | [`Adsorbent`](@ref) `particle_porosity` |
| ``\gamma,\ \eta`` | Heat capacity ratio and efficiency of pumps and compressors | – | [`Process`](@ref) |
| ``\lambda`` | Pressure ramp rate | 1/s | step or [`Process`](@ref) `ramp_rate` |
| ``\mu`` | Gas viscosity | Pa·s | [`Gas`](@ref) `viscosity` |
| ``\rho_s,\ \rho_w`` | Density of the adsorbent and of the wall | kg/m³ | [`Adsorbent`](@ref) `density`, [`Wall`](@ref) `density` |
| ``\tau`` | Tortuosity | – | [`Adsorbent`](@ref) `tortuosity` |
