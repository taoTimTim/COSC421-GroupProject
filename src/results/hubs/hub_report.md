# Hub Electrode Analysis - Research Question 1

Generated from R analysis outputs in `src/results/hubs`.

Date: 2025-11-26

---

**Summary**

This report summarizes the hub electrode analysis comparing meditation and thinking states for Alpha and Beta frequency bands. The analysis used weighted Phase Lag Index (wPLI) connectivity matrices thresholded at multiple densities (10%, 15%, 20%, 25%), with hub metrics averaged across these levels to reduce single-threshold bias. Centrality measures computed per electrode:

- Degree (unweighted)
- Strength (weighted)
- Betweenness (normalized)

Outputs used to create this report:

- `all_hub_metrics_multi_density.csv` — complete electrode metrics for all conditions (averaged across 10/15/20/25% density thresholds)
- `region_summary.csv` — aggregated region-level statistics (avg_strength, avg_degree, avg_betweenness, hub_count)
- `top5_*.csv` — top 5 hub electrodes per condition (six files)
- PDF visualizations in `plots/`:
  - `top5_hubs_per_condition.pdf`
  - `hub_region_distribution.pdf`
  - `meditation_vs_thinking_strength.pdf`

---

## Methods (brief)

- Input: 64A-64 wPLI average matrices per condition (alpha_med1, alpha_med2, alpha_thinking, beta_med1, beta_med2, beta_thinking).
- Thresholding: Proportional density thresholding applied at 10%, 15%, 20%, and 25%; metrics were averaged across these densities for stability.
- Graph construction: `igraph::graph_from_adjacency_matrix` with weighted edges.
- Centrality: `degree`, `strength` (weights), and `betweenness` (normalized).
- Region mapping: electrodes mapped to `frontal`, `central`, `parietal`, `occipital`, `temporal`, or `other`.

---

## Key Findings (high level)

- Alpha band: Occipital electrodes (`Oz`, `Iz`, `O1`) consistently rank highest in strength across meditation and thinking, indicating posterior hubs in alpha.
- Beta band: Hub locations are more distributed, with central electrodes (e.g., `CP1`) and frontal electrodes (`AF3`, `AFz`) appearing among top hubs in meditation.
- Region-level aggregates show higher avg_strength in occipital regions for alpha, while beta med/think strengths are lower and more distributed across regions.

---

## Region Summary

The table below is taken from `region_summary.csv` and shows aggregated region metrics.

| band | state | region | avg_strength | avg_degree | avg_betweenness | hub_count |
|---:|---|---|---:|---:|---:|---:|
"alpha" | "meditation" | "occipital" | 2.34959905892955 | 31.8125 | 0.0189452124935996 | 8
"alpha" | "meditation" | "other" | 0.707245458717064 | 10.5685483870968 | 0.00841963563087392 | 62
"alpha" | "meditation" | "frontal" | 0.668910377806363 | 9.88636363636364 | 0.00172811059907834 | 22
"alpha" | "meditation" | "central" | 0.540760139052733 | 8.16666666666667 | 0.00619417420492689 | 18
"alpha" | "meditation" | "parietal" | 0.536358192248964 | 8.30357142857143 | 0.0113287250384025 | 14
"alpha" | "meditation" | "temporal" | 0.296399732021051 | 5.125 | 0.00928059395801331 | 4
"alpha" | "thinking" | "occipital" | 2.24283430175711 | 31.0625 | 0.0259856630824373 | 4
"alpha" | "thinking" | "other" | 0.710800046351585 | 10.8629032258065 | 0.00806038683249921 | 31
"alpha" | "thinking" | "frontal" | 0.618521744752147 | 9.34090909090909 | 0.00145463855141274 | 11
"alpha" | "thinking" | "parietal" | 0.546375581281976 | 8.60714285714286 | 0.00972862263184844 | 7
"alpha" | "thinking" | "central" | 0.529384722965331 | 8.16666666666667 | 0.00702622745633498 | 9
"alpha" | "thinking" | "temporal" | 0.224500285029086 | 4 | 0.00902457757296467 | 2
"beta" | "meditation" | "occipital" | 0.241620030352871 | 20.4375 | 0.0284338197644649 | 8
"beta" | "meditation" | "central" | 0.152533594109661 | 12.1111111111111 | 0.0129288274449565 | 18
"beta" | "meditation" | "other" | 0.134730857262186 | 11.1653225806452 | 0.0115868060717176 | 62
"beta" | "meditation" | "frontal" | 0.126903905732256 | 10.3295454545455 | 0.00560908625424754 | 22
"beta" | "meditation" | "parietal" | 0.0810389185438882 | 7.28571428571429 | 0.0159004461999854 | 14
"beta" | "meditation" | "temporal" | 0.0228956412731583 | 2 | 0.00444828469022017 | 4
"beta" | "thinking" | "occipital" | 0.170753918658302 | 17.9375 | 0.025857654889913 | 4
"beta" | "thinking" | "central" | 0.144303928315642 | 13.8611111111111 | 0.0102406554019457 | 9
"beta" | "thinking" | "other" | 0.109954191449021 | 11.25 | 0.0147044249541648 | 31
"beta" | "thinking" | "parietal" | 0.0920427027504426 | 9.85714285714286 | 0.0213042206129764 | 7
"beta" | "thinking" | "frontal" | 0.0749792556359186 | 7.90909090909091 | 0.00161755806917097 | 11
"beta" | "thinking" | "temporal" | 0.0185418023551061 | 2.125 | 0.00083205325140809 | 2

---

## Top 5 hubs (per condition)

### Alpha - Meditation 1 (top5)

| electrode | degree | strength | betweenness | region |
|---|---:|---:|---:|---|
| Oz | 40.75 | 3.32767027459102 | 0.0135688684075781 | occipital |
| Iz | 40.25 | 3.12295023008275 | 0.0165130568356375 | occipital |
| P9 | 39.25 | 2.97162597020484 | 0.0295698924731183 | other |
| O1 | 40.00 | 2.9641440950287 | 0.0232974910394265 | occipital |
| PO7 | 37.75 | 2.61534747775696 | 0.0702764976958525 | parietal |

### Alpha - Meditation 2 (top5)

| electrode | degree | strength | betweenness | region |
|---|---:|---:|---:|---|
| Oz | 39.75 | 3.01574684127048 | 0.0345622119815668 | occipital |
| Iz | 38.75 | 2.84042827194754 | 0.0193292370711726 | occipital |
| P9 | 37.75 | 2.64068509431368 | 0.0421146953405018 | other |
| O1 | 35.00 | 2.33206284129532 | 0.0281618023553507 | occipital |
| PO7 | 33.00 | 2.09414565796147 | 0.0544034818228367 | parietal |

### Alpha - Thinking (top5)

| electrode | degree | strength | betweenness | region |
|---|---:|---:|---:|---|
| Oz | 39.75 | 3.0790107781804 | 0.0318740399385561 | occipital |
| Iz | 39.75 | 2.91009682757789 | 0.0211213517665131 | occipital |
| P9 | 36.75 | 2.56746860822984 | 0.0405785970302099 | other |
| O1 | 37.00 | 2.5236198034563 | 0.0450588837685612 | occipital |
| PO7 | 33.50 | 2.16696047402058 | 0.0542754736303123 | parietal |

### Beta - Meditation 1 (top5)

| electrode | degree | strength | betweenness | region |
|---|---:|---:|---:|---|
| CP1 | 33.50 | 0.385097377999183 | 0.0417306707629288 | central |
| Oz | 35.00 | 0.368978368984516 | 0.0431387608806964 | occipital |
| C1 | 30.75 | 0.341639126759444 | 0.036226318484383 | other |
| POz | 31.75 | 0.322419278392229 | 0.0783410138248848 | other |
| Iz | 30.00 | 0.312145501064619 | 0.0508192524321557 | occipital |

### Beta - Meditation 2 (top5)

| electrode | degree | strength | betweenness | region |
|---|---:|---:|---:|---|
| CP1 | 36.25 | 0.572371039404078 | 0.0801331285202253 | central |
| AF3 | 32.50 | 0.505201582157251 | 0.0472350230414747 | frontal |
| AFz | 31.00 | 0.475471792637632 | 0.0403225806451613 | other |
| C1 | 25.50 | 0.390381945834533 | 0.0311059907834101 | other |
| Oz | 23.25 | 0.333420022755894 | 0.0160010240655402 | occipital |

### Beta - Thinking (top5)

| electrode | degree | strength | betweenness | region |
|---|---:|---:|---:|---|
| POz | 39.25 | 0.407102043403761 | 0.090757808499744 | other |
| CP1 | 35.75 | 0.402633467610243 | 0.0459549411162314 | central |
| C1 | 32.00 | 0.351351012577394 | 0.0561955965181772 | other |
| Oz | 32.25 | 0.312628061519889 | 0.0614439324116743 | occipital |
| CP3 | 27.00 | 0.285995086269846 | 0.0083205325140809 | central |

---

## Visualizations

Included visualizations (embedded PNGs regenerated from CSVs):

![Top-5 hubs per condition](plots/top5_hubs_per_condition.png)

![Hub region distribution](plots/hub_region_distribution.png)

![Meditation vs Thinking strength comparison](plots/meditation_vs_thinking_strength.png)


---

