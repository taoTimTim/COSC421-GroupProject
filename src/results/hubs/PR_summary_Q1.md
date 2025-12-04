**PR Title:** Research Question 1 - Hub electrode analysis (branch `Q1`)

**Author / Assignee:** Aaditya (assigned)

**TL;DR**
- Adds R analysis and reporting for Research Question 1 (hub electrode identification).
- `R/hub_analysis.R` now averages hub metrics across 10/15/20/25% density thresholds.
- Reporting helpers live in `src/results/hubs/{regen_plots.R,build_report_pdf.R,hub_report.md}`.
- `.gitignore` ignores generated outputs (CSV, PDF, PNG, plots/top5); generated files remain on disk locally but stay untracked.

**What changed (high level)**
- New R analysis: `R/hub_analysis.R` loads averaged wPLI CSVs, thresholds each network at 10/15/20/25% density, averages degree/strength/betweenness across those levels, and exports combined CSVs plus per-condition top-5 lists.
- Reporting: `src/results/hubs/regen_plots.R` (generate PNGs from CSVs) and `src/results/hubs/build_report_pdf.R` (assemble `hub_report_combined.pdf` from markdown and images).
- Report source: `src/results/hubs/hub_report.md` (methods, tables, interpretation) now documents the multi-density approach.
- Repo housekeeping: `.gitignore` updated to ignore generated hub outputs; generated artifacts were removed from the index and remain only locally.

**Files touched in this PR**
- Modified: `.gitignore` — ignores R session files and generated hub outputs under `src/results/hubs`.
- Added: `R/hub_analysis.R`
- Added: `src/results/hubs/regen_plots.R`
- Added: `src/results/hubs/build_report_pdf.R`
- Added: `src/results/hubs/hub_report.md`
- Added: `.vscode/launch.json`, `.vscode/settings.json` (workspace convenience settings)

**What is intentionally tracked vs ignored**
- Tracked: analysis and report source files (`R/hub_analysis.R`, `src/results/hubs/*.R`, `src/results/hubs/hub_report.md`) and helper scripts — required to reproduce outputs.
- Ignored (and untracked): generated CSVs, PNGs, PDFs, and `top5/` CSVs under `src/results/hubs`. `.gitignore` entries include:
  - `/src/results/hubs/hub_report_combined.pdf`
  - `/src/results/hubs/all_hub_metrics_multi_density.csv`
  - `/src/results/hubs/region_summary.csv`
  - `/src/results/hubs/*.png`
  - `/src/results/hubs/*.pdf`
  - `/src/results/hubs/top5/`
  - `/src/results/hubs/plots/`

**Verification notes (what I checked)**
- Confirmed branch `Q1` contains the new scripts and the `.gitignore` change.
- Regenerated hub outputs with the multi-density averaging, refreshed plots, and rebuilt the combined PDF locally using `Rscript` (R 4.5.1).

**How to reproduce regenerated outputs locally (commands)**
1. Run hub analysis to generate CSVs:
```powershell
cd "C:\Users\Aaditya\COSC421-GroupProject"
Rscript --vanilla R\hub_analysis.R
```
2. Regenerate PNGs from the CSVs (creates `src/results/hubs/plots/*.png`):
```powershell
cd "C:\Users\Aaditya\COSC421-GroupProject"
Rscript --vanilla src\results\hubs\regen_plots.R
```
3. Rebuild the combined PDF (creates `src/results/hubs/hub_report_combined.pdf`):
```powershell
cd "C:\Users\Aaditya\COSC421-GroupProject"
Rscript --vanilla src\results\hubs\build_report_pdf.R
```

**Reviewer checklist**
- [ ] Confirm `R/hub_analysis.R` meets the analysis requirements (multi-density thresholding, centrality metrics, region mapping).
- [ ] Confirm `src/results/hubs/hub_report.md` drafts are satisfactory (methods, interpretation reflect multi-density approach).
- [ ] Run the commands above to verify CSVs, PNGs, and combined PDF regenerate correctly.
- [ ] Optionally approve removal of generated artifacts from git history (if desired, we can remove them from earlier commits using git-filter-repo, but that is not necessary now).

**Next actions (recommended)**
- Merge PR `Q1` when you're ready; the branch now contains the multi-density hub analysis, reporting scripts, and excludes large generated artifacts.

**Notes / Caveats**
- `.gitignore` prevents future commits of generated outputs but does not remove them from historical commits.

----
Generated on: 2025-12-02
