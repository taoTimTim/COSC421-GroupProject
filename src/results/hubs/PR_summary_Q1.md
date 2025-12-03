**PR Title:** Research Question 1 — Hub electrode analysis (branch `Q1`)

**Author / Assignee:** Aaditya (assigned)

**TL;DR**
- Adds R analysis and reporting for Research Question 1 (hub electrode identification).
- Includes `R/hub_analysis.R` (analysis), `src/results/hubs/{regen_plots.R,build_report_pdf.R,hub_report.md}` (report source and helpers).
- `.gitignore` updated to stop tracking generated outputs (CSV, PDF, PNG, plots/top5); generated files were untracked and remain on disk.

**What changed (high level)**
- New R analysis: `R/hub_analysis.R` — loads averaged wPLI CSVs, thresholds at 15% density, computes per-node degree/strength/betweenness, exports CSVs and per-condition top-5 lists.
- Reporting: `src/results/hubs/regen_plots.R` (generate PNGs from CSVs) and `src/results/hubs/build_report_pdf.R` (assemble `hub_report_combined.pdf` from markdown and images).
- Report source: `src/results/hubs/hub_report.md` (methods, tables, interpretation).
- Repo housekeeping: `.gitignore` updated to ignore generated hub outputs; generated artifacts were removed from the index and a commit pushed so future runs won't re-add them.

**Files touched in this PR**
- Modified: `.gitignore` — now ignores R session files and generated hub outputs under `src/results/hubs`.
- Added: `R/hub_analysis.R`
- Added: `src/results/hubs/regen_plots.R`
- Added: `src/results/hubs/build_report_pdf.R`
- Added: `src/results/hubs/hub_report.md`
- Added: `.vscode/launch.json`, `.vscode/settings.json` (workspace convenience settings)

**What is intentionally tracked vs ignored**
- Tracked (kept in repo): analysis and report source files (`R/hub_analysis.R`, `src/results/hubs/*.R`, `src/results/hubs/hub_report.md`) and helper scripts — these are required to reproduce the outputs.
- Ignored (and now untracked): generated CSVs, PNGs, PDFs, and `top5/` CSVs under `src/results/hubs`. `.gitignore` entries:
  - `/src/results/hubs/hub_report_combined.pdf`
  - `/src/results/hubs/all_hub_metrics_15pct.csv`
  - `/src/results/hubs/region_summary.csv`
  - `/src/results/hubs/*.png`
  - `/src/results/hubs/*.pdf`
  - `/src/results/hubs/top5/`
  - `/src/results/hubs/plots/`

**Verification notes (what I checked)**
- Confirmed branch `Q1` contains the new scripts and the `.gitignore` change.
- I untracked generated outputs and pushed that commit (`cf8abd9`). Generated files still exist on disk locally.
- Attempted to run regeneration scripts from PowerShell, but `Rscript` wasn't available in PATH in this environment so I could not execute them here. The scripts are present and should run on your machine or in an environment with R installed.

**How to reproduce regenerated outputs locally (commands)**
1. Regenerate PNGs from CSVs (creates `src/results/hubs/plots/*.png`):
```powershell
cd "C:\Users\Aaditya\COSC421-GroupProject"
Rscript --vanilla src\results\hubs\regen_plots.R
```
2. Rebuild the combined PDF (creates `src/results/hubs/hub_report_combined.pdf`):
```powershell
cd "C:\Users\Aaditya\COSC421-GroupProject"
Rscript --vanilla src\results\hubs\build_report_pdf.R
```
3. If `Rscript` is not on PATH, open R and run:
```r
source("src/results/hubs/regen_plots.R")
source("src/results/hubs/build_report_pdf.R")
```

**Reviewer checklist**
- [ ] Confirm `R/hub_analysis.R` meets the analysis requirements (thresholding, centrality metrics, region mapping).
- [ ] Confirm `src/results/hubs/hub_report.md` drafts are satisfactory (methods, interpretation).
- [ ] Run the two R commands above locally to verify PNGs and combined PDF regenerate correctly.
- [ ] Optionally approve removal of generated artifacts from git history (if desired, we can remove them from earlier commits using git-filter-repo, but that is not necessary now).

**Next actions (recommended)**
- I can run the regeneration here if you add `Rscript` to PATH or allow running in the R interactive terminal. Tell me which you prefer.
- Otherwise, merge PR `Q1` when you're ready; the branch now contains the analysis and report source and excludes large generated artifacts.

**Notes / Caveats**
- The `.gitignore` change prevents future re-commits of generated outputs but does not remove them from historical commits (they remain in git history). If history cleanup is required, that must be done explicitly.

----
Generated on: 2025-12-02
