# monad_omicskit P0-P2 Repair Design

## Scope

Fix the confirmed P0-P2 defects from the Kimi K3 review with explicit data
semantics rather than compatibility shims:

- eliminate formula injection by renaming selected data columns to safe internal
  names before constructing R formulas;
- make time-dependent ROC load its runtime dependencies;
- require an explicit gene identifier type and map SYMBOL, ENSEMBL, or ENTREZID;
- define calibration `groups` as the requested number of groups;
- require matrix type where count versus normalized semantics affect analysis;
- use Poisson GSVA only for count matrices;
- collapse duplicate raw-count rows by sum and normalized rows by highest mean;
- aggregate duplicate GSEA ENTREZ mappings deterministically;
- transform count matrices before PCA and validate heatmap annotations;
- reject non-numeric matrix cells immediately and accept legitimate low-count
  integer matrices;
- base EPV warnings on candidate predictors and keep plots inside configured
  output devices.

## API Decisions

- Enrichment requires `id_type`: `SYMBOL`, `ENSEMBL`, or `ENTREZID`.
- PCA and heatmap require `matrix_type`: `counts` or `normalized`.
- Calibration `groups` is the desired number of calibration groups; the
  implementation derives observations per group for `rms::calibrate`.
- Invalid or ambiguous inputs fail at the MCP boundary instead of being guessed.

## Verification

Add security and scientific-semantics regression tests, then run MCP smoke
tests, the full Python suite, Python compilation, R parsing, and real R execution
for installed optional analysis packages.
