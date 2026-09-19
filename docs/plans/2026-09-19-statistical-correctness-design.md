# General statistical correctness repair

Scope approved by the user: repair systemic defects identified in the review, without dataset-specific rules, gene lists, thresholds or target result counts.

## Decisions

- Pass the requested adjusted-P cutoff to DESeq2 `results(alpha=...)`. Preserve missing statistical values, explicitly distinguish unavailable results, and report prefilter counts separately.
- Obtain complete GSEA results, exclude invalid statistical rows, and apply explicit raw-P/BH-adjusted-P thresholds plus an optional q-value cutoff. Save the full and significant tables separately; derive the count, preview and plot from the significant table. Record a configurable fixed seed and use serial execution for reproducibility. Display NES on the GSEA plot.
- Normalize the complete count matrix with edgeR TMM and transform to logCPM before selecting heatmap genes or variable PCA genes. This avoids fitting a dispersion trend just to visualize small matrices and uses an existing dependency. Preserve externally transformed matrices as supplied; record the actual transformation and effective library sizes. Center PCA without rescaling each gene to unit variance.
- Strengthen the Skill's general report checks: reopen tables, verify counts and plot thresholds, establish source methods and contrast directions, and calculate regression in the direction of its plotted axes. Do not bake the reviewed study into the Skill.

## Alternatives considered

Prompt-only changes would leave incorrect backend counts and data transformations intact. Automatically using VST for every count matrix would introduce dispersion fitting and failure cases for small input matrices; callers can still supply validated VST/rlog matrices as `normalized`. A general normalization fix plus explicit contracts is the narrower implementation.

## Implementation and verification

1. Repair R backends and CLI schemas/dependency declarations.
2. Update the Skill, references and version together.
3. Add synthetic statistical regression tests: direct DESeq2 agreement at different alpha values, missing statistics, full-matrix normalization before gene selection, nominal-only/invalid/no-hit GSEA results and repeated-seed reproducibility.
4. Run the existing Python/R suite, build wheel/sdist/connector, verify an isolated wheel install and CLI analyses, and visually inspect changed plots.
5. Update the local WorkBuddy installation after checks pass, verify source/runtime parity and record exact version and test evidence. Existing analysis outputs remain untouched.

References: [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#p-values-and-adjusted-p-values), [edgeR guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).
