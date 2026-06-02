# Collection Report

Collection subagent: `Dirac`, model override `gpt-5.5`, reasoning effort `xhigh`, service tier `priority`.

The subagent inspected the repository read-only and did not edit files. It reported that `handoff-pro/prompt-000/context` already contains the right core bundle and that the copied sources match the current repository sources.

## Included Core Bundle

| Original path | Included path | Collection finding |
|---|---|---|
| `paper/proximal_bcd_graph_realization.tex` | `context/paper/proximal_bcd_graph_realization.tex` | Authoritative paper draft. Contains the DGP objective, proximal convergence framework, and confirmed fragile Powell-type cycle. |
| `codes/powell_type_dgp_cycle.py` | `context/codes/powell_type_dgp_cycle.py` | Exact rational verifier for the confirmed fragile cycle. It verifies a legal tied cycle, not stability. |
| `review/review-01.md` | `context/review/review-01.md` | Full write-up of the stronger \(K=1\), \(K_5\) candidate. Useful but contains the known overclaim about open basin. |
| `codes/powell_stable_dgp_candidate.py` | `context/codes/powell_stable_dgp_candidate.py` | Computational source for formal \(K_5\) checks. Does not test the later \(J^3\) branch failure and does not prove true global branch invariance. |
| `context/branch_invariance_failure.md` | `context/branch_invariance_failure.md` | Concise diagnostic note for the failed \(K_5\) branch-selection claim. |
| `context/research_status.md` | `context/research_status.md` | Compact status summary for orientation. |
| `references/powell1973search.md` | `context/references/powell1973search.md` | Markdown/OCR source for Powell 1973. Useful for the target notion of Powell-type cycling and stability. |
| `paper/references.bib` | `context/paper/references.bib` | Citation metadata. |

## Skipped Material

The subagent recommended skipping:

- `paper/proximal_bcd_graph_realization.md`, because the `.tex` appears authoritative and more current.
- Planar experiment generators, tables, CSVs, and figures, because they are not central to the Pro task.
- Broad DGP reference markdown files, unless the Pro model needs additional literature context.
- Binary PDFs, unless the manual browser submission specifically benefits from them.

## Prompt Caveats Added

The prompt explicitly states:

- The confirmed \(K=1\), \(K_4\) example is valid but fragile.
- The \(K=1\), \(K_5\) candidate has correct formal tied-cycle checks but its claimed branch-selection invariance fails at the \(J^3\) test.
- A valid answer must prove genuine global-minimizer branch selection on an open set, or else prove a rigorous obstruction for a clearly defined class.
