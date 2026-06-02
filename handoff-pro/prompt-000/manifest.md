# Handoff Manifest

Package: `handoff-pro/prompt-000`

Purpose: manual submission to a Pro Extended model for the DGP/Powell-type cycling research question.

This package prepares material only; do not automate browser submission as part of this handoff.

## Primary Prompt

| File | Role |
|---|---|
| `prompt.md` | Main prompt to paste first. It states the problem, known examples, failed attempt, required deliverables, and acceptance criteria. |

## Context Files

| Original path | Copied path | Why included | Caveat |
|---|---|---|---|
| `paper/proximal_bcd_graph_realization.tex` | `context/paper/proximal_bcd_graph_realization.tex` | Main paper draft. Defines the DGP objective, proximal BCD framework, convergence theorem, and confirmed fragile Powell-type counterexample. | Full LaTeX source; model should focus on the counterexample and BCD/DGP definitions. |
| `paper/references.bib` | `context/paper/references.bib` | Bibliographic entry for Powell 1973 and related convergence references. | Only needed for citation metadata. |
| `codes/powell_type_dgp_cycle.py` | `context/codes/powell_type_dgp_cycle.py` | Exact rational verifier for the confirmed \(K=1\), \(K_4\) fragile two-cycle. | Verifies the fragile example only; not an open-basin construction. |
| `codes/powell_stable_dgp_candidate.py` | `context/codes/powell_stable_dgp_candidate.py` | Verifies formal tied subproblem identities, nonzero gradients, branch Jacobian, and initial cone checks for the failed \(K_5\) candidate. | Does not prove invariant branch selection; should be read with `context/branch_invariance_failure.md`. |
| `references/powell1973search.md` | `context/references/powell1973search.md` | Markdown/OCR conversion of Powell's 1973 cycling paper. Useful for the target notion of Powell-type cycling and stability. | OCR has formatting/noise artifacts. The PDF exists locally but is not included because the Markdown is easier to submit. |
| `review/review-01.md` | `context/review/review-01.md` | Detailed write-up of the stronger \(K=1\), \(K_5\) candidate. | Contains an overclaim: the branch-selection cone is not invariant, so the open-basin conclusion is invalid. |
| none, generated for this handoff | `context/research_status.md` | Self-contained summary of the current state of the research and known conclusions. | Authored summary, not an independent source. |
| none, generated for this handoff | `context/branch_invariance_failure.md` | Explains why the \(K_5\) candidate's claimed open basin fails at the \(J^3\) branch-selection test. | Numeric check; sufficient to flag the overclaim, not a complete impossibility proof. |

## Notes

| File | Role |
|---|---|
| `notes/collection-report.md` | Collection phase report. |
| `notes/review-report.md` | Final review phase report and edits applied. |

## Recommended Submission Order

1. Paste or upload `prompt.md`.
2. Upload or paste `context/research_status.md`.
3. Upload or paste `context/branch_invariance_failure.md`.
4. Upload the paper source and verifier scripts if the interface allows attachments.
5. Include `context/review/review-01.md` only with the explicit caveat that its open-basin conclusion is false.
