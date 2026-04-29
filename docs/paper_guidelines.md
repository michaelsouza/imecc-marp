# Optimization Paper Checklist

A structured guideline for evaluating the quality of scientific papers in the optimization field.

> **Tags:** `MUST` = dealbreaker if missing · `GOOD` = strong quality signal · `WATCH` = red flag · `BONUS` = exceptional

---

## 01 · Problem Definition

- [ ] The optimization problem is formally stated (objective function, constraints, variables) `MUST`
- [ ] The problem class is identified (convex, non-convex, combinatorial, stochastic, etc.) `MUST`
- [ ] Assumptions are made explicit (smoothness, Lipschitz continuity, convexity, etc.) `MUST`
- [ ] The scope and limitations of the formulation are acknowledged `GOOD`

---

## 02 · Contributions & Motivation

- [ ] Contributions are listed explicitly in the Introduction `MUST`
- [ ] The gap in existing literature is clearly identified `MUST`
- [ ] The motivation is backed by real-world applicability or theoretical need `GOOD`
- [ ] Claims in the abstract match the actual results in the paper `WATCH`

---

## 03 · Proposed Method

- [ ] The algorithm is presented as formal pseudocode `MUST`
- [ ] Intuition behind the method is explained before the formal derivation `GOOD`
- [ ] Computational complexity (time and space) is analyzed `GOOD`
- [ ] Hyperparameters and their sensitivity are discussed `GOOD`
- [ ] The method is compared conceptually to related approaches `GOOD`

---

## 04 · Theoretical Analysis

- [ ] Convergence guarantees are proven (or absence is justified) `MUST`
- [ ] Convergence rate is stated (e.g., O(1/k), O(1/k²), linear, superlinear) `GOOD`
- [ ] Proofs are self-contained or clearly reference supporting lemmas `GOOD`
- [ ] Optimality conditions or bounds are derived `BONUS`
- [ ] Theoretical results are connected back to practical implications `BONUS`

---

## 05 · Experiments & Validation

- [ ] Benchmarks use standard, well-known datasets or problem instances `MUST`
- [ ] Baselines are competitive and fairly tuned (not strawmen) `MUST`
- [ ] Ablation studies isolate the effect of each component `GOOD`
- [ ] Statistical significance is reported (mean ± std over multiple runs) `GOOD`
- [ ] Scalability experiments are included (varying problem size) `GOOD`
- [ ] Wall-clock time or computational cost is reported `GOOD`
- [ ] Failure cases or limitations are discussed honestly `WATCH`

---

## 06 · Reproducibility

- [ ] All hyperparameters are reported in full detail `MUST`
- [ ] Code is publicly available (GitHub, supplementary material, etc.) `GOOD`
- [ ] Random seeds and initialization strategies are specified `GOOD`
- [ ] Hardware and software environment are described `BONUS`

---

## 07 · Writing & Presentation

- [ ] Notation is consistent and defined before use `MUST`
- [ ] Figures and tables are self-contained with descriptive captions `GOOD`
- [ ] Related work is comprehensive and fairly cited `GOOD`
- [ ] The paper follows a recognized venue template (NeurIPS, ICML, IEEE, etc.) `GOOD`
- [ ] Conclusions accurately reflect results without overclaiming `WATCH`

---

## Scoring Guide

| Checked (%) | Verdict |
|---|---|
| < 30% | ⚠️ Weak paper — significant gaps in quality |
| 30–54% | 🔍 Below average — proceed with caution |
| 55–71% | 📊 Acceptable — solid but room for improvement |
| 72–87% | ✅ Good paper — well-structured and credible |
| ≥ 88% | ⭐ Excellent — meets high scientific standards |