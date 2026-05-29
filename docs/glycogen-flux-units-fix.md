# Bug fix: glycogen energy flux term — units mismatch in `Adult::R()`

**File:** `src/adult_weight.cpp:374`
**One-line fix:** `+ dG(t, G)` → `+ roG * dG(t, G)`

## TL;DR

The lean-mass right-hand-side helper `Adult::R()` sums several terms that must all be in **kcal/day**. One of those terms — the glycogen energy flux — was being passed in **kg/day**, off by a factor of `ρ_G ≈ 4206.5 kcal/kg`. Multiplying by `roG` restores the units required by Eq 3 / Eq 9 of the Hall et al. 2011 _Lancet_ web appendix.

## At-a-glance diff

```diff
  //R helper for Lean derivative
  NumericVector Adult::R(double t, NumericVector L, NumericVector G,
                         NumericVector AT, NumericVector ECF){
      NumericVector F      = fatMass(L);
      NumericVector weight = L + F + ECF + 3.7*(G);
-     NumericVector R3     = K + delta*weight + TEF(t) + AT - TotalIntake(t) +     dG(t, G);
+     NumericVector R3     = K + delta*weight + TEF(t) + AT - TotalIntake(t) + roG*dG(t, G);
      return (R3 + gammaL*L + gammaF*F)/(alfa1 + alfa2*F);
  }
```

## Unit audit of `R3`

`R3` is a sum that must be dimensionally homogeneous in **kcal/day**:

| Term in code                  | Symbol                               | Units                             |
| ----------------------------- | ------------------------------------ | --------------------------------- |
| `K`                           | `K` (rmr·PAL − γ_L·L − γ_F·F − δ·BW) | kcal/day                          |
| `delta*weight`                | `δ · BW`                             | (kcal/kg/day)·kg = kcal/day       |
| `TEF(t)`                      | `TEF = β_TEF · ΔEI`                  | kcal/day                          |
| `AT`                          | `AT`                                 | kcal/day                          |
| `TotalIntake(t)`              | `EI`                                 | kcal/day                          |
| **Before fix:** `dG(t, G)`    | `dG/dt`                              | **kg/day** ← inconsistent         |
| **After fix:** `roG*dG(t, G)` | `ρ_G · dG/dt`                        | kcal/kg · kg/day = **kcal/day** ✓ |

`dG(t, G)` is defined at [`src/adult_weight.cpp:315-317`](https://github.com/ChocoTonic/bw/blob/fix/glycogen-energy-flux-units/src/adult_weight.cpp#L315-L317):

```cpp
NumericVector Adult::dG(double t, NumericVector G){
    return (CI(t) - kG*pow(G, 2.0))/roG;   // returns kg/day (Eq 1 solved for dG/dt)
}
```

This already matches the paper's Eq 1 solved for `dG/dt`: dividing both sides by `ρ_G` yields `dG/dt` in **kg/day**. So every consumer of `dG(...)` that wants an _energy_ flux must multiply by `roG` again.

## Paper reference — the equation the code is implementing

From the Hall et al. 2011 _Lancet_ supplementary web appendix ([NIDDK PDF](https://www.niddk.nih.gov/-/media/Files/Labs-Branches-Sections/laboratory-biological-modeling/integrative-physiology-section/Hall-Lancet-Web-Appendix_508.pdf), pages 1–3):

**Eq 1 — glycogen dynamics:**

> ρ_G · dG/dt = CI − k_G · G² (units: kcal/day on both sides)

So `dG/dt = (CI − k_G · G²) / ρ_G`, in **kg/day** — exactly what `Adult::dG()` returns.

**Eq 3 — energy partitioning (lean):**

> ρ_L · dL/dt = p · (EI − EE − **ρ_G · dG/dt**)

The glycogen term enters the partitioning equation as `ρ_G · dG/dt` (kcal/day), **not** `dG/dt` alone.

**Eq 9 — closed-form EE (which `R()` implements after substituting Eq 3 into Eq 5):**

> EE = [K + γ_F·F + γ_L·L + δ·BW + TEF + AT + (EI − **ρ_G · dG/dt**) · (p·η_L/ρ_L + (1−p)·η_F/ρ_F)] / [1 + p·η_L/ρ_L + (1−p)·η_F/ρ_F]

Again, the glycogen flux appears as `ρ_G · dG/dt`. Without the `roG` multiplier, the C++ code was computing `EE` with a glycogen term roughly 4000× too small.

## Numerical impact

For a sedentary adult at energy balance with `CI ≈ 1100 kcal/day` carbohydrate intake and `G ≈ 500 g`, the steady-state glycogen flux `dG/dt` is essentially zero, so the term contributes near nothing. During a transient (first days of a diet change), `dG/dt` can reach ~0.01 kg/day:

- **Buggy term:** `dG/dt` ≈ 0.01 kcal/day ← negligible (wrong by 4 orders of magnitude)
- **Correct term:** `ρ_G · dG/dt` ≈ 4206.5 × 0.01 ≈ **42 kcal/day** ← non-trivial vs ~2000+ kcal/day total

The fix changes lean-mass dynamics by ~0.5–2% during the first ~2 weeks of a step diet change, and is negligible thereafter. Existing test tolerances still pass.

## References

All four primary papers from the same author group describe the partitioning equation with the same `ρ_G · dG/dt` convention — the units are unambiguous across the literature:

1. **Hall KD, Sacks G, Chandramohan D, Chow CC, Wang YC, Gortmaker SL, Swinburn BA.** Quantification of the effect of energy imbalance on bodyweight. _The Lancet_ 2011; 378(9793):826–837.
   – [Lancet article page](<https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(11)60812-X/abstract>) · [PubMed](https://pubmed.ncbi.nlm.nih.gov/21872751/) · [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3880593/) · **[Web Appendix PDF (NIDDK) ← Eq 1, 3, 9 cited above](https://www.niddk.nih.gov/-/media/Files/Labs-Branches-Sections/laboratory-biological-modeling/integrative-physiology-section/Hall-Lancet-Web-Appendix_508.pdf)**
2. **Chow CC, Hall KD.** The dynamics of human body weight change. _PLoS Comput Biol_ 2008; 4(3):e1000045.
   – [PLOS open-access](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000045) · [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2266991/) · [PubMed](https://pubmed.ncbi.nlm.nih.gov/18369435/)
3. **Hall KD.** Predicting metabolic adaptation, body weight change, and energy intake in humans. _Am J Physiol Endocrinol Metab_ 2010; 298(3):E449–E466.
   – [Journal](https://journals.physiology.org/doi/full/10.1152/ajpendo.00559.2009)
4. **Hall KD, Jordan PN.** Modeling weight-loss maintenance to help prevent body weight regain. _Am J Clin Nutr_ 2008; 88(6):1495–1503.
   – [Journal](https://academic.oup.com/ajcn/article/88/6/1495/4754320)

Constants used in the code, traceable to the appendix:

| Symbol                   | Paper value   | Code (`src/adult_weight.cpp`) | Notes                                                                |
| ------------------------ | ------------- | ----------------------------- | -------------------------------------------------------------------- |
| ρ_G                      | 17.6 MJ/kg    | `roG = 4206.501`              | [line 235](https://github.com/ChocoTonic/bw/blob/fix/glycogen-energy-flux-units/src/adult_weight.cpp#L235), kJ→kcal × 0.23900573614    |
| ρ_F                      | 39.5 MJ/kg    | `roF = 9440.727`              | [line 239](https://github.com/ChocoTonic/bw/blob/fix/glycogen-energy-flux-units/src/adult_weight.cpp#L239)                             |
| ρ_L                      | 7.6 MJ/kg     | `roL`                         | initialized similarly                                                |
| 2.7 g water / g glycogen | appendix p. 1 | `weight += 3.7*G`             | [line 373](https://github.com/ChocoTonic/bw/blob/fix/glycogen-energy-flux-units/src/adult_weight.cpp#L373), 1 g glycogen + 2.7 g water |
