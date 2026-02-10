# Defying the Degree Barrier

**Arithmetic Shielding and Bimodal Effective Degree in the Titan Polynomial Family**

> Ruqing Chen — GUT Geoservice Inc., Montreal  
> Contact: ruqing@hotmail.com

---

## Overview

This repository contains the data, source code, and figures accompanying the paper:

**"Defying the Degree Barrier: Arithmetic Shielding and Bimodal Effective Degree in the Titan Polynomial Family"**

We study the Titan polynomial family Q_q(n) = n^q − (n−1)^q for prime q, and demonstrate that:

1. **Unified root count**: ω_q(p) = gcd(q, p−1) − 1 for all primes p (including p = q, by Fermat's Little Theorem).
2. **Arithmetic shielding**: The singular series S(f) grows logarithmically with q, systematically counteracting the 1/d degree barrier.
3. **Sophie Germain bifurcation**: S(f) splits into two tracks — a "shielded" track (2q+1 composite) and a "penalized" track (2q+1 prime, penalty ≈ 1/2).
4. **Bateman–Horn fit**: All 18 datasets match the heuristic prediction to within 0.3–3.1%.

## Key Result

| q | Degree d | S(f) | d_eff | Amplification |
|---|----------|------|-------|---------------|
| 3 (SG) | 2 | 3.36 | 0.59 | 3.3× |
| 47 | 46 | 8.68 | 5.30 | 8.6× |
| 83 (SG) | 82 | 4.94 | 16.60 | 4.9× |
| **167** | **166** | **14.24** | **11.65** | **14.2×** |

A degree-166 polynomial produces **14.2× more primes** than a generic polynomial of the same degree.

## Repository Structure

```
titan-polynomial-prime-shielding/
├── README.md
├── LICENSE
├── requirements.txt
├── paper/
│   ├── paper.tex                    # LaTeX source (arXiv-ready)
│   ├── paper.pdf                    # Compiled PDF
│   ├── figure1_bifurcation.png      # Figure for LaTeX compilation
│   └── figure2_collapse.png         # Figure for LaTeX compilation
├── data/
│   ├── summary_table.csv            # Main results table (18 rows)
│   ├── prime_counts_by_N.csv        # π_Q(N) at 50 values of N per q
│   ├── singular_series_components.csv # S(f) decomposition
│   └── penalty_table.csv            # Sophie Germain penalty factors
├── figures/
│   ├── figure1_bifurcation.pdf      # Sophie Germain bifurcation
│   ├── figure2_collapse.pdf         # Effective degree collapse
│   ├── figure3_accuracy.pdf         # Bateman–Horn accuracy + amplification
│   ├── figure4_verification.pdf     # 6-panel B-H verification
│   └── figure5_convergence.pdf      # Convergence + penalty limit
└── scripts/
    ├── compute_primes.py            # Generate prime values of Q_q(n)
    ├── compute_singular_series.py   # Compute S(f) via Euler product
    ├── generate_figures.py          # Reproduce all 5 figures
    └── verify_bateman_horn.py       # Full B-H verification pipeline
```

## Quick Start

```bash
# Install dependencies
pip install numpy scipy matplotlib

# Compute singular series for all 18 exponents
python scripts/compute_singular_series.py

# Generate all figures
python scripts/generate_figures.py

# Compute primes for a single exponent (e.g., q=47, N=10^6)
python scripts/compute_primes.py --q 47 --N 1000000

# Full Bateman–Horn verification
python scripts/verify_bateman_horn.py
```

> **Note**: Full computation to N = 10^8 requires significant time. The pre-computed summary data in `data/` allows reproduction of all figures and analyses without rerunning the prime search.

## Data Description

### `data/summary_table.csv`
The central results table with 18 rows (one per prime exponent q):
- `q`: prime exponent
- `is_sophie_germain`: whether 2q+1 is prime
- `degree`: d = q − 1
- `singular_series`: S(f), computed via Euler product over first 100,000 primes
- `BH_predicted`: Bateman–Horn predicted count = S(f) · Li(10^8) / d
- `actual_count`: verified prime count π_Q(10^8)
- `ratio`: actual / predicted
- `d_eff`: effective degree = d / S(f)

### `data/prime_counts_by_N.csv`
Prime counts π_Q(N) sampled at 50 logarithmically spaced values of N from 10^4 to 10^8, for each of the 18 exponents. Used to generate Figures 4 and 5a.

> **Full prime lists**: The raw prime lists (~44 million primes, ~200 MB) are not included due to size. Use `scripts/compute_primes.py` to regenerate them, or contact the author.

## Mathematical Core

**Theorem 1 (Unified Root Count).** For any prime q and any prime p:

    ω_q(p) = gcd(q, p−1) − 1

**Proof sketch.** The congruence n^q ≡ (n−1)^q (mod p) reduces to u^q ≡ 1 (mod p) where u = n(n−1)^{−1}. This has gcd(q, p−1) solutions in (Z/pZ)*. Excluding u = 1 (impossible) gives the result. For p = q, Fermat's Little Theorem gives Q_q(n) ≡ 1 (mod q) for all n, so ω = 0.

**Sophie Germain Penalty.** When 2q+1 is prime, the first obstruction at p = 2q+1 contributes factor (q+2)/(2q) → 1/2.

## Citation

If you use this data or code, please cite:

```bibtex
@article{chen2025titan,
  title={Defying the Degree Barrier: Arithmetic Shielding and Bimodal 
         Effective Degree in the Titan Polynomial Family},
  author={Chen, Ruqing},
  year={2025},
  note={arXiv preprint}
}
```

## License

MIT License. See [LICENSE](LICENSE).
