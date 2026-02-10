#!/usr/bin/env python3
"""
generate_figures.py — Reproduce all figures from the paper.

Reads pre-computed data from ../data/ and generates:
  - figure1_bifurcation.pdf   (Sophie Germain bifurcation of S(f))
  - figure2_collapse.pdf      (Effective degree collapse)
  - figure3_accuracy.pdf      (Bateman-Horn accuracy + amplification)
  - figure4_verification.pdf  (6-panel B-H verification curves)
  - figure5_convergence.pdf   (Convergence + penalty limit)

Usage:
    python generate_figures.py
    python generate_figures.py --outdir ../figures
"""

import argparse
import csv
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from math import gcd, log
from scipy.integrate import quad
from scipy.optimize import curve_fit

# ── Helpers ──────────────────────────────────────────────────────────

def sieve(limit):
    s = [True] * (limit + 1)
    s[0] = s[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            for j in range(i*i, limit + 1, i):
                s[j] = False
    return [i for i in range(2, limit + 1) if s[i]]

PRIMES = sieve(1_400_000)[:100_000]

def li(x):
    if x <= 2: return 0
    r, _ = quad(lambda t: 1/log(t), 2, x)
    return r

def singular_series(q, primes=PRIMES):
    c = 1.0
    for p in primes:
        g = gcd(q, p - 1)
        c *= (p - g + 1) / (p - 1)
    return c

def singular_series_running(q, primes):
    c = 1.0
    result = []
    for p in primes:
        g = gcd(q, p - 1)
        c *= (p - g + 1) / (p - 1)
        result.append(c)
    return result

def is_sophie_germain(q):
    p = 2 * q + 1
    if p < 2: return False
    return all(p % d != 0 for d in range(2, int(p**0.5) + 1))

def power_law(x, a, b):
    return a * x**b

# ── Load data ────────────────────────────────────────────────────────

def load_summary(datadir):
    """Load summary_table.csv and return structured dict."""
    data = {}
    with open(os.path.join(datadir, 'summary_table.csv')) as f:
        reader = csv.DictReader(f)
        for row in reader:
            q = int(row['q'])
            data[q] = {
                'is_sg': row['is_sophie_germain'] == 'True',
                'd': int(row['degree']),
                'S': float(row['singular_series']),
                'pred': int(row['BH_predicted']),
                'actual': int(row['actual_count']),
                'ratio': float(row['ratio']),
                'd_eff': float(row['d_eff']),
            }
    return data

def load_counts(datadir):
    """Load prime_counts_by_N.csv."""
    Ns = []
    counts = {}
    with open(os.path.join(datadir, 'prime_counts_by_N.csv')) as f:
        reader = csv.reader(f)
        header = next(reader)
        q_vals = [int(h.replace('q=','')) for h in header[1:]]
        for q in q_vals:
            counts[q] = []
        for row in reader:
            Ns.append(int(row[0]))
            for i, q in enumerate(q_vals):
                counts[q].append(int(row[i+1]))
    return np.array(Ns), counts, q_vals

# ── Style ────────────────────────────────────────────────────────────

def setup_style():
    plt.rcParams.update({
        'font.size': 11, 'font.family': 'serif',
        'mathtext.fontset': 'dejavuserif',
        'axes.labelsize': 13, 'axes.titlesize': 14,
        'legend.fontsize': 10, 'xtick.labelsize': 10,
        'ytick.labelsize': 10, 'axes.linewidth': 1.0,
        'figure.dpi': 300, 'savefig.dpi': 300,
        'axes.grid': True, 'grid.alpha': 0.2, 'grid.linestyle': '--',
    })

# ── Figure 1: Bifurcation ───────────────────────────────────────────

def make_figure1(data, outdir):
    q_vals = sorted(data.keys())
    sg_q = [q for q in q_vals if data[q]['is_sg']]
    nsg_q = [q for q in q_vals if not data[q]['is_sg']]

    nsg_d = np.array([data[q]['d'] for q in nsg_q])
    nsg_S = np.array([data[q]['S'] for q in nsg_q])
    sg_d = np.array([data[q]['d'] for q in sg_q])
    sg_S = np.array([data[q]['S'] for q in sg_q])

    p_nsg, _ = curve_fit(power_law, nsg_d, nsg_S, p0=[1, 0.5])
    p_sg, _ = curve_fit(power_law, sg_d[sg_d > 3], sg_S[sg_d > 3], p0=[1, 0.3])
    d_range = np.linspace(2, 180, 200)

    fig, ax = plt.subplots(figsize=(11, 7.5))
    ax.plot(d_range, power_law(d_range, *p_nsg), '-', color='#5DADE2', lw=2, alpha=0.5)
    ax.plot(d_range, power_law(d_range, *p_sg), '-', color='#F1948A', lw=2, alpha=0.5)
    ax.fill_between(d_range, power_law(d_range, *p_nsg)*0.88,
                    power_law(d_range, *p_nsg)*1.12, color='#5DADE2', alpha=0.05)
    ax.fill_between(d_range, power_law(d_range, *p_sg)*0.85,
                    power_law(d_range, *p_sg)*1.15, color='#F1948A', alpha=0.05)
    ax.scatter(nsg_d, nsg_S, c='#2471A3', s=130, marker='D',
              edgecolors='#1B4F72', lw=1.2, zorder=5)
    ax.scatter(sg_d, sg_S, c='#E74C3C', s=130, marker='o',
              edgecolors='#922B21', lw=1.2, zorder=5)
    ax.axhline(y=1.0, color='#AAA', ls=':', lw=1.2, alpha=0.5)

    label_cfg = {3:(-3,-0.8,8), 5:(-3,-0.8,8), 7:(3,0.5,8), 11:(3,-0.8,8),
                 13:(3,0.5,8), 19:(3,0.5,9), 23:(3,-0.8,9), 47:(3,0.5,9),
                 83:(3,-0.8,9), 167:(3,0.5,9.5)}
    for q in q_vals:
        d = data[q]['d']; S = data[q]['S']; is_sg = data[q]['is_sg']
        color = '#922B21' if is_sg else '#1B4F72'
        if q in label_cfg:
            dx, dy, fs = label_cfg[q]
            ax.annotate(f'q={q}', (d, S), xytext=(d+dx, S+dy),
                       fontsize=fs, color=color, fontweight='bold',
                       arrowprops=dict(arrowstyle='-', color=color, lw=0.6, alpha=0.4)
                       if abs(dx) > 4 else None)

    # SG penalty arrow for q=83
    q = 83; d = data[q]['d']; S = data[q]['S']
    np_ = power_law(d, *p_nsg)
    ax.annotate('', xy=(d, S), xytext=(d, np_),
               arrowprops=dict(arrowstyle='->', color='#C0392B', lw=2, ls='--'))
    ax.text(d+4, (S+np_)/2, 'SG penalty\n ×0.51', fontsize=9,
            color='#C0392B', fontweight='bold', va='center')

    ax.set_xlabel('Degree  d = q − 1', fontsize=14)
    ax.set_ylabel(r'Singular Series  $\mathfrak{S}(f)$', fontsize=14)
    ax.set_title('Figure 1:  The Sophie Germain Bifurcation',
                fontsize=16, fontweight='bold', pad=12)
    ax.set_xlim(-5, 180); ax.set_ylim(0, 16)
    handles = [
        plt.Line2D([], [], marker='D', color='#2471A3', markersize=8,
                   markeredgecolor='#1B4F72', ls='',
                   label=r'Non-SG: $\mathfrak{S} \approx %.1f \cdot d^{%.2f}$' % (p_nsg[0], p_nsg[1])),
        plt.Line2D([], [], marker='o', color='#E74C3C', markersize=8,
                   markeredgecolor='#922B21', ls='',
                   label=r'SG: $\mathfrak{S} \approx %.1f \cdot d^{%.2f}$' % (p_sg[0], p_sg[1]))
    ]
    ax.legend(handles=handles, fontsize=11, loc='upper left', framealpha=0.9)
    fig.tight_layout(pad=1.5)
    fig.savefig(os.path.join(outdir, 'figure1_bifurcation.pdf'), bbox_inches='tight')
    plt.close()
    print("  ✓ figure1_bifurcation.pdf")

# ── Figure 2: Collapse ──────────────────────────────────────────────

def make_figure2(data, outdir):
    q_vals = sorted(data.keys())
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.plot([0, 18], [0, 18], '--', color='#CCC', lw=1.2, zorder=0)

    for track, marker, c, ec, label in [
        (False, 'D', '#2471A3', '#1B4F72', 'Non-SG (shielded)'),
        (True,  'o', '#E74C3C', '#922B21', 'Sophie Germain (penalized)')]:
        ds = [data[q]['d'] for q in q_vals if data[q]['is_sg'] == track]
        des = [data[q]['d_eff'] for q in q_vals if data[q]['is_sg'] == track]
        ax.scatter(ds, des, c=c, s=140, marker=marker, edgecolors=ec,
                  lw=1.2, zorder=5, label=label)

    main_labels = {31:(5,0.3), 37:(5,-0.8), 41:(5,0.4), 43:(5,-0.8),
                   47:(5,0.4), 53:(5,0.4), 61:(5,-0.8), 71:(5,0.4),
                   83:(-5,1.8), 167:(-8,1.5)}
    for q in q_vals:
        if q in main_labels:
            d = data[q]['d']; de = data[q]['d_eff']; is_sg = data[q]['is_sg']
            color = '#922B21' if is_sg else '#1B4F72'
            xo, yo = main_labels[q]
            ax.annotate(f'q={q}', (d, de), xytext=(d+xo, de+yo),
                       fontsize=9, color=color, fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.15',
                                facecolor='#FADBD8' if is_sg else '#D6EAF8',
                                edgecolor='#ccc', alpha=0.85),
                       arrowprops=dict(arrowstyle='->', color=color, lw=0.8)
                       if q in [83, 167] else None)

    # Inset for q=3..23
    axins = ax.inset_axes([0.55, 0.05, 0.43, 0.50])
    axins.set_xlim(-1.5, 24); axins.set_ylim(-0.5, 6.2)
    axins.set_title('Zoom: q = 3 to 23', fontsize=9, fontweight='bold', color='#555')
    axins.grid(True, alpha=0.15, ls='--')
    axins.plot([0, 6.2], [0, 6.2], '--', color='#CCC', lw=0.8, zorder=0)

    for q in q_vals:
        if q > 23: continue
        d = data[q]['d']; de = data[q]['d_eff']; is_sg = data[q]['is_sg']
        m = 'o' if is_sg else 'D'
        c = '#E74C3C' if is_sg else '#2471A3'
        ec = '#922B21' if is_sg else '#1B4F72'
        axins.scatter([d], [de], c=c, s=70, marker=m, edgecolors=ec, lw=1, zorder=5)

    ins_cfg = {3:(1.2,-0.50), 5:(-3.8,0.25), 7:(1.2,0.25), 11:(-5.0,-0.45),
               13:(1.2,-0.45), 17:(-5.0,0.25), 19:(1.2,0.25), 23:(1.2,-0.45)}
    for q in [3,5,7,11,13,17,19,23]:
        d = data[q]['d']; de = data[q]['d_eff']; is_sg = data[q]['is_sg']
        color = '#922B21' if is_sg else '#1B4F72'
        xo, yo = ins_cfg[q]
        axins.annotate(f'q={q}', (d, de), xytext=(d+xo, de+yo),
                      fontsize=7.5, color=color, fontweight='bold',
                      arrowprops=dict(arrowstyle='-', color=color, lw=0.3, alpha=0.4))
    axins.set_xlabel('d', fontsize=8)
    axins.set_ylabel(r'$d_{\mathrm{eff}}$', fontsize=8)
    axins.tick_params(labelsize=7)

    rect = Rectangle((-1.5, -0.5), 25.5, 6.7, lw=1, edgecolor='#888',
                     facecolor='#F8F9FA', alpha=0.2, ls='--', zorder=1)
    ax.add_patch(rect)

    ax.set_xlabel('Physical Degree  d = q − 1', fontsize=14)
    ax.set_ylabel(r'Effective Degree  $d_{\mathrm{eff}} = d\, /\, \mathfrak{S}(f)$', fontsize=14)
    ax.set_title('Figure 2:  Effective Degree Collapse', fontsize=16, fontweight='bold', pad=10)
    ax.set_xlim(-5, 180); ax.set_ylim(-1, 21)
    ax.legend(fontsize=11, loc='upper left', framealpha=0.9)
    fig.savefig(os.path.join(outdir, 'figure2_collapse.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ figure2_collapse.pdf")

# ── Figure 3: Accuracy + Amplification ──────────────────────────────

def make_figure3(data, outdir):
    q_vals = sorted(data.keys())
    LI = li(1e8)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    x = np.arange(len(q_vals))
    colors = ['#E74C3C' if data[q]['is_sg'] else '#2471A3' for q in q_vals]
    ratios = [data[q]['ratio'] for q in q_vals]

    ax1.bar(x, ratios, color=colors, edgecolor='#555', lw=0.4, alpha=0.85)
    ax1.axhline(y=1.0, color='#333', ls='-', lw=1.2)
    ax1.fill_between([-0.5, len(q_vals)-0.5], 0.99, 1.01, color='#2ECC71', alpha=0.12)
    ax1.set_xticks(x); ax1.set_xticklabels([str(q) for q in q_vals], fontsize=7, rotation=45)
    ax1.set_ylabel('Actual / Predicted', fontsize=12); ax1.set_ylim(0.96, 1.015)
    ax1.set_title('Figure 3a: Bateman–Horn Accuracy', fontweight='bold', fontsize=13)
    for i, r in enumerate(ratios):
        ax1.text(i, r+0.0015, f'{r:.3f}', ha='center', fontsize=5.5, fontweight='bold')
    ax1.legend(handles=[Patch(facecolor='#2471A3', label='Non-SG'),
                        Patch(facecolor='#E74C3C', label='SG')], fontsize=9)

    naive = np.array([LI/(q-1) for q in q_vals])
    amp = np.array([data[q]['actual'] for q in q_vals]) / naive
    ax2.bar(x, amp, color=colors, edgecolor='#555', lw=0.4, alpha=0.85)
    ax2.scatter(x, [data[q]['S'] for q in q_vals], color='#F39C12', s=35, zorder=5,
               marker='*', label='Theoretical S(f)')
    ax2.set_xticks(x); ax2.set_xticklabels([str(q) for q in q_vals], fontsize=7, rotation=45)
    ax2.set_ylabel('Amplification vs generic', fontsize=12)
    ax2.set_title('Figure 3b: Amplification over Generic', fontweight='bold', fontsize=13)
    for i in range(len(q_vals)):
        ax2.text(i, amp[i]+0.25, f'{amp[i]:.1f}×', ha='center', fontsize=5.5, fontweight='bold')
    ax2.legend(fontsize=9)

    fig.tight_layout(pad=1.5)
    fig.savefig(os.path.join(outdir, 'figure3_accuracy.pdf'), bbox_inches='tight')
    plt.close()
    print("  ✓ figure3_accuracy.pdf")

# ── Figure 4: 6-panel verification ──────────────────────────────────

def make_figure4(data, Ns, counts, outdir):
    selected = [3, 23, 47, 83, 71, 167]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5))
    axes = axes.flatten()

    for idx, q in enumerate(selected):
        ax = axes[idx]
        d = data[q]['d']; S = data[q]['S']; is_sg = data[q]['is_sg']
        color = '#E74C3C' if is_sg else '#2471A3'
        pi_v = np.array(counts[q])
        bh_v = np.array([S * li(n) / d for n in Ns])
        na_v = np.array([li(n) / d for n in Ns])
        ratio = data[q]['ratio']

        ax.plot(Ns, pi_v, '-', color=color, lw=2.5, label='Actual', zorder=3)
        ax.plot(Ns, bh_v, '--', color='#27AE60', lw=2, label=f'B-H (S={S:.2f})')
        ax.plot(Ns, na_v, ':', color='#999', lw=1.5, label='Generic (S=1)')
        ax.fill_between(Ns, na_v, bh_v, alpha=0.08, color='#27AE60')

        sg_tag = " [SG]" if is_sg else ""
        ax.set_title(f'q={q}{sg_tag}, d={d}, Ratio={ratio:.4f}',
                    fontweight='bold', fontsize=9.5)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.legend(fontsize=6.5, loc='upper left')
        if idx >= 3: ax.set_xlabel('N', fontsize=10)
        if idx % 3 == 0: ax.set_ylabel('Prime count', fontsize=10)

    fig.suptitle('Figure 4: Bateman–Horn Verification (6 selected exponents)',
                fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(outdir, 'figure4_verification.pdf'), bbox_inches='tight')
    plt.close()
    print("  ✓ figure4_verification.pdf")

# ── Figure 5: Convergence + Penalty ─────────────────────────────────

def make_figure5(data, outdir):
    q_vals = sorted(data.keys())
    sg_q = [q for q in q_vals if data[q]['is_sg']]
    cc = ['#DC143C', '#E67E22', '#2471A3', '#9B59B6', '#1B2A4A']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    for i, q in enumerate([3, 23, 47, 83, 167]):
        run = singular_series_running(q, PRIMES[:2000])
        ax1.plot(range(1, len(run)+1), run, '-', lw=1.5, color=cc[i],
                label=f'q={q}' + (' [SG]' if q in [3,23,83] else ''))
        for j, p in enumerate(PRIMES[:2000]):
            if (p - 1) % q == 0:
                ax1.plot(j+1, run[j], 'v', color=cc[i], ms=7, zorder=6)
                break
    ax1.set_xlabel('Number of primes in product', fontsize=11)
    ax1.set_ylabel('Running S(f)', fontsize=11)
    ax1.set_title('Figure 5a: Singular Series Convergence', fontweight='bold', fontsize=12)
    ax1.set_xscale('log'); ax1.legend(fontsize=9)

    sg_ext = [2,3,5,11,23,29,41,53,83,89,113,131,173,179,191,233,239,251,281,293]
    ax2.plot(sg_ext, [(q+2)/(2*q) for q in sg_ext], 'o-', color='#E74C3C', ms=4, lw=1.2)
    ax2.axhline(y=0.5, color='#333', ls='--', lw=1.5, label='Limit: 1/2')
    for q in sg_q:
        if q in sg_ext:
            ax2.scatter([q], [(q+2)/(2*q)], c='#E74C3C', s=100,
                       edgecolors='black', lw=1.5, zorder=6)
            ax2.annotate(f'q={q}', (q, (q+2)/(2*q)), textcoords="offset points",
                        xytext=(8, 4), fontsize=8, fontweight='bold')
    ax2.set_xlabel('Sophie Germain prime q', fontsize=11)
    ax2.set_ylabel('Penalty factor (q+2)/(2q)', fontsize=11)
    ax2.set_title('Figure 5b: SG Penalty → 1/2', fontweight='bold', fontsize=12)
    ax2.legend(fontsize=10)

    fig.tight_layout(pad=1.5)
    fig.savefig(os.path.join(outdir, 'figure5_convergence.pdf'), bbox_inches='tight')
    plt.close()
    print("  ✓ figure5_convergence.pdf")

# ── Main ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Generate all figures for the paper')
    parser.add_argument('--outdir', type=str, default='../figures',
                        help='Output directory (default: ../figures)')
    parser.add_argument('--datadir', type=str, default='../data',
                        help='Data directory (default: ../data)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    setup_style()

    print("Loading data...")
    data = load_summary(args.datadir)
    Ns, counts, q_vals = load_counts(args.datadir)

    print(f"Generating figures in {args.outdir}/")
    make_figure1(data, args.outdir)
    make_figure2(data, args.outdir)
    make_figure3(data, args.outdir)
    make_figure4(data, Ns, counts, args.outdir)
    make_figure5(data, args.outdir)
    print("\nAll figures generated.")

if __name__ == '__main__':
    main()
