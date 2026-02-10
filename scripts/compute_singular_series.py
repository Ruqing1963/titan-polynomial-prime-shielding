#!/usr/bin/env python3
"""
compute_singular_series.py — Compute the singular series S(f) for Titan polynomials

The singular series is:
    S(f) = prod_p [ (1 - w(p)/p) / (1 - 1/p) ]

where w(p) = gcd(q, p-1) - 1 by Theorem 1.

Usage:
    python compute_singular_series.py
    python compute_singular_series.py --q 167
    python compute_singular_series.py --num_primes 200000
"""

import argparse
from math import gcd, log

def sieve_of_eratosthenes(limit):
    """Return list of primes up to limit."""
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return [i for i in range(2, limit + 1) if is_prime[i]]

def compute_singular_series(q, primes):
    """
    Compute S(Q_q) using Euler product over given primes.
    
    For each prime p:
      w(p) = gcd(q, p-1) - 1
      local_factor = (1 - w(p)/p) / (1 - 1/p) = (p - w(p)) / (p - 1)
    
    Returns (S, boost_product, obstruction_product, details)
    """
    S = 1.0
    boost = 1.0
    obstruction = 1.0
    first_obstruction = None
    
    for p in primes:
        g = gcd(q, p - 1)
        w = g - 1
        factor = (p - w) / (p - 1)
        S *= factor
        
        if w == 0:
            boost *= factor  # = p/(p-1) > 1
        else:
            obstruction *= factor  # < 1
            if first_obstruction is None:
                first_obstruction = (p, w, factor)
    
    return S, boost, obstruction, first_obstruction

def is_sophie_germain(q):
    """Check if q is a Sophie Germain prime (2q+1 also prime)."""
    p = 2 * q + 1
    if p < 2:
        return False
    for d in range(2, int(p**0.5) + 1):
        if p % d == 0:
            return False
    return True

def main():
    parser = argparse.ArgumentParser(description='Compute singular series for Titan polynomials')
    parser.add_argument('--q', type=int, nargs='*', default=None,
                        help='Prime exponent(s). Default: all 18 from the paper')
    parser.add_argument('--num_primes', type=int, default=100000,
                        help='Number of primes in Euler product (default: 100000)')
    args = parser.parse_args()

    # Default: all 18 exponents from the paper
    q_list = args.q if args.q else [3,5,7,11,13,17,19,23,31,37,41,43,47,53,61,71,83,167]
    
    # Generate primes
    # Need enough primes; nth prime ~ n*ln(n), so for 100k primes need ~1.4M
    limit = max(args.num_primes * 15, 2000000)
    primes = sieve_of_eratosthenes(limit)[:args.num_primes]
    print(f"Using {len(primes):,} primes (up to {primes[-1]:,}) in Euler product\n")
    
    # Header
    print(f"{'q':>4} {'SG?':>4} {'d':>4} {'S(f)':>10} {'d_eff':>8} {'Boost':>10} "
          f"{'Obstr':>10} {'1st_obstr':>12} {'Penalty':>8}")
    print("-" * 85)
    
    for q in q_list:
        d = q - 1
        sg = is_sophie_germain(q)
        S, B, O, fo = compute_singular_series(q, primes)
        d_eff = d / S
        
        fo_str = f"p={fo[0]}" if fo else "none"
        pen_str = f"{(q+2)/(2*q):.3f}" if sg else "-"
        
        print(f"{q:>4} {'Yes' if sg else 'No':>4} {d:>4} {S:>10.4f} {d_eff:>8.2f} "
              f"{B:>10.4f} {O:>10.4f} {fo_str:>12} {pen_str:>8}")
    
    print()
    print("Legend:")
    print("  S(f)  = singular series (Euler product)")
    print("  d_eff = effective degree = d / S(f)")
    print("  Boost = product of boost factors p/(p-1) for shielded primes")
    print("  Obstr = product of obstruction factors for p ≡ 1 (mod q)")
    print("  Penalty = (q+2)/(2q) for Sophie Germain primes")

if __name__ == '__main__':
    main()
