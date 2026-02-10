#!/usr/bin/env python3
"""
compute_primes.py — Compute prime values of Titan polynomials Q_q(n) = n^q - (n-1)^q

Usage:
    python compute_primes.py --q 47 --N 1000000
    python compute_primes.py --q 47 --N 100000000 --output primes_q47.txt

For the full dataset (N=10^8), expect runtimes of minutes to hours per exponent
depending on hardware.
"""

import argparse
import sys
from math import gcd, isqrt
from sympy import isprime as sympy_isprime

def Q(n, q):
    """Compute Q_q(n) = n^q - (n-1)^q."""
    return n**q - (n - 1)**q

def is_prime_miller_rabin(n):
    """Deterministic Miller-Rabin for n < 3.317×10^24 (first 13 bases)."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    # Write n-1 = 2^r * d
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    # Test with deterministic bases sufficient for n < 3.317e24
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def compute_titan_primes(q, N, verbose=True):
    """Find all n in [1, N] such that Q_q(n) is prime."""
    primes = []
    report_interval = max(1, N // 100)
    for n in range(1, N + 1):
        val = Q(n, q)
        if val > 1 and is_prime_miller_rabin(val):
            primes.append(n)
        if verbose and n % report_interval == 0:
            pct = 100 * n / N
            print(f"\r  q={q}: {pct:.0f}% done, {len(primes)} primes found", end='', flush=True)
    if verbose:
        print(f"\r  q={q}: 100% done, {len(primes)} primes found")
    return primes

def main():
    parser = argparse.ArgumentParser(description='Compute prime values of Q_q(n)')
    parser.add_argument('--q', type=int, required=True, help='Prime exponent q')
    parser.add_argument('--N', type=int, default=10**6, help='Upper bound for n (default: 10^6)')
    parser.add_argument('--output', type=str, default=None, help='Output file (default: stdout summary)')
    args = parser.parse_args()

    q = args.q
    N = args.N

    # Verify q is prime
    if not sympy_isprime(q):
        print(f"Error: q={q} is not prime.")
        sys.exit(1)

    print(f"Computing primes of Q_{q}(n) for n in [1, {N:,}]...")
    print(f"  Degree: {q-1}")
    print(f"  Q_{q}(n) = n^{q} - (n-1)^{q}")
    print()

    primes_n = compute_titan_primes(q, N)

    print(f"\nResult: pi_Q({N:,}) = {len(primes_n):,}")
    print(f"  Smallest n: {primes_n[0] if primes_n else 'none'}")
    print(f"  Largest n:  {primes_n[-1] if primes_n else 'none'}")

    if args.output:
        with open(args.output, 'w') as f:
            for n in primes_n:
                f.write(f"{n}\n")
        print(f"\nPrime n-values written to {args.output}")

if __name__ == '__main__':
    main()
