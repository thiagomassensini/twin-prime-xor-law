# GitHub Copilot Instructions

---
## ⚠️ CRITICAL RULES
**NEVER USE EMOJIS IN ANY CODE, COMMENTS, OR OUTPUT**

**Author Information:**
- Name: Thiago Fernandes Motta Massensini Silva
- Email: thiago.massensini@gmail.com
- Repository: twin-prime-xor-law
---

## Project Context

This repository contains research on the binary structure of twin primes through XOR operations and their connection to elliptic curves via 2-adic valuations.

### Main Results

1. **Twin Prime XOR Identity**: For twin prime pairs (p, p+2) with p > 3:
   ```
   p ⊕ (p+2) = 2^(v₂(p+1)+1) - 2
   ```
   where v₂(n) is the 2-adic valuation (highest power of 2 dividing n).

2. **Elliptic Curve λ₂-Invariant Formula**:
   ```
   λ₂(E) = v₂(c₂)
   ```
   where c₂ is the Tamagawa number at p=2.

3. **Empirical Distribution**: k-values follow geometric law P(k) = 2^(-k)

## Code Guidelines

### Mathematical Notation
- Use `v2` or `v_2` for 2-adic valuation
- Use `xor_val` or `xor_result` for XOR operations
- Use `lambda2` for Iwasawa λ₂-invariant
- Use `k` for v₂(p+1) in twin prime contexts

### C++ Code (Performance-Critical)
When generating C++ code for prime computations:
- Use `uint64_t` for primes up to 10^12
- Use `__builtin_ctzll()` for count-trailing-zeros (fast v₂ computation)
- Prefer deterministic Miller-Rabin with witnesses {2, 3, 5, 7, 11, 13, 17}
- Use OpenMP for parallelization: `#pragma omp parallel for`
- Include GMP for arbitrary precision: `#include <gmp.h>`

Example pattern:
```cpp
// Fast 2-adic valuation
int v2(uint64_t n) {
    return __builtin_ctzll(n);
}

// XOR identity verification
uint64_t xor_identity(uint64_t p) {
    int k = v2(p + 1);
    return (1ULL << (k + 1)) - 2;
}
```

### Python/SageMath Code (Elliptic Curves)
When working with elliptic curves:
- Use SageMath's `EllipticCurve()` constructor
- Access Tamagawa numbers via `.tamagawa_number(2)`
- Use `.cremona_label()` for standard curve naming
- For BSD checks, use `.lseries().dokchitser()`

Example pattern:
```python
# Lambda2 computation
def compute_lambda2(E):
    c2 = E.tamagawa_number(2)
    return c2.valuation(2)

# Validate formula
E = EllipticCurve([0, 0, 1, -7, 6])  # 5077a1
assert compute_lambda2(E) == E.modular_degree().valuation(2)
```

### LaTeX Documentation
When generating LaTeX:
- Use `\vTwo` for v₂ operator
- Use `\xor` for XOR symbol (⊕)
- Use `\lambda_2` for Iwasawa invariant
- Reference theorems as `Theorem~\ref{thm:main}`
- Use `\eqref{eq:xor-identity}` for equation references

### Data Validation
All computational results should include:
- Sample size (e.g., "1,004,364,744 twin prime pairs")
- Range (e.g., "p < 10^12")
- Statistical tests (chi-squared with df and p-value)
- 100% verification rate confirmation

### Performance Benchmarks
When discussing performance:
- Report throughput in pairs/second or curves/second
- Include hardware specs (CPU model, cores, RAM)
- Separate parsing time from computation time
- Report validation accuracy as percentage

## File Structure

- `twin_primes_xor_clean.tex`: Twin primes focus (477 lines)
- `unified_paper_final.tex`: Combined twin primes + elliptic curves (463 lines)
- `twin_prime_xor_paper.tex`: Initial simplified version (140 lines)
- `gitignore/`: Directory for temporary files (auto-ignored by git)

## Mathematical Context

### Key Concepts
- **2-adic valuation v₂(n)**: Exponent of highest power of 2 dividing n
- **Twin primes**: Prime pairs (p, p+2)
- **XOR (⊕)**: Bitwise exclusive-or operation
- **Iwasawa λ₂-invariant**: Growth rate of 2-primary Selmer groups
- **Tamagawa number c_p**: Order of component group at prime p
- **Mordell-Weil rank**: Number of independent rational points on E(ℚ)

### Important Invariants
- XOR and λ₂ are **correlated** (ρ = +0.3744)
- λ₂ and rank are **independent** (ρ = 0.023)
- Both follow approximate geometric distributions

## Research Standards

### Citation Requirements
Always include:
- Zenodo DOI: 10.5281/zenodo.17520242 (twin primes dataset)
- GitHub: https://github.com/thiagomassensini/twin-prime-xor-law
- Cremona database for elliptic curves

### Validation Standards
- Twin primes: Deterministic Miller-Rabin primality testing
- Elliptic curves: SageMath built-in verification
- Statistical tests: Chi-squared with 95% confidence level
- Report any exceptions or failures explicitly

## Common Operations

### Twin Prime XOR Computation
```python
def twin_prime_xor(p):
    """Compute XOR of twin prime pair (p, p+2)"""
    k = (p + 1).valuation(2)  # 2-adic valuation
    return 2^(k+1) - 2
```

### Distribution Analysis
```python
def analyze_k_distribution(pairs):
    """Verify geometric distribution P(k) = 2^(-k)"""
    k_counts = Counter(k for p, k in pairs)
    chi2 = sum((obs - exp)^2 / exp for obs, exp in zip(observed, expected))
    return chi2 < 23.685  # Critical value at 95%, df=14
```

### Elliptic Curve Validation
```python
def validate_lambda2_formula(E):
    """Verify λ₂ = v₂(c₂)"""
    c2 = E.tamagawa_number(2)
    lambda2_expected = c2.valuation(2)
    # Compare with actual Iwasawa invariant computation
    return lambda2_expected
```

## Notes for Copilot

- Prioritize mathematical correctness over code brevity
- Always validate formulas against known test cases
- Include docstrings with mathematical notation
- Reference paper theorems when implementing algorithms
- Maintain consistency with existing codebase notation
- When in doubt about mathematical definitions, refer to Hardy & Wright or Silverman

## Language Preferences

- Primary: English (mathematical/code)
- Secondary: Portuguese (user communication)
- LaTeX: Standard mathematical notation
- Code comments: English with mathematical symbols
