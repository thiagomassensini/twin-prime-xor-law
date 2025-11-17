# Twin Prime XOR Law

Computational validation of the Twin Prime XOR Identity theorem.

## Overview

This repository contains the implementation and validation for the theorem:

For any twin prime pair (p, p+2) with p > 3:
```
p XOR (p+2) = 2^(v_2(p+1)+1) - 2
```
where v_2(n) denotes the 2-adic valuation of n.

## Key Results

- Validated over 1 billion twin prime pairs (p < 10^12)
- 100% verification rate
- Chi-squared test: 20.40 (df=14, p < 0.05)
- Confirms geometric distribution P(k) = 2^(-k)

## Contents

- `twin_prime_miner_ultimate.cpp` - High-performance twin prime generator
- `twin_prime_validator.cpp` - Validation suite with OpenMP parallelization
- `ultra_v5.cpp` - Enhanced implementation
- `unified_paper_final.tex` - Research paper LaTeX source
- `unified_paper_final.pdf` - Compiled paper

## Data Availability

- **Zenodo DOI**: [10.5281/zenodo.17629124](https://doi.org/10.5281/zenodo.17629124)
- **Complete Dataset** (11 CSV files): https://tprime.massensini.com.br/
- **OSF Project**: https://osf.io/bkgme/

## Author

Thiago Fernandes Motta Massensini Silva  
Email: thiago.massensini@gmail.com  
ORCID: [0009-0002-3415-9805](https://orcid.org/0009-0002-3415-9805)

## License

- Code: CC BY 4.0
- Paper: CC BY 4.0

## Citation

```
Silva, T. F. M. M. (2025). Binary Structure of Twin Primes and Connection 
to Iwasawa Lambda-Invariants. DOI: 10.5281/zenodo.17629124
```

## References

Hardy, G. H., & Wright, E. M. (2008). An Introduction to the Theory of Numbers (6th ed.). Oxford University Press.

Crandall, R., & Pomerance, C. (2005). Prime Numbers: A Computational Perspective (2nd ed.). Springer.
Discovery, proof and large-scale computational validation of a new binary law for twin primes. Includes a high-performance C++ miner capable of scanning up to billions of twin prime pairs and confirming the exact XOR identity  𝑝 ⊕ ( 𝑝 + 2 ) = 2 𝑣 2 ( 𝑝 + 1 ) + 1 − 2 p⊕(p+2)=2 v 2 (p+1)+1 −2.
