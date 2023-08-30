
# Sparse-USR-EigMem

Sparse-USR-EigMem is a collection of MATLAB codes designed to perform EigMem Underdetermined Sparse Source Recovery (EigMem-USR). It addresses the challenge of underdetermined blind source separation (UBSS) by recovering sparse source signals in an underdetermined system.

## Overview

In the UBSS problem, we have the following components:

- **X**: Mixture matrix (m*T) - Known
- **A**: Mixing matrix (m*n) - Known and singular
- **S**: Source matrix (n*T) - Unknown
- **n > m**: More sources than sensors (underdetermined BSS)
- **S matrix is sparse**: Each column of S is a k-sparse vector, meaning that there are k <= (m-1) active sources at each time instant. In other words, the sparsity level (k) or the norm-zero of source vectors (s_q) is the number of non-zero elements. For example, if S=[2 0 3 5;3 0 0 -1], S(:,1) is a 3-sparse vector, and S(:,2) is a 2-sparse vector.

## Aim

The primary goal is to find the source matrix, S, which corresponds to Underdetermined Source Recovery (USR).

## Usage

To use this code, follow these steps:

1. Clone this repository to your local machine.
2. Open MATLAB.
3. Run the "MainEigMem.m" script included in this repository.

## Authors

- Ehsan Eqlimi
- Bahador Makkiabadi

## Copyright

This code is copyright protected by Ehsan Eqlimi and Bahador Makkiabadi, who are affiliated with the Department of Medical Physics and Biomedical Engineering at Tehran University of Medical Sciences (TUMS) in Tehran, Iran.

## Date

This code was developed between May 2016 and January 2017.

## Contact

For any inquiries or questions, you can contact the authors at the following email addresses:

- Ehsan.Eqlimi@ugent.be
- Ehsun.Eghlimi@gmail.com
- Eghlimi@razi.tums.ac.ir

## Reference

Please cite the following paper (and future papers) if you use or benchmark this code:

- [E. Eqlimi, B. Makkiabadi, N. Samadzadehaghdam, H. Khajehpour, F. Mohagheghian, and S. Sanei, “A novel underdetermined source recovery algorithm based on k-sparse component analysis,” Circuits, Systems, and Signal Processing, vol. 38, no. 3, pp. 1264–1286, 2019.](link_to_paper)
```

Replace "link_to_paper" with the actual link to your reference paper. This updated README provides a clear and structured overview of your project and its usage.
