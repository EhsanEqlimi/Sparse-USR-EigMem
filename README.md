# Sparse-USR-EigMem
Matlab codes to run EigMem Underdetermined Sparse Source Recovery (EigMem-USR)

Run MainEigMem and enjoy!!

Let X = AS: Underdetermined Blind Source Seperation (UBSS) Problem

X= mixture matrix (m*T): known

A= mixing matrix (m*n) : known and singular

S= source matrix (n*T) : unknown

n>m: more source than sensors (underdetermined BSS)

S matrix is sparse: each column of S is a k-sparse vector, i.e., there is

k<=(m-1) active sources at each time instant. In other words, sparsity
level (k) or norm-zero of source vectors (s_q) is the number of non-zero
elemnets. For example, assume S=[2 0 3 5;3 0 0 -1], S(:,1) is 3-sparse
and S(:,2) is 2-sparse vectors.

Aim: Find  source matrix, S--> Underdetermined Source Recovery (USR)

Witten by Ehsan Eqlimi, @TUMS, Tehra, Iran

Copyright @ Ehsan Eqlimi and Bahador Makkiabadi

Department of Medical Physics and Biomedical Enfineering,

Tehran University of Medical Sciences (TUMS), Tehran, Iran

Date: May 2016 - Jan 2017

E-mail: Ehsan.Eqlimi@ugent.be, Ehsun.Eghlimi@gmail.com,

Eghlimi@razi.tums.ac.ir
**************************************************************************
Reference:
[1] E. Eqlimi, B. Makkiabadi, N. Samadzadehaghdam, H. Khajehpour,
F. Mohagheghian, and S. Sanei, “A novel underdetermined source
recovery algorithm based on k-sparse component analysis,” Circuits,
Systems, and Signal Processing, vol. 38, no. 3, pp. 1264–1286, 2019.

Please cite the above paper (and future papers) in case of any usage or
benchmarking.
