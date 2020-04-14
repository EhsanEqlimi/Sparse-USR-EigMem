%**************************************************************************
% Main and demo m-file to run EigMem Underdetermined Sparse Source Recovery
%**************************************************************************
% Let X = AS: Underdetermined Blind Source Seperation (UBSS) Problem
% X= mixture matrix (m*T): known
% A= mixing matrix (m*n) : known and singular
% S= source matrix (n*T) : unknown
% n>m: more source than sensors (underdetermined BSS)
% S matrix is sparse: each column of S is a k-sparse vector, i.e., there is
% k<=(m-1) active sources at each time instant. In other words, sparsity
% level (k) or norm-zero of source vectors (s_q) is the number of non-zero
% elemnets. For example, assume S=[2 0 3 5;3 0 0 -1], S(:,1) is 3-sparse
% and S(:,2) is 2-sparse vectors.
% Aim: Find  source matrix, S--> Underdetermined Source Recovery (USR)
%**************************************************************************
% Witten by Ehsan Eqlimi, @TUMS, Tehran, Iran
% Copyright @ Ehsan Eqlimi and Bahador Makkiabadi
% Department of Medical Physics and Biomedical Enfineering,
% Tehran University of Medical Sciences (TUMS), Tehran, Iran
% Date: May 2016 - Jan 2017
% E-mail: Ehsan.Eqlimi@ugent.be, Ehsun.Eghlimi@gmail.com,
%Eghlimi@razi.tums.ac.ir
%**************************************************************************
% Reference:
% [1] E. Eqlimi, B. Makkiabadi, N. Samadzadehaghdam, H. Khajehpour,
% F. Mohagheghian, and S. Sanei, “A novel underdetermined source
% recovery algorithm based on k-sparse component analysis,” Circuits,
% Systems, and Signal Processing, vol. 38, no. 3, pp. 1264–1286, 2019.

%Please cite the above paper (and future papers) in case of any usage or
% benchmarking.
%%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%
clc
clear
close all
%% Preliminaries
m=3; % The number of sensors/channels a.k.a the ambient dimension
n=5; % The number of sources
k=m-1; % The number of acitive source in each time point (nonzero elements
% in each column of S) or the dimension of each subspace
r=m-k; % The dimension of the normal space (normal vector/orthogonal complement) of each subspace
c=nchoosek(n,k); % The number of all possible r-dim subspaces
T= 500; % The Number of data poins a.k.a the number of samples (time points)
Sigma=0;%1e-4 % Parameter controls the standard deviation of normal noise over the
% zero sources (non-strict saprsity) a.k.a source noise!
AMode=1; % k-SCA Condition for A are satisfied
IterNum=500;% The number of iteration to generate a good singular mixing matrix
RankTh=0.1; % to generate a good singular mixing matrix
Orth=0; %if Orth=1 A is orth (not applicable!)
EigMem=1; % Switch for EigMem method ((not applicable!))
SNRIn=[300,60,40,30]; % SNRs in dB (additive noise)
%% SCA Mixing (Simulted Sources)
% Mixing Mode:
% 1- PermkSCA: Uniform and Permuted SCA (if Sigma!=0 => Noisy PermkSCA)
% 2- kSCA:Uniform and Organized SCA (if Sigma!=0 => Noisy kSCA)
% 3- MSCA(Multiple SCA): Uniform and noisy SCA but 0<k<m.
% 4-kSCANoisy=kSCA in noisy case
MixingMode='PermkSCA';% {'PermkSCA'kSCA';kSCANoisy;'MSCA'};
N=zeros(1,c);
N(1,:)=ceil(T/c); % N is the number of the subspaces in k-SCA mixing mode
Nk=zeros(c,k);
for j=1:k
    for i=1:nchoosek(n,j)
        Nk(i,j)=ceil(ceil(T/k)/nchoosek(n,j));% Nk showes the number of each subspace in MSCA mode
    end
end
%% Sparse Component Mixing, X=AS+E simulation based k-SCA
[X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,Sigma,IterNum,RankTh,AMode,MixingMode);
% additive noise
for j=1:length(SNRIn)
    for i=1:size(X,1)
        SumSquareX=sum(X(i,:).^2);
        [SigmaIn(j,i),VarIn(j,i)]=FnSNR2Sigma(SumSquareX,SNRIn(j),size(X,2));
        Noise{j}(i,:)=SigmaIn(j,i).*randn(1,size(X,2));
        XNoisy{j}(i,:)=X(i,:)+ Noise{j}(i,:);
        SNRIn_2(j,i)=10*log10(SumSquareX/(sum( Noise{j}(i,:).^2)));
    end
end
%% EigMem Underdetermined Soure Recovery based on k-SCA assumptions
for i=1:length(SNRIn)
    Xin=XNoisy{i};
    % EigMem Algorithm (Alg. 1 in the paper) for all modes except MSCA
    if ~strcmp(MixingMode,'MSCA')
        [EigMemSubspaceInds,EigMemLabelCluster]=FnEigMembershipSubspaceClustering(Xin,A,k,MixingMode);
        % Clustering error
        Error_in_Clustering(i) = sum(Labels(:) ~= EigMemLabelCluster(:)) / length(Labels);
        disp(['Clustering error for input SNR_' num2str(SNRIn(i)) '=' num2str(Error_in_Clustering(i))])
        
        % MSCA mode: in progress.... Please wait for future versions
    else
        [EigMemSubspaceInds,EigMemLabelCluster]=FnEigMembershipSubspaceClusteringMSCA(Xin,A,k);
    end
    %% Underdtermined Source Recovery for all modes except MSCA
    SHatEigMem=FnSparseSourceRecovery(X,EigMemSubspaceInds,A);
    % In very noisy case, you may use this function based on MMSE
    % SHatEigMem_MMSE=FnSparseSourceRecovery_MMSE(X,EigSubspaceInds,A,mean(Var(i,:)));
    %% Evaluation, Source Seperation SNR
    for q=1:n
        SNREigMem(q)=FnSNRCalc(S(q,:),SHatEigMem(q,:));
    end
    MeanSNR(q)=mean(SNREigMem);
    disp(['Sepearion SNR for input SNR_'  num2str(SNRIn(i)) '=' num2str(MeanSNR(q))])
end







