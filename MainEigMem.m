%**************************************************************************
% Main and demo m-file to run EigMem Underdetermined Sparse Source Recovery
%**************************************************************************
% Copyright @ Ehsan Eqlimi and Bahador Makkiabadi
% Department of Medical Physics and Biomedical Enfineering,
% Tehran University of Medical Sciences (TUMS), Tehran, Iran
% Date: May 2016 - Jan 2017
% E-mail: Ehsan.Eqlimi@ugent.be, Ehsun.Eghlimi@gmail.com,Eghlimi@razi.tums.ac.ir
%**************************************************************************
% Reference:
% [1] E. Eqlimi, B. Makkiabadi, N. Samadzadehaghdam, H. Khajehpour,
% F. Mohagheghian, and S. Sanei, “A novel underdetermined source
% recovery algorithm based on k-sparse component analysis,” Circuits,
% Systems, and Signal Processing, vol. 38, no. 3, pp. 1264–1286, 2019.

%Please cite the above papers in case of any usage or benchmarking.
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
Var=1e-5; % Parameter controls the variance of normal noise over the
% zero sources (non-strict saprsity) a.k.a source noise!
AMode=1; % k-SCA Condition for A are satisfied
IterNum=500;% The number of iteration to generate a good mixing matrix
RankTh=0.1; % to generate a good mixing matrix
Orth=0; %if Orth=1 A is orth
EigMem=1; % Switch for EigMem
SNRIn=[40,50,60];
%% SCA Mixing (Simulted Sources
% Mixing Mode:
% 1- PermkSCA: Uniform and Permuted SCA (if Var!=0 => Noisy PermkSCA)
% 2- kSCA:Uniform and Organized SCA (if Var!=0 => Noisy kSCA)
% 3- MSCA(Multiple SCA):Uniform and noisy SCA
MixingMode='PermkSCA';% {'PermkSCA'kSCA';'MSCA'};
N=zeros(1,c);
N(1,:)=ceil(T/c); % N is the number of the subspaces in k-SCA mixing mode
Nk=zeros(c,k);
for j=1:k
    for i=1:nchoosek(n,j)
        Nk(i,j)=ceil(ceil(T/k)/nchoosek(n,j));% Nk showes the number of each subspace in MSCA mode
    end
end
%% Sparse Component Mixing
[X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,Var,IterNum,RankTh,AMode,Orth,MixingMode);
for j=1:length(SNRIn)
    for i=1:size(X,1)
        SumSquareX=sum(X(i,:).^2);
        [SigmaIn(j,i),VarIn(j,i)]=FnSNR2Sigma(SumSquareX,SNRIn(j),size(X,2));
        Noise{j}(i,:)=SigmaIn(j,i).*randn(1,size(X,2));
        XNoisy{j}(i,:)=X(i,:)+ Noise{j}(i,:);
        SNRIn_2(j,i)=10*log10(SumSquareX/(sum( Noise{j}(i,:).^2)));
    end
end
%% EigMem Underdetermined Soure Recovery based on k-SCA assumtions
for i=1:length(SNRIn)
    Xin=XNoisy{i};
    [EigSubspaceInds,EigLabelCluster]=FnEigMembershipSubspaceClustering(Xin,A,k,MixingMode);
    Error_in_Clustering = sum(Labels(:) ~= EigLabelCluster(:)) / length(Labels);
    %% Underdtermined Source Recovery
    SHatEigMem=FnSparseSourceRecovery(X,EigSubspaceInds,A);
    % Noisy case
    %     SHatEigMem_MMSE=FnSparseSourceRecovery_MMSE(X,EigSubspaceInds,A,mean(Var(i,:)));
    %% Source Seperation SNR
    for i=1:n
        SNREigMem(i)=FnSNRCalc(S(i,:),SHatEigMem(i,:));
    end
    MeanSNR=mean(SNREigMem)
    
    % %     %% Source Seperation SNR
    % %     for i=1:n
    % %         SNREigMem_MMSE(i)=FnSNRCalc(S(i,:),SHatEigMem_MMSE(i,:));
    % %     end
    % %     MeanSNR_MMSE=mean(SNREigMem_MMSE)
    
end







