% Sparse Componnet Mixing function for simulation: FnSparseComponentMixing
%% Inputs:
% m: The number of sensors a.k.a the ambient dimension
% n: The number of sources
% k: The number of acitive source in each time point (nonzero elements
% in each col. of S) or the dimension of each subspace
% N: The cardinal of each subspace in kSCA mode
% Nk: the cardinal of  each subspace in MSCA mode
% Orth:%if Orth=1 A is orth
% MixingMode: Sparse Component Mixing Mode:
%'kSCA'=> Uniform SCA (k is same for all columns)
%'MSCA'=> Multiple SCA (0<k<m)
%'PermkSCA'=> Permuted and uniform SCA (same as 'kSCA' but permuted and captures the source noise)
%'kSCANoisy'=> Permuted and uniform SCA (same as 'kSCA' but captures the source noise)

% Sigma: we called Var^S and Sigma^S in the following paper
% i.e., the standard deviation over inactive elemnts of each source vector
% and refers to the noise during source inactive state, i.e. how much inactive components are
% really zero. Indeed, the k active sources in each time instant are a zero-mean Gaussian
% random variable with unit variance and n ?k inactive sources in each time instant are
% a zero-mean Gaussian random variable with variance Sigma^S << 1.

%IterNum,RankTh,AMode: If Amode==1 then we try to find the best
%ill-conditioned and singular mixing matrix based on k_rank_EE function
%% Outputs:
% X: The mixtures / mixed data/observed data/data points : m*T
% S: The source matrix/The sparse componenets/The weights in linaer
% combination of the columns of A a basis of the subspaces:n*T
% A: The mixing matrix/ the basis of subspaces/ The dictionary: m*n
% Labels: ground-truth for subspace clustering:1*T,Label(i) is a number
% between 1 and the total number of the subspaces
% SubspaceInds: each row of P shows the indices of the columns of A that
% has generated each subspace.P is c*k in kSCA mode which c is the total
% number of subspace and is a cell in MSCA mode.
% SubspaceNum: The total number of the subspaces.
% *************************************************************************
% Ehsan Eqlimi, 26 Bahman 1394/ Winter 2016
% Edit 1: Add a parameter to control the variance of gaussian noise over the
% inactive(nonzeros) sources: 9 Xordad 1394 and Add PermkSCA and kSCANoisy
% *************************************************************************
% Witten by Ehsan Eqlimi, @TUMS, Tehran, Iran
% Copyright @ Ehsan Eqlimi and Bahador Makkiabadi
% Department of Medical Physics and Biomedical Enfineering,
% Tehran University of Medical Sciences (TUMS), Tehran, Iran
% Date: May 2016 - Jan 2017
% E-mail: Ehsan.Eqlimi@ugent.be, Ehsun.Eghlimi@gmail.com,
% Eghlimi@razi.tums.ac.ir
%**************************************************************************
% Reference:
% [1] E. Eqlimi, B. Makkiabadi, N. Samadzadehaghdam, H. Khajehpour,
% F. Mohagheghian, and S. Sanei, �A novel underdetermined source
% recovery algorithm based on k-sparse component analysis,� Circuits,
% Systems, and Signal Processing, vol. 38, no. 3, pp. 1264�1286, 2019.

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

function [X,S,OrthA,A,Labels,SubspaceInds,SubspaceNum]=FnSparseComponentMixing(m,n,k,T,N,Nk,Sigma,IterNum,RankTh,AMode,MixingMode)
if AMode
    BestKa=0;
    Iter=0;
    while BestKa<RankTh && Iter<IterNum
        A=FnColNormalizer(randn(m,n));
        [R,Ka]=k_rank_EE(A);
        if Ka>BestKa
            BestKa=Ka;
            BestA=A;
        end
        Iter=Iter+1;
    end
    A=BestA;
else
    A=FnColNormalizer(randn(m,n));
end
P=nchoosek(1:n,k);
c=length(P);
OrthA=[];
X=[];
Labels=[];
S=[];
% % % A=[0.5525 0.3919 0.5707 0.3934 0.6904;0.6863 -0.6066 0.5166 0.8634
% -0.6007;0.4730 0.6917 -0.6383 -0.3158 0.4032];%Whashizava example
switch MixingMode
    case 'kSCA' %
        SubspaceInds=P;
        for i=1:c
            %%a=orth(A(:,P(i,:)));
            a=(A(:,P(i,:)));
            OrthA=[OrthA a];
            s=randn(k,N(i));
            sWithZeros=zeros(n,N(i));
            sWithZeros(P(i,:),:)=s;
            IdxActive = zeros(n,1);
            IdxActive(P(i,:))=1;
            IdxActive=logical(IdxActive);
            sWithZeros(~IdxActive,:)=Sigma*randn(n-k,N(i));
            S=[S sWithZeros];
            x=a*s;
            X=[X x];
            % S=[S s];
            Labels=[Labels i*ones(1,N(i))];
        end
        
    case 'PermkSCA'
        SubspaceInds=P;
        S=zeros(n,T);
        for i=1:T
            IdxActive = zeros(n,1);
            %IdxRandPerm = randperm(k);
            Labels(i) = randperm(size(P,1),1);
            IdxRandPerm=P(Labels(i),:);
            IdxActive(IdxRandPerm) = 1;
            IdxActive = logical(IdxActive);
            S(IdxActive,i) = randn(k,1);
            
            S(~IdxActive,i) = Sigma*randn(n-k,1);
        end
        X=A*S;
        
    case 'kSCANoisy'
        SubspaceInds=P;
        for i=1:c
            %a=orth(A(:,P(i,:)));
            %a=(A(:,P(i,:)));
            %OrthA=[OrthA a];
            s=randn(k,N(i));
            sWithZeros=zeros(n,N(i));
            sWithZeros(P(i,:),:)=s;
            IdxActive = zeros(n,1);
            IdxActive(P(i,:))=1;
            IdxActive=logical(IdxActive);
            sWithZeros(~IdxActive,:)=Sigma*randn(n-k,N(i));
            S=[S sWithZeros];
            %= x=a*s;
            %X=[X x];
            %S=[S s];
            Labels=[Labels i*ones(1,N(i))];
        end
        %for j=1:n
        %S(j,:)=S(j,:)-mean(S(j,:));
        %end
        
        X=A*S;
    case 'MSCA'
        Labels=0;
        for j=1:k
            Added(j)=max(Labels);
            kj=j;
            Pj=nchoosek(1:n,kj);
            SubspaceInds{j}=Pj;
            cj=length(Pj);
            for i=1:cj
                a=orth(A(:,Pj(i,:)));
                OrthA=[OrthA a];
                s=randn(kj,Nk(i,j));
                sWithZeros=zeros(n,Nk(i,j));
                sWithZeros(Pj(i,:),:)=s;
                S=[S sWithZeros];
                x=a*s;
                X=[X x];
                %S=[S s];
                Labels=[Labels (i+Added(j))*ones(1,Nk(i,j))];
            end
        end
        Labels=Labels(2:end);
        
    case 'GMM' %Equation 8 of Movahedi's paper (in progress! :-))
        p=k/n;
        z=rand(n,T);
        S=(z>p).*randn(n,T)*std_noise+(z<=p).*(randn(n,T)*std_source);
end
SubspaceNum=max(Labels);