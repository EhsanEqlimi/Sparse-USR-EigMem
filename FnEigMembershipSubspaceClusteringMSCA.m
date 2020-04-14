% This function imlements Algorithm 1 in the following paper
% Note:This function works in MSCA mode!!!
% Please see Fig.1 in the paper
% This function is main function of paper (EigMem idea) for MSCA mode
%**************************************************************************
%Input variables:
%1) X: mixture matrix (m*T)
%2) A: mixing matrix (n*T)
%3) k: sparsity level
%4) Mixing mode: not involving in accuracy
%**************************************************************************
%output variables:
%1) EigMemSubspaceInds: h^t in the following apaper: a cell 
%**************************************************************************
% Witten by Ehsan Eqlimi,@TUMS,Tehran, Iran
% Copyright Ehsan Eqlimi and Bahador Makkiabadi
% @Department of Medical Physics and Biomedical Enfineering,
% Tehran University of Medical Sciences (TUMS), Tehran, Iran
% Date: May Oct. 2014
% E-mail: Ehsan.Eqlimi@ugent.be,Ehsun.Eghlimi@gmail.com,
% Eghlimi@razi.tums.ac.ir
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

function [EigSubspaceInds,EstimatedLabelCluster,EigenVectors,EigenValues]=FnEigMembershipSubspaceClusteringMSCA(X,A,k)
% Size of mixing matrix
[m,n]=size(A);
% km is the number of the active sources in each instant if the mixing mode
% is not Multiple SCA i.e. the mixing mode are PermkSCA for example.
km=k;
% k is the maximum paossible number of active sources in each instant i.e.
% the number of mixtures minus 1.
k=m-1;
% T is the number of the mixtures.
T=size(X,2);
% SubspaceGenInds are the indices that generate the subspaces i.e. the all
% possible k columns k columns of the columns of A.
SubspacesGenInds=nchoosek(1:n,k);
% sel=SubspacesGenInds;
% Add 05/08/2016 Ehsan
% P_km is the indices to generate the subspaces when we know about the
% number of active source in each instant i.e. uniform SCA (PermkSCA and kSCANoisy) and
% not Mutiple SCA (MSCA) mode.
P_km=nchoosek(1:n,km);
% Here we find the number of nondisjoint number(k-nondisjoint SCA)
for k_i=1:k-1 % the number of non disjoint is k-1 beacuse k=1 is disjoint
    if k-k_i>1
        Diff=sum(diff(SubspacesGenInds(:,1:end-k_i))');
    else
        Diff=(diff(SubspacesGenInds(:,1:end-k_i))');
    end
    Value=find(Diff);
    Num_k_NDisjoint(k_i)=Value(1);
end

for t=1:T
    x=X(:,t);
    for j=1:size(SubspacesGenInds,1)
        F=[x A(:,SubspacesGenInds(j,:))];
        Cov=F*F';
        [V,D]=eig(Cov);
        %Or: SVD
        %[U,S,V]=svd(F); % S is sqrt of D
        %%
        Eig1(j)=D(1);
        EigenVector(:,j)=V(:,1);
    end
    [Value2,Index2]=sort(Eig1);
    
    for k_i=1:k
        if k_i==k
            Inds=[];
            Inds=SubspacesGenInds(Index2(1:size(SubspacesGenInds,1)));
        else
            Inds=SubspacesGenInds(Index2(1:Num_k_NDisjoint(k-k_i)),:); % perhaps selcet m-k_i(and not m-1-k) 0.6.012017,Ehsan
        end
        
        [p,v]=hist(Inds(:),1:n);
        [val,idx]=sort(p);
        rsel_eig_hist=idx(end-k_i:end);
        EigSubspaceInds{k_i}(:,t)=rsel_eig_hist';
        %         EigenVectors(:,t)=EigenVector(:,Index2(1));
        
    end
    % %     if ~strcmp(MixingMode,'MSCA')
    % %         EigSubspaceInds= EigSubspaceInds{km};
    % %     end
    
    %     for i=1:size(P_km,1)
    %         if P_km(i,:)==sort(EigSubspaceInds(:,t)');
    %             EstimatedLabelCluster(t)=i;
    %         end
    %     end
    EstimatedLabelCluster=[];
end
