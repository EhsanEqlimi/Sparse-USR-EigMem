% Is written by Ehsan Eqlimi. 18 Azar 95, 09 Dec 2016 for Multiple SCA case
% This function developes Algorithm 1 in paper
% Copyright @ Ehsan Eqlimi and Bahador Makkiabadi
% Department of Medical Physics and Biomedical Enfineering,
% Tehran University of Medical Sciences (TUMS), Tehran, Iran
% FunNetMap Group
% Date: May 2016 - Jan 2017
% E-mail: Eghlimi@razi.tums.ac.ir
%**************************************************************************
% Reference:
% Paper: A New Underdetermined Source Recovery Algorithm  based
% on $k$-Sparse Component Analysis, E.Eqlimi and B.Makkiabadi
%%
function [EigMemSubspaceInds,EigMemLabelCluster]=FnEigMembershipSubspaceClustering(X,A,k,MixingMode)
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
SubspacesGenInds=nchoosek(1:n,k); % This is "I" in paper
% sel=SubspacesGenInds;
% Add 05/08/2016 Ehsan
% P_km is the indices to generate the subspaces when we know about the
% number of active source in each instant i.e. uniform SCA and
% not Mutiple SCA mode.
P_km=nchoosek(1:n,km);
% Here we find the number of nondisjoint number(k-nondisjoint SCA)

for k_i=1:k-1 % the number of non disjoint subspaces is k=m-1-1 (beacuse k=m-1 is not consiedred)
    if k-k_i>1
        Diff=sum(diff(SubspacesGenInds(:,1:end-k_i))');
    else
        Diff=(diff(SubspacesGenInds(:,1:end-k_i))');
    end
    Value=find(Diff);
    Num_k_NDisjoint(k-k_i)=Value(1); %Num_k_NDisjoint is "h" in paper
end
for t=1:T
    x=X(:,t);
    for j=1:size(SubspacesGenInds,1)
        F=[x A(:,SubspacesGenInds(j,:))];
        Cov=F*F';
        [V,D]=eig(Cov);
        % % %             %% Or: SVD
        % % %             [U,S,V]=svd(F); % S is sqrt of D
        %%
        Eig1(j)=D(1); % "v" in paper
        EigenVector(:,j)=V(:,1);
    end
    [Value,Index]=sort(Eig1);
    if strcmp(MixingMode,'PermkSCA') || strcmp(MixingMode,'kSCA')
        if km==k
            EigMemLabelCluster(t)=Index(1);
            EigMemSubspaceInds(:,t)= SubspacesGenInds(Index(1),:);
        else
            Inds=SubspacesGenInds(Index(1:Num_k_NDisjoint(km)),:); % perhaps selcet m-k_i(and not m-1-k) 0.6.012017,Ehsan
            [P,V]=hist(Inds(:),1:n);
            [Val,Idx]=sort(P,'descend');
            EigSubspaceInds(:,t)=Idx(1:km);
            for i=1:size(P_km,1)
                if P_km(i,:)==sort(EigSubspaceInds(:,t)');
                    EigMemLabelCluster(t)=i;
                    EigMemSubspaceInds(:,t)=P_km(i,:);
                end
            end
            
        end  
    end
end
