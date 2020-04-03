
% Is written by Ehsan Eqlimi. 18 Azar 95, 09 Dec 2016 for Multiple SCA case
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

function [EigSubspaceInds,EstimatedLabelCluster,EigenVectors,EigenValues]=FnEigMembershipSubspaceClustering(X,A,k,MixingMode)
[m,n]=size(A);
km=k;
k=m-1;

T=size(X,2);
SubspacesGenInds=nchoosek(1:n,k);sel=SubspacesGenInds;
% Add 05/08/2016 Ehsan
for y=1:k
    km=y;
P_km=nchoosek(1:n,km);
for lrank=1:k-1
    if k-lrank>1
        dff=sum(diff(sel(:,1:end-lrank))');
    else
        dff=(diff(sel(:,1:end-lrank))');
    end
    val=find(dff);
    cnt_lrank(lrank)=val(1);
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
        Eig1(j)=D(1);
        
        %member2(j)=abs(V(:,1)'*A(:,j));
        EigenVector(:,j)=V(:,1);
    end
    
    %[sel(idx2(1),:);sel(idx2(2),:);sel(idx2(3),:)]
    if km==m-1
        [EigenValues(t),Idx]=min(Eig1);
        % % %         EstimatedLabelCluster(t)=Idx;
        EigSubspaceInds(:,t)=SubspacesGenInds(Idx,:)';
        EigenVectors(:,t) =EigenVector(:,Idx);
        %         sel(idx,:),sel(i,:)
        %         [Val2,Idx2]=sort(Eig1);
        %     [sel(idx2(1),:) val2(1)']
    else
        lowrank=k-km;
        [val2,idx2]=sort(Eig1);
        [sel(idx2(1:cnt_lrank(lowrank)),:) val2(1:cnt_lrank(lowrank))'];
        vv=sel(idx2(1:cnt_lrank(lowrank)),:);
        [p,v]=hist(vv(:),1:n);
        [val,idx]=sort(p);
        rsel_eig_hist=idx(end-k+lowrank+1:end);
        EigSubspaceInds(:,t)=rsel_eig_hist';
        EigenVectors(:,t)=EigenVector(:,idx2(1));
        % Edit By Ehsan 05/08/2016
        %  EstimatedLabelCluster(t)=SubspacesGenInds(idx2(1),1:km);
        % Diff=P_km-repmat(EigSubspaceInds(:,t)',size(P_km,1),1);
        
        
        
        %[sel(idx2(1),:);sel(idx2(2),:);sel(idx2(3),:)]
        %         if MixingMode=='kSCA'
        %             [SubspacesGenInds(Idx2(1),:) Val2(1)'];
        %         else
        %             [SubspacesGenInds(Idx2(1:cnt_lrank(lowrank)),:) val2(1:cnt_lrank(lowrank))'];
        %
    end
    for i=1:size(P_km,1)
        if P_km(i,:)==sort(EigSubspaceInds(:,t)');
            EstimatedLabelCluster(t)=i;
        end
    end
end
end
