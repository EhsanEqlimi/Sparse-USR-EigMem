% This function estimates the source matrix using EigMem
% output(SubspaceInds)
% Each source point could be recovered
% by pseudo-inverse minimum mean squared error
% This function imlements Eq (8) in the following paper for tackling noisy case
%**************************************************************************
% Input variable: 
% 1)X : mixture matrix (m*T)
% 2)A: mixing matrix (m*n)
% 3)SubspaceInds: subspaces indices of each data point (X) 
% 4) NoisePwr: The power of the additive noise (E in X=AS+E)
%**************************************************************************
% Output variable (important ones) 
% 1) Shat: estimated source matrix (n*T) based on pseudo-inverse minimum 
% mean squared error (MMSE)
%**************************************************************************
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

function SHat=FnSparseSourceRecovery_MMSE(X,SubspaceInds,A,NoisePwr)
SHat=zeros(size(A,2),size(X,2));
nR=size(X,1);
for i=1:size(X,2)
    H=A(:,SubspaceInds(:,i));
    R=X(:,i);
    SHat(SubspaceInds(:,i),i)=H'*inv(H*H'+NoisePwr*eye(nR))*R;
end
