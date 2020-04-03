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

function S = FnSparseRecovery(X,A,TypeSubspace,Method,lambda)
T=size(X,2);
m=size(X,1);
n=size(A,2);
for i = 1:T
    x = X(:,i);
    

    % L1 optimization using CVX
    if TypeSubspace == 1
        if ( strcmp(Method , 'Lasso') )
            cvx_begin;
            cvx_precision high
            variable s(N-1,1);
            minimize( norm(s,1) + lambda * norm(A * s  - x) );
            subject to
            sum(s) == 1;
            cvx_end;
        elseif ( strcmp(Method , 'L1Perfect') )
            cvx_begin;
            cvx_precision high
            variable s(N-1,1);
            minimize( norm(s,1) );
            subject to
            A * s  == x;
            sum(s) == 1;
            cvx_end;
        elseif ( strcmp(Method , 'L1Noisy') )
            cvx_begin;
            cvx_precision high
            variable s(N-1,1);
            minimize( norm(s,1) );
            subject to
            norm( A * s  - x ) <= lambda;
            sum(s) == 1;
            cvx_end;
        elseif ( strcmp(Method , 'L1ED') )
            cvx_begin;
            cvx_precision high
            variable s(N-1+D,1);
            minimize( norm(c,1) );
            subject to
            [A eye(m)] * s == x;
            sum(s(1:T-1)) == 1;
            cvx_end;
        end
    else
        if ( strcmp(Method , 'Lasso') )
            cvx_begin;
            cvx_precision high
            variable s(n,1);
            minimize( norm(A * s  - x) );
            cvx_end;
        elseif ( strcmp(Method , 'L1Perfect') )
            cvx_begin;
            cvx_precision high
            variable s(n,1);
            minimize( norm(s,1) );
            subject to
            A * s  == x;
            cvx_end;
        elseif ( strcmp(Method , 'L1Noisy') )
            cvx_begin;
            cvx_precision high
            variable s(n,1);
            minimize( norm(s,1) );
            subject to
            norm( A * s  - x ) <= lambda;
            cvx_end;
        elseif ( strcmp(Method , 'L1ED') )
            cvx_begin;
            cvx_precision high
            variable s(T-1+m,1);
            minimize( norm(c,1) );
            subject to
            [A eye(m)] * s  == x;
            cvx_end;
        end
    end
    S(:,i)=s;
end
    