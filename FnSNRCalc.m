% This function calcultes the seperation SNR (dB)
% Equations (9) and (10) in the
% following paper)
% the y-axis quantity of the figures in the following paper
%**************************************************************************
% Input variables: 
% 1)S : original source matrix (n*T)
% 2)Shat : estimated source matrix (n*T)
%**************************************************************************
% Output variable: 
% SNR (dB) 
%**************************************************************************
% Output variable (important ones) 
% 1) min_krank: minimum k-rank
% 2) mineig: minimum eigenvalue of m*m-1 submatrices
%**************************************************************************
% Witten by Ehsan Eqlimi, @TUMS, Tehran, Iran
% Copyright @ Ehsan Eqlimi and Bahador Makkiabadi
% Department of Medical Physics and Biomedical Enfineering,
% Tehran University of Medical Sciences (TUMS), Tehran, Iran
% Date: May 2016 - Jan 2017
% E-mail: Ehsan.Eqlimi@ugent.be, Ehsun.Eghlimi@gmail.com,Eghlimi@razi.tums.ac.ir
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
function SNR=FnSNRCalc(S,Shat)
SNR=10*log10((norm(S))^2/(norm(Shat-S))^2);