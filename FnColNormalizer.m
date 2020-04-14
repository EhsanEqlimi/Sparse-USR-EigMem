% Normalize every cloumn of matrix X to [0,1] range 
% In fact, X is normalized to unit length vector, i.e.,
% $\textbf{x}_i=\frac{\textbf{x}_i}{\left\|\textbf{x}_i\right\|}$,
% where $\left\| \sbullet \right\|$ stands for the norm $\ell_2$.
%*************************************************************************
% Syntax: X=FnNormalizer(X)
%*************************************************************************
% Input variable: 
% 1)X : Input Matrix or column vector (e.g. chennel*time matrix)
%*************************************************************************
% Output variable:
% 1)X: The matrix after normalization
%*************************************************************************
% Written by Ehsan Eqlimi @TUMS, Tehran, Iran
% Copyright @ Ehsan Eqlimi and Bahador Makkiabadi
% Department of Medical Physics and Biomedical Enfineering,
% Tehran University of Medical Sciences (TUMS), Tehran, Iran
% Date: May 2014
% E-mail: Ehsan.Eqlimi@ugent.be,Ehsun.Eghlimi@gmail.com,
% Eghlimi@razi.tums.ac.ir
%*************************************************************************
function X=FnColNormalizer(X)
for i=1:size(X,2)
    X(:,i)=X(:,i)/norm(X(:,i));
end