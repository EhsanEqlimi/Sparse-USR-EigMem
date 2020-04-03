% Normalize every cloumn of X to [0,1] range
%==========================================================================
% Syntax: X=FnNormalizer(X)
%==========================================================================
% Input variable: 
% 1)X : Input Matrix or column vector.The columns are observation and the
% rows are 
%==========================================================================
% Output variable:
% 1)X: The matrix after normalization
%==========================================================================
% Ehsan Eqlimi (Ehsun.Eghlimi@gmail.com)
% Date: 20-October-2014
%==========================================================================
function X=FnColNormalizer(X)

for i=1:size(X,2)
    X(:,i)=X(:,i)/norm(X(:,i));
end