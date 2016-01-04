%%**************************************************
% compute A*z for 
% A = A_U * A_V' - Sparse_Z;
% 
% Atz = matvec(z,A_U,A_V,Sparse_Z); 
%
%%**************************************************

  function Atz = Atxz(z); 

  global X spZ
  
  Atz = X.V * (X.U' * z) + (z' * spZ)';
