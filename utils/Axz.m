%%**************************************************
% compute A*z for 
% A = A_U * A_V' - Sparse_Z;
% 
% Az = matvec(z,A_U,A_V,Sparse_Z); 
%
%%**************************************************
%%
  function Az = Axz(z); 

  global X spZ
  
  Az = X.U * (X.V' * z)+spZ * z;