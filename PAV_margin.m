function z = PAV_margin(y, wt, eps, winit);

n = length(y);
if (size(y,1)~=1)
    assert(size(y,2)==1,'y should be a vector');
    y=y';
end
 
% build \tilde{eps}

tilde_eps = 0:(n-1);
tilde_eps =  -1 *eps * tilde_eps;
ybar = y + tilde_eps;

costs = zeros(1, n)  ;
costs(1) = wt;
costs(n) = -1.0 *wt;

z = c_pav_margin(ybar, costs, winit, 1000);

z = z - tilde_eps;

if (size(z,1)==1)
    z=z';
end



%ybar = [1, 11, 8, 7, 6, 7, 13, 12, 11, 8, 9, 10, 4, 8]
%costs = [0, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

%
%n = n-2;
%a = -1:n;
%a =  eps * a;
%ybar = y + a;
%
%costs = [0, ones(1, n-1)]  ;
%costs = wt * costs;
%z = c_pav_margin(ybar, costs);


% test 
% addpath('pav')
% ybar = [1, 11, 8, 7, 6, 7, 13, 12, 11, 8, 9, 10, 4, 8]
% PAV_margin(ybar, 1,1)