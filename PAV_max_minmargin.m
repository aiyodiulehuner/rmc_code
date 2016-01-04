function z = PAV_max_minmargin(y, c, zinit);

n = length(y);

if (size(y,1)~=1)
    assert(size(y,2)==1,'y should be a vector');
    y=y';
end

% build \tilde{eps}
epsinit=min(diff(zinit));
z = c_pav_max_minmargin_suriya(y, c, epsinit, 1000);

z=z';



% ybar = [1, 11, 8, 7, 6, 7, 13, 12, 11, 8, 9, 10, 4, 8]
% c= 10


% test
% addpath('pav')
% ybar = [1, 11, 8, 7, 6, 7, 13, 12, 11, 8, 9, 10, 4, 8]
% PAV_max_minmargin(ybar, 1, ybar)