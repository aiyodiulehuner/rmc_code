 function [y,eps,stat]=RetargetColScores_margin(x,eps0,PAV_QP)
    eps0=0*eps0;
    PAV_QP.f=-[x;eps0];
    %disp(PAV_QP);
    %PAV_QP.x0=[x;eps0];
    [z,stat.fval,stat.exitflag]=quadprog(PAV_QP);
    y=z(1:length(x));
    eps=z(length(x)+1:end);
end
