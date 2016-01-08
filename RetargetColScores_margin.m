 function [y,eps,stat]=RetargetColScores_margin(x,eps,Jcol)
    eps0=0*eps0;
    f=-x;
    %fprintf('new')
    %disp(PAV_QP);
    
    [z,stat.fval,stat.exitflag,stat.output,stat.lambda]=quadprog(PAV_QP.H(1:2500,1:2500),f,...
        PAV_QP.Aineq(1:2450,1:2500),PAV_QP.bineq(1:2450),[],[],[],[],x,PAV_QP.options);
    y=z(1:length(x));
    eps=z(length(x)+1:end);
end
