 function [y,eps,stat]=RetargetColScores_margin(x,eps0,PAV_QP)
    
    PAV_QP.f=-[x;zeros(length(eps0),1)];%/(norm(x)^2);
    PAV_QP.H=PAV_QP.H;%/(norm(x)^2);
    %disp(PAV_QP);
    %PAV_QP.x0=[x;eps0];
    [z,stat.fval,stat.exitflag,stat.output,stat.lambda]=quadprog(PAV_QP);
    y=z(1:length(x));
    eps=z(length(x)+1:end);
end
