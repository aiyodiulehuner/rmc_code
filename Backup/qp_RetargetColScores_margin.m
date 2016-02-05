<<<<<<< HEAD
 function [y,eps,stat]=RetargetColScores_margin(x,eps0,PAV_QP)
    
    PAV_QP.f=-[x;zeros(length(eps0),1)];%/(norm(x)^2);
    PAV_QP.H=PAV_QP.H;%/(norm(x)^2);
    %disp(PAV_QP);
    %PAV_QP.x0=[x;eps0];
    [z,stat.fval,stat.exitflag,stat.output,stat.lambda]=quadprog(PAV_QP);
=======
 function [y,eps,stat]=RetargetColScores_margin(x,eps,Jcol)
    eps0=0*eps0;
    f=-x;
    %fprintf('new')
    %disp(PAV_QP);
    
    [z,stat.fval,stat.exitflag,stat.output,stat.lambda]=quadprog(PAV_QP.H(1:2500,1:2500),f,...
        PAV_QP.Aineq(1:2450,1:2500),PAV_QP.bineq(1:2450),[],[],[],[],x,PAV_QP.options);
>>>>>>> 0b55ce12e45ebfa688c25c3e8066084b116cb4a2
    y=z(1:length(x));
    eps=z(length(x)+1:end);
end
