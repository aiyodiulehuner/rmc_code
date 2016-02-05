function PAV_QP=qpparams_pav_margin(Jcol)
    n=Jcol(end);d2=length(Jcol)-1;
    PAV_QP.H=[speye(n,n+d2);sparse(d2,n+d2)];
    PAV_QP.beq=1;    
    PAV_QP.Aeq=sparse(ones(d2,1),((n+1):(n+d2)),ones(d2,1));
    PAV_QP.Aineq=sparse(0,n+d2);
    njsum=0;
    for j=1:d2
        ind = (Jcol(j)+1):Jcol(j+1);
        nj=length(ind);
        if nj>0
            Dj=spdiags([ones(nj-1,1),-ones(nj-1,1)],[njsum,njsum+1],nj-1,n+d2);
            Dj=Dj+sparse((1:nj-1),(n+j)*ones(nj-1,1),ones(nj-1,1),nj-1,n+d2);
            PAV_QP.Aineq=[PAV_QP.Aineq;Dj];
            njsum=njsum+nj;
        end
    end
    PAV_QP.Aineq=[PAV_QP.Aineq;[sparse(d2,n),-speye(d2)]];
    PAV_QP.bineq=zeros(size(PAV_QP.Aineq,1),1);
    PAV_QP.lb=[];
    PAV_QP.ub=[];
    PAV_QP.solver='quadprog';
    PAV_QP.options=optimset('Display','off','Algorithm','interior-point-convex','TolFun',1e-16);
end