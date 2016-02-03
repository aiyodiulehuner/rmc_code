function rho=spearman_rho(x,y)
    xrank=tiedrank(x,0);
    yrank=tiedrank(y,0);
    n=length(x);nconst=n*(n^2-1)/6;
    rho=1-(sum((xrank-yrank).^2)/nconst);
end