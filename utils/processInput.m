%inputs ii,jj,yy and returns csc format: ii is indices, Jcol is indptr start at 0
%ii(Jcol(j)+1:Jcol(j+1)) is such that yy(Jcol(j)+1:Jcol(j+1)) is sorted min:max

function [yys,iis,Jcol]=processInput(ii,jj,yy)
    assert(length(ii)==length(yy));assert(issorted(jj));
    Jcol=compJcol(jj); %csc sparse format ii is indices, jcol is indptr
    iis=ii;yys=yy;
    for j=1:(length(Jcol)-1)
        Yj=yy(Jcol(j)+1:Jcol(j+1));
        [y,IX]=sort(Yj,'ascend');
        Ij=ii(Jcol(j)+1:Jcol(j+1));
        iis(Jcol(j)+1:Jcol(j+1))=Ij(IX);
        yys(Jcol(j)+1:Jcol(j+1))=y;
    end