function [y,ii]=block_sort(y,ii,blk)

for i=1:length(blk)
    ind=blk{i}(1):blk{i}(2);
    [yt,ix]=sort(y(ind),'ascend');
    y(ind)=yt;
    it=ii(ind);
    ii(ind)=it(ix);
end