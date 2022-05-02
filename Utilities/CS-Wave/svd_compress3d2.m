function [ res, res2, comp_mtx ] = svd_compress3d2( in, in2, num_svd, flip_on )
% svd coil compression for 3d data
% assumes that coil axis is the 4th dimension

if nargin < 3
    flip_on = 0;
end

mtx_size = size(in(:,:,:,1));
mtx_size2 = size(in2(:,:,:,1));
 
temp = reshape(in, prod(mtx_size), []);
temp2 = reshape(in2, prod(mtx_size2), []);

[v,d] = eig(temp'*temp);

if flip_on
    v = flipdim(v,2);
    comp_mtx = v(:,1:num_svd);    
else
    comp_mtx = v(:,end-num_svd+1:end);
end

res = reshape(temp * comp_mtx, [mtx_size, num_svd]);
res2 = reshape(temp2 * comp_mtx, [mtx_size2, num_svd]);

end