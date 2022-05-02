function [Dcompressed] = decompress_data(Data, Svd_matrix)

[sx,sy,sz,~] = size(Data);
sc = size(Svd_matrix,1);
temp2 = reshape(Data,prod([sx,sy,sz]),[]);
Dcompressed = reshape(temp2 * Svd_matrix' , [sx,sy,sz,sc]);


end