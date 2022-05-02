function vart = prodPSFx(vart,PSF)

    
    
    numdims = length(size(PSF));
    if numdims == 3 % no PSF per 3D coil  
        
        %for i = size(vart,4)
        %    vart(:,:,:,i) = vart(:,:,:,i) .* PSF(:,:,:);
        %end

        vart = vart.*repmat(PSF,[1,1,1,size(vart,4)]);    

    elseif numdims == 4 % PSF per 3D coil
%         for i = size(vart,4)
%             res(:,:,:,i) = vart(:,:,:,i) .* PSF(:,:,:,i);
%         end
        vart = vart.*PSF;
    else
        error('!!! Error in PSF .* data computation. Please check');
    end
    
end