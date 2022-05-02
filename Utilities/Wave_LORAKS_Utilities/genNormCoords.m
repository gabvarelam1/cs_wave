function [normalised_kspace_coordinates] = genNormCoords(Kspace_Matrix_Size, varargin)
%%% This function generates the actual normalised and centered positions
%%% for general pursposes. This function only requires the matrix size.

%%% gvm 2019-10-29 
% UNDER DEVELOPMENT

switch numel(Kspace_Matrix_Size)
    case 1 %1D%
    case 2 %2D%
        X = Kspace_Matrix_Size(1);
        Y = Kspace_Matrix_Size(2);
        
        % Checking for even or odd numbers
        if mod(X,2) == 0
            x_range = -ceil((X-1)/2) : (ceil((X-1)/2)-1);
        else
            x_range = -ceil((X-1)/2) : ceil((X-1)/2);
        end
        
        if mod(Y,2) == 0
            y_range = -ceil((Y-1)/2) : (ceil((Y-1)/2)-1);
        else
            y_range = -ceil((Y-1)/2) : ceil((Y-1)/2);
        end
        
        
        max_all = max([max(max(x_range)),max(y_range)]);
        
        x_range = x_range /  max_all;
        y_range = y_range /  max_all;
        
        [kx,ky] = meshgrid(y_range,x_range);
        normalised_kspace_coordinates = [kx(:) , ky(:)];
        
    case 3 %3D%
        X = Kspace_Matrix_Size(1);
        Y = Kspace_Matrix_Size(2);
        Z = Kspace_Matrix_Size(3);
        
        % Checking for even or odd numbers
        if mod(X,2) == 0
            x_range = -ceil((X-1)/2) : (ceil((X-1)/2)-1);
        else
            x_range = -ceil((X-1)/2) : ceil((X-1)/2);
        end
        
        if mod(Y,2) == 0
            y_range = -ceil((Y-1)/2) : (ceil((Y-1)/2)-1);
        else
            y_range = -ceil((Y-1)/2) : ceil((Y-1)/2);
        end
        
        if mod(Z,2) == 0
            z_range = -ceil((Z-1)/2) : (ceil((Z-1)/2)-1);
        else
            z_range = -ceil((Z-1)/2) : ceil((Z-1)/2);
        end
        
        
        max_all = max([max(max(x_range)),max(y_range),max(z_range)]);
        
        x_range = x_range /  max_all;
        y_range = y_range /  max_all;
        z_range = z_range /  max_all;
        
        [kx,ky,kz] = meshgrid(y_range,x_range, z_range);
        normalised_kspace_coordinates = [kx(:) , ky(:) , kz(:)];
end






end