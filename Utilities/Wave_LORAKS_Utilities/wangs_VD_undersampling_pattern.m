function [USP] = wangs_VD_undersampling_pattern( Matrix_Size , Undersampling_factor , Rho_b , Rho_a , varargin)
%%% This function generates a 2D or 3D undersampling pattern based on
%%% Wang's equation for undersampling which is the following:

%%% SPk = exp( - rho_a * ( sqrt(kx.^2 + ky.^2)/N ).^ rho_b );

%%% where rho_a, rho_b are parameters, kx,ky are kspace coordinates; N is
%%% the size of kspace (example in 2D).

%%% (Wang, 2018) Accelerating quantitative susceptibility imaging
%%% acquisition using compressed sensing
% According to supplementary material
% Rho_a = 1.8
% Rho_b = [3.6,7.6]

%%% In another (Wang, 2018), they have the same equation but they flip rho's positions xD

%%% UNDER DEVELOPMENT
%%% WHAT I NEED TO ADD:
%   1. Calibration lines
%   2. Partial Fourier
%   3. "Preserve" undersampling factor as correct as possible
%   4. 3D or even more
%   5. Include amount of Coils

    % built-in functions
    
    VEC = @(v_) v_(:);
    %%%%%%%%%%
    % Parameters from varargin{'Partial Fourier', A , 'Calibration Lines'}
    % Default Values are empty ...
    PF = [];
    CL = [];
    COILS = [];
    PIx = [];
    PIy = [];
    dPIx = [];
    dPIy = [];

    if ~isempty(varargin)
        valid_argnames = {'Partial Fourier','Calibration Lines','Coils',...
                          'Parallel Imaging', 'Dephased PI'};
        [args_specified,args_location] = ismember(valid_argnames,varargin(1:2:end));
        
        if args_specified(1) == true
            PF = varargin{(args_location(1)*2-1)+1};
        end
        if args_specified(2) == true
            CL = varargin{(args_location(2)*2-1)+1};
        end
        if args_specified(3) == true
            COILS = varargin{(args_location(3)*2-1)+1};
        end
        if args_specified(4) == true
            PIxy = varargin{(args_location(4)*2-1)+1};
            switch PIxy
                case '2x1'
                    PIx = 2; PIy = 1;
                case '1x2'
                    PIx = 1; PIy = 2;
                case '4x1'
                    PIx = 4; PIy = 1;
                case '1x4'
                    PIx = 1; PIy = 4;
                case '2x2'
                    PIx = 2; PIy = 2;
                case '4x4'
                    PIx = 4; PIy = 4;
                otherwise
                    disp('Not Implemented yet');
            end
            dPIx = 0;
            dPIy = 0;
            dPIxy = 2;
            if args_specified(5) == true
                % From Breuer, 2006
                if ~(PIx == 1)
                    old_PIx = PIx;
                    old_PIy = PIy;
                    temp = PIx*PIy;
                    PIx = 1;
                    PIy = temp;
                    disp(strcat('Replacing R_d1=',num2str(old_PIx),'and R_d2=',num2str(old_PIy),...
                        'to',' R_d1=',num2str(old_PIx),'and R_d2=',num2str(old_PIy),'and shifting=',num2str(dPIxy)));
                end
                 disp('Working in implementation...let me know if it works! :p');
            end
        end
    end
    %%%%%%%%%%
    % NUMBER OF DIMENSIONS TO WORK WITH KSPACE DATA
    if numel(Matrix_Size) == 2
        flag_operation = 1;
        
        disp('Under-sampling for 2D k-space data...');
        
        SPk_nD = @(kspace_ , rho_a_ , rho_b_) exp( - rho_a_ * ...
        (sqrt( kspace_(:,:,1).^2 + kspace_(:,:,2).^2) /  2 ).^ rho_b_ );
        %sqrt(prod([size(kspace_,1),size(kspace_,2)]))
        pdf_nD = @(matrix_ ,mu_ , sigma_) mu_ + sigma_ .* (rand(size(matrix_(:,:,1)))-0.5);
    else
        flag_operation = 0;
        disp(strcat('number of dimensions from kspace data is',num2str(ndims(Kspace_data)),'and it is not implemented yet...'));
        
    end
    
    if flag_operation == 1
        %%% DO VARIABLE DENSITY FUNCTION
        
        % built-in functions
        
        VEC = @(a_) a_(:);
        
        kspace = genNormCoords(Matrix_Size);
        kspace = cat(3,reshape(kspace(:,1),Matrix_Size),reshape(kspace(:,2),Matrix_Size));
%         set_range_odd = @(a_) -ceil((a_-1)/2):ceil((a_-1)/2);
%         set_range_even = @(a_) -ceil((a_-1)/2):(ceil((a_-1)/2-1));
%         set_max = @(a_) max([ceil((a_(1)-1)/2),ceil((a_(2)-1)/2)]);
%         
%         ref = [];
%         if mod(Matrix_Size(1) , 2) == 0
%             if mod(Matrix_Size(2) , 2) == 0
%                 blup1 = set_range_even(Matrix_Size(2));
%                 blup2 = set_range_even(Matrix_Size(1));
%                 ref = [blup1(end),blup2(end)];
%                 [ky,kx] = meshgrid(blup1,blup2);
%             else
%                 blup1 = set_range_odd(Matrix_Size(2));
%                 blup2 = set_range_even(Matrix_Size(1));
%                 ref = [blup1(end),blup2(end)];
%                 [ky,kx] = meshgrid(blup1,blup2);
%             end
%         else
%             if mod(Matrix_Size(2) , 2) == 0
%                 
%                 blup1 = set_range_even(Matrix_Size(2));
%                 blup2 = set_range_odd(Matrix_Size(1));
%                 ref = [blup1(end),blup2(end)];
%                 [ky,kx] = meshgrid(blup1,blup2);
%             else
%                 blup1 = set_range_odd(Matrix_Size(2));
%                 blup2 = set_range_odd(Matrix_Size(1));
%                 ref = [blup1(end),blup2(end)];
%                 [ky,kx] = meshgrid(blup1,blup2);
%             end
%         end
%         
%         max_ = set_max(ref);
%         kspace = cat(3,kx/max_,ky/max_);
        
        USF = Undersampling_factor;
        rho_a  = Rho_a;
        rho_b = Rho_b;
        
        
        variable_density_undersampling_pattern = SPk_nD(kspace,rho_a, rho_b);
        %mu = mean(variable_density_undersampling_pattern(:));
        %sigma = std(variable_density_undersampling_pattern(:));
        mu = 0;
        sigma = 1;
        random_pdf = pdf_nD(kspace,mu,sigma);
        
        samples = ceil(numel(variable_density_undersampling_pattern)/USF);
        [M,N,P] = size(variable_density_undersampling_pattern);
        
        [usf_sorted, positions] = sort(VEC(variable_density_undersampling_pattern.*random_pdf),'descend');
        us_pattern = zeros(numel(variable_density_undersampling_pattern),1);
        us_pattern(positions(1:samples)) = 1;
        
        USP = squeeze(reshape(us_pattern,[M,N,P]));
        
        %%% DO EXTRA PROCEDURES
        % order:
        % parallel imaging
        % calibration lines
        % partial fourier
        % number of coils
        if ~isempty(PIx) % One of the two is enough
            % PARALLEL IMAGING - STANDARD
            % 
            PIxy_matrix = zeros(size(USP));
            PIxy_matrix((1+dPIx):PIx:end,(1+dPIy):PIy:end) = 1;
            
            USP = USP.*PIxy_matrix;
        end
        
        if ~isempty(CL)
            % CALIBRATION LINES
            % MAKE SURE THAT THE CALIBRATION LINES ARE ACQUIRED
            
            [X,Y] = size(USP);
            midX = ceil((X-1)/2);
            midY = ceil((Y-1)/2);
            
            if mod(X,2) == 0
                x_range = (midX -  CL/2) : ( midX + CL/2 -1) ;
            else
                x_range = (midX -  CL/2) : ( midX + CL/2) ;
            end
        
            if mod(Y,2) == 0
                y_range = (midY  -  CL/2) : ( midY + CL/2 -1) ;
            else
                y_range = (midY  -  CL/2) : ( midY + CL/2) ; % -ceil((Y-1)/2) Again, from left to right
            end
            
            USP(x_range, y_range) = 1;
           
            % REAL QUESTION: ONLY LEFT TO RIGHT?
        end
        
        if ~isempty(PF)
            % PARTIAL FOURIER
            % ERASE PHASE ENCODING LINES FROM LEFT TO RIGHT
            
            zz = round(PF*size(USP,2)); 
            USP(:,1:zz) = 0;
            
            % REAL QUESTION: ONLY LEFT TO RIGHT?
            % REAL QUESTION: ROUND OR CEIL???
        end
        
        combined_acceleration_factor = numel(USP)/sum(VEC(USP));
        disp(strcat('combined Acceleration Factor:',num2str(combined_acceleration_factor)));
        
        if ~isempty(COILS)
            dimss = size(USP);
            vec_dimss = zeros(numel(dimss)+1,1);
            vec_dimss(1:numel(dimss)) = 1;
            vec_dimss(end) = COILS;
            USP = repmat(USP , vec_dimss');
        end
        
        
    else
        %%% SOMETHING IS WRONG OR NOT IMPLEMENTED YET ... SKIP
        % not implemented yet
        USP = [];
        
    end

end