function [Selected_Level, Sparsity_per_level, Coefficients_per_level] = nD_wavelet_level_selector(Signal, Wavelet_Name , varargin)
%%% this function computes the ideal level of decomposition which maximizes
%%% the asymptotic sparsity from wavelets for n-th dimensional data.

%%% Example [level, Sparsity_per_level] = nD_wavelet_level_selector(I, wname);
%%% Input:
%%% (1) I: Data {1,2,3}D please ahaha
%%% (2) wname: mother wavelet name. Example 'db4'.
%%% (3) [Optinal] Mask: generates and idea if any part of the signal can be
%%%                     ignored. Otherwise, the algorithm assumes the
%%%                     complete Data as relevant for the algorithm.
%%% Output:
%%% (1) level: ideal level for wavelet decomposition
%%% (2) Sparsity_per_level: sparsity per level according to the Gini Index
%%% (3) [Optional] Coefficients_per_1D: a cell which size is according to
%%%                                     number dimensions from data and contains cells which contain
%%%                                     the coefficient per level of the different dimensions.  

%%%% by gab, 2019

data_dim = define_number_dimensions_from_signal(Signal);

%%% 2. Preset for Algorithm
if ~isempty(varargin)
    data = varargin{1} .* Signal;
else
    data = Signal;
end

wavelet_name = Wavelet_Name;

switch data_dim
    case 1
        [Selected_Level, Sparsity_per_level, Coefficients_per_level] = run_algorithm_1D(data,wavelet_name);
    case 2
        [Selected_Level, Sparsity_per_level, Coefficients_per_level] = run_algorithm_2D(data,wavelet_name);
    case 3
        [Selected_Level, Sparsity_per_level, Coefficients_per_level] = run_algorithm_3D(data,wavelet_name);
end



end

function [number_dimensions] = define_number_dimensions_from_signal(Signal)
    %%% Determine number of dimensions and check
    [a,b,c,d] = size(Signal);

    %%% 1. Determine the number of dimensions in the data
    %%% not the best procedure, but it helps if people uses a reasonably
    %%% structure for the input data...
    if ( a == 1 ) || ( b==1 )
        %%% one dimensional data
        number_dimensions = 1;
    elseif (c == 1) && (d == 1)
        %%% two dimensional data
        number_dimensions = 2;
    elseif (d == 1) 
        %%% three dimensional data
        number_dimensions = 3;
    elseif d > 1
        disp('Error ... impossible to characterize the data in nD_wavelet_level_selector function...');
        number_dimensions = [];
        return;
    end
end

function [appropriate_level, sparsity_per_level, coefficients_per_level] = run_algorithm_1D(data,wavelet_name)
    %%% Algorithm
    appropriate_level = 1;
    n_dims = 1;
    flag = 1;
    good_sparsity_level = 0.75;
    %
    while flag
        next_level = appropriate_level + 1;
        [coeffs,book] = wavedec(data,next_level,wavelet_name);
         
        % book contains: {AN,DN, DN-1, DN-2, ... , DN-N , Length of Signal}

        gi = zeros(next_level,1);
        coeffs_per_level = cell(size(gi));

        gi_low = Gini_Index(nD_extract_coefficients_from_wavelets(coeffs, book, n_dims , 1, wavelet_name, [] ));
        %coeffs_per_level{1} = coeffs( (1) : sum(book(1)) ); 
        for iter = (next_level+1):-1:2
            level_coeffs = nD_extract_coefficients_from_wavelets(coeffs, book, n_dims , iter, wavelet_name, [] );
            gi(iter-1) = Gini_Index(level_coeffs);
            coeffs_per_level{iter-1} = level_coeffs;
        end

        gi = flipud(gi);

        for iter = 1:(length(gi)-1)
            % condition to stop
            % Sj <= Tr && Sj+1 > Tr
            if  (gi(iter+1) - gi(iter)) < 0 && gi(iter) > good_sparsity_level
                flag = 0;
                break;
            end
        end

        if flag ~= 0
             % increase j
            appropriate_level = appropriate_level + 1;
        end

    end
    %appropriate_level = appropriate_level;
    sparsity_per_level = [gi_low ; flipud(gi(1:end-1))];
    coefficients_per_level = cell(appropriate_level+1,1);
    coefficients_per_level{1} = nD_extract_coefficients_from_wavelets(coeffs, book, n_dims , 1, wavelet_name, [] );
    for iter = 2:(appropriate_level+1)
        coefficients_per_level{iter} = coeffs_per_level{iter};
    end
end

function [appropriate_level, sparsity_per_level, coefficients_per_level] = run_algorithm_2D(data,wavelet_name)

    appropriate_level = 1;
    flag = 1;
    dtype = {'H','V','D'};
    n_dims = 2;
    good_sparsity_level = 0.6;
    %

    while flag
    
        next_level = appropriate_level + 1;
        
        [coeffs,book] = wavedec2(data,next_level,wavelet_name);
         
        % book contains: {AN,DN, DN-1, DN-2, ... , DN-N , Length of Signal}
        %Sj = zeros(appropriate_level+1,1);
        gi_mean = zeros(appropriate_level+1,1);
        coeffs_per_level = cell(size(gi_mean));
        
        gi_x = zeros(next_level,1);
        gi_y = zeros(next_level,1);
        gi_xy = zeros(next_level,1);
    %   Sparsity Values

    %    Sj(1) = [sj_funct( nD_extract_coefficients_from_wavelets(coeffs, book, n_dims , 1, 'both',[]) )];
        coeffs_low = nD_extract_coefficients_from_wavelets(coeffs, book, n_dims , 1, 'both',[]);
        gi_low = Gini_Index(coeffs_low); 
        for iter = (next_level+1):-1:2
             for iter2 = 1:3
                switch iter2
                    case 1
                        level_coeffs_x = nD_extract_coefficients_from_wavelets(coeffs, book, n_dims , iter , 'both', dtype{iter2});
                        gi_x(iter-1) = Gini_Index(level_coeffs_x);
                    case 2
                        level_coeffs_y = nD_extract_coefficients_from_wavelets(coeffs, book, n_dims , iter , 'both', dtype{iter2});
                        gi_y(iter-1) = Gini_Index(level_coeffs_y);
                   otherwise
                        level_coeffs_d = nD_extract_coefficients_from_wavelets(coeffs, book, n_dims , iter , 'both', dtype{iter2});
                        gi_xy(iter-1) = Gini_Index(level_coeffs_d);
                 end
             end
             %gi_mean(iter-1) = mean([gi_x(iter-1),gi_y(iter-1),gi_xy(iter-1)]);
             coeffs_per_level{iter-1} = [level_coeffs_x,level_coeffs_y,level_coeffs_d];
             gi_mean(iter-1) = Gini_Index(coeffs_per_level{iter-1});
         end

        gi_mean = flipud(gi_mean);

        for iter = 1:(size(gi_mean)-1)
            % condition to stop
            % Sj <= Tr && Sj+1 > Tr
            % if gi(iter) <= Tr && gi(iter+1) > Tr
            if (gi_mean(iter+1) - gi_mean(iter)) < 0 && gi_mean(iter) > good_sparsity_level
                flag = 0;
                break;
            end
        end

        if flag ~= 0
             % increase j
            appropriate_level = appropriate_level + 1;
        end
        
%         if appropriate_level == 10
%             disp('appropriate_level');
%         end

    end
    [coeffs,book] = wavedec2(data,appropriate_level,wavelet_name);
    coeffs_low = nD_extract_coefficients_from_wavelets(coeffs, book, n_dims , 1, 'both',[]);
    gi_low = Gini_Index(coeffs_low); 
    gi_mean = flipud(gi_mean(1:end-1));
    appropriate_level = appropriate_level;
    sparsity_per_level = [gi_low ; gi_mean];
    coefficients_per_level = cell(next_level,1);
    coefficients_per_level{1} = coeffs_low;
    for iter = 2:(next_level)
        coefficients_per_level{iter} = coeffs_per_level{iter};
    end
end

function [appropriate_level, sparsity_per_level, coefficients_per_level] = run_algorithm_3D(data,wavelet_name)

    appropriate_level = 1;
    flag = 1;
    good_sparsity_level = 0.45;
    %

    while flag
    
        [wstruct] = wavedec3(data,appropriate_level+1,wavelet_name);
        
        %Nl = appropriate_level + 2; 
        % book contains: {AN,DN, DN-1, DN-2, ... , DN-N , Length of Signal}
        %Sj = zeros(appropriate_level+1,1);
        gi_mean = zeros(appropriate_level+1,1);
        directions_per_level = 7;
        gi_3D_per_level = zeros(directions_per_level,1);
    %   Sparsity Values

    %    Sj(1) = [sj_funct( nD_extract_coefficients_from_wavelets(coeffs, book, n_dims , 1, 'both',[]) )];
        gi_mean(1) = Gini_Index( vec( wstruct.dec{1} ) ); 
        for iter = 2:appropriate_level+1
             for iter2 = 1:7
                gi_3D_per_level(iter2) = Gini_Index( vec( wstruct.dec{ ...
                                                         directions_per_level*(iter-2) + iter2 } ) );
             end
             gi_mean(iter) = mean(gi_3D_per_level);
         end

        gi_mean = flipud(gi_mean(2:end));

        for iter = 1:(size(gi_mean)-1)
            % condition to stop
            % Sj <= Tr && Sj+1 > Tr
            % if gi(iter) <= Tr && gi(iter+1) > Tr
            if (gi_mean(iter+1) - gi_mean(iter)) < 0 && gi_mean(iter) > good_sparsity_level
                flag = 0;
                break;
            end
        end

        if flag ~= 0
             % increase j
            appropriate_level = appropriate_level + 1;
        end
        
        if appropriate_level > 8
            warning('enforcing appropriate_level');
            appropriate_level = 3;
            flag = 0;
            break;
            
        end

    end
    appropriate_level = appropriate_level-1;
    sparsity_per_level = gi_mean;
    coefficients_per_level = [];
end