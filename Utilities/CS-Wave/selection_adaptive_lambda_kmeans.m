function [list_lambdas] = selection_adaptive_lambda_kmeans(x0 , wname, N, signal_dims, flag_levels)
%%% this functions determines for the initial zero-filled image the
%%% approximate values for the threshold in each decomposition level
    

    switch signal_dims
        case 1
            list_lambdas =  helper1(x0 , wname, N, flag_levels);
        case 2
            list_lambdas =  helper2(x0 , wname, N, flag_levels);
        case 3
            list_lambdas =  helper3(x0 , wname, N, flag_levels);
        otherwise
            disp('Unknown number of dimensions to define a wavelet proximal operator');
            list_lambdas = [];
    end

end

function List_of_Lambdas = helper1(signal, wname, N, flag_levels )
    
    
    set_range = @(Book_,Extraction_Number_)(1 + sum(Book_(1:(Extraction_Number_-1))) : sum(Book_(1:Extraction_Number_)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [signal, r] = randshift(signal);
    % dwtmode('ppd', 0);

    s = size(signal);
    
    n_dims = 1;
    
    List_of_Lambdas = zeros(1+N,prod(s(2:end))); % app + detail levels
    
    for i = 1:prod(s(2:end))

        if strcmp(flag_levels,'multiple')
             % extract high-frequency wavelet coefficients level-wise for kmeans
            List_of_Lambdas = zeros(1+N,prod(s(2:end)));

            [coeffs_real, book_real] = wavedec(real(signal(:,i)),N,wname);

            for rj = 2:N+1

                rj_data = ...
                    nD_extract_coefficients_from_wavelets(coeffs_real, book_real, n_dims , rj , wname, []);

                th_Lambda = threshold_based_on_kmeans_v2(rj_data);

                List_of_Lambdas(rj,i) = th_Lambda;
                %coeffs_real( set_range( book_real , rj )) = SoftThresh( rj_data, th_Lambda );
            end
        else
            % extract all high-frequency wavelet coefficients for kmeans
            rj_data = [];
            for rj = 2:N+1

                rj_data = [ rj_data , 
                    nD_extract_coefficients_from_wavelets(coeffs_real, book_real, n_dims , rj , wname, []) ];
            end
            
            th_Lambda = threshold_based_on_kmeans_v2(rj_data);
            List_of_Lambdas(2:end,i) = th_Lambda;
              
        end

    end

end

function [Approx_factor, Fact_1 , Fact_2, Fact_3, Fact_4, Init_index] = define_factors(Extraction_Number, Detail_Type)

    Init_index = 2; 
    switch Extraction_Number
        case 1
            % just approx coefficients
            Fact_1 = 1;
            Fact_2 = 1;
            Init_index = 1;
            Approx_factor = 0;
        case 2
            % boundary between approx and detail coefficients
            Fact_1 = 0;
            Fact_2 = 3;
            Approx_factor = 1;
        otherwise
            % detail coefficients
            Fact_1 = 3;
            Fact_2 = 3;
            Approx_factor = 1;
    end
    
    if strcmp(Detail_Type,'H')
        Fact_3 = 0;
        Fact_4 = 1;
    elseif strcmp(Detail_Type,'V')
        Fact_3 = 1;
        Fact_4 = 2;
    elseif strcmp(Detail_Type,'D')
        Fact_3 = 2;
        Fact_4 = 3;        
    else
        disp('Unknown behaviour for L1 proximal operator');
        Fact_3 = 0;
        Fact_4 = 0;  
    end

end

function List_of_Lambdas = helper2(signal, wname, N, flag_levels)
    
    s = size(signal);
    
    n_dims = 2;
    dtype = {'H','V','D'};
    
    List_of_Lambdas = zeros(1+N,prod(s(3:end)));

    for i = 1:prod(s(3:end))
        
            [coeffs_real, book_real] = wavedec2(real(signal(:,:,i)),N,wname);

            rj_data = ...
                    nD_extract_coefficients_from_wavelets(coeffs_real, book_real, n_dims , 1 , 'both', []);

            th_Lambda = threshold_based_on_kmeans_v2(rj_data);

            List_of_Lambdas(1,i) = th_Lambda; 
        
        if strcmp(flag_levels,'multiple')
            
            for rj = 2:N+1
                rj_per_level = [];
                for rj2 = 1:3

                    rj_data = nD_extract_coefficients_from_wavelets(coeffs_real, book_real, n_dims , rj, 'both', dtype{rj2});

                    rj_per_level = [rj_per_level, rj_data];

                end
                th_Lambda = threshold_based_on_kmeans_v2(rj_per_level);

                List_of_Lambdas(rj,i) = th_Lambda;
            end

        else
        % extract all high-frequency wavelet coefficients for kmeans

            [coeffs_real, book_real] = wavedec2(real(signal(:,:,i)),N,wname);

            rj_data = ...
                    nD_extract_coefficients_from_wavelets(coeffs_real, book_real, n_dims , 1 , 'both', []);

            th_Lambda = threshold_based_on_kmeans_v2(rj_data);

            List_of_Lambdas(1,i) = th_Lambda; 

            rj_data_total = [];
            for rj = 2:N+1
                for rj2 = 1:3

                    rj_data = nD_extract_coefficients_from_wavelets(coeffs_real, book_real, n_dims , rj, 'both', dtype{rj2});

                    rj_data_total = [rj_data_total, rj_data];

                end

            end
            th_Lambda = threshold_based_on_kmeans_v2(rj_data_total);
            List_of_Lambdas(2:end,i) = th_Lambda;
        end
    end

end

function  List_of_Lambdas = helper3(signal, wname, N, flag_levels)
    
    s = size(signal);

    List_of_Lambdas = zeros( N*7 + 1 , 1 );

    for i = 1:prod(s(4:end))

        if strcmp(flag_levels,'multiple')
        
            real_stwave = wavedec3(real(signal(:,:,:,i)),N,wname);

            %Nreal = length(real_stwave.dec);
            for rj = 2:(N+1)
                
                    rj_data_total = [];
                    for jk = 1 : 7
                        rj_data = vec(real_stwave.dec{1+7*(rj-2)+jk});
                        rj_data_total = [rj_data_total,rj_data];
                    end

                    th_Lambda = threshold_based_on_kmeans_v2(rj_data_total);

                    List_of_Lambdas(1+7*(rj-2)+1 : 1+7*(rj-2)+7,i) = th_Lambda; 

            end
            
        else
            
            real_stwave = wavedec3(real(signal(:,:,:,i)),N,wname);

            Nreal = length(real_stwave.dec);
            rj_data = [];
            for rj = 2:Nreal

                rj_data = [ rj_data , vec(real_stwave.dec{rj}) ];

            end
            
            th_Lambda = threshold_based_on_kmeans_v2(rj_data);
            List_of_Lambdas(2:end,i) = th_Lambda; 
            
        end

    end

end