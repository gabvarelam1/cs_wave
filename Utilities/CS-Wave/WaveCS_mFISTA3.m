function [Yn, YnSets, OptKit] = WaveCS_mFISTA3(y,S,C,x,WaveCSparams,WaveletName,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% REFERENCES: 

%%% 1. "General phase regularized reconstruction using phase cycling"
%%% F. Ong, J. Cheng, M. Lustig, 2018.

%%% 2."Automatic determination of the regularization weighting for
%%% Wavelet-based compressed sensing reconstructions" 
%%% Authors: G. Varela-Mattatall, C. Baron, R. Menon, 2020.

%%% Note: Equations labeling are from reference 2.

%%% FUNCTION'S OBJECTIVE 

%%% This function solves, 

%%%     x_opt = min x     ||Ax-y||_2^2 + \lambda|| Wx ||_1,      eq.[2]

%%% using a proximal gradient descent algorithm of the style

%%%     x_n+1 = prox_operator_{\lambda||Wx||_1}{x_n - ....
%%%                             alpha \nabla{Ax_n - y}}          eq.[3]

%%% where the proximal operator is the soft-thresholding function,S{*},
%%%     prox_operator_{\lambda||Wx||_1} = W^{-1}S_\lambda{Wx}.   eq.[4]

%%% This was implemented as a "practical" monotonic FISTA algorithm

%%% Function parameters and examples

%%% how to call this function?
%%% [Yn, Convergence_rate, Level , Lambda] = complex_mFISTA2(Kspace,Fourier_op,Sensitivity_Coils,Chemical_shift,StartingPoint,WaveletName,'Lambda',Arg1,'Levels',Arg2)

%%%%%%%%%%%%%%%%
%%%   INPUT
%%%

%%%
%%% mandatory
%%%

%%% Kspace,
%%% Sampling Pattern,
%%% Sensitivity Coil Maps,
%%% StartingPoint,
%%% WaveCSparameters,
%%% WaveletName

%%%
%%% optional
%%% 

%%% 'Lambda',
%%% Arg1,
%%% 'Levels' 
%%% Arg2

%%%%%%%%%%%%%%%
%%% OUTPUT
%%%

%%% Yn, 
%%% Convergence_rate, 
%%% Level
%%% Lambda
%%%

%%%%%%%%%%%%%%%%%%%%%%
%%% EXAMPLES - OPTIONS
%%%

%%% [...] = complex_mFISTA2(y,F,S,C,x,'db4','Lambda','single','Levels','automatic')
%%% [...] = complex_mFISTA2(y,F,S,C,x,'db4','Lambda','single','Levels',4)
%%% [...] = complex_mFISTA2(y,F,S,C,x,'db4','Lambda','multiple','Levels','automatic')
%%% [...] = complex_mFISTA2(y,F,S,C,x,'db4','Lambda','multiple','Levels',4)


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% --- VERSION 2.0 ---  
%%% by GVM on May 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% This was implemented as a "practical" monotonic FISTA algorithm
    
    
    % 0. Setting optional parameters
    LAMBDA = [];
    LEVELS = [];
    FIXED_LAMBDA = [];

    if ~isempty(varargin)
        valid_argnames = {'Lambda','Levels','Fixed Lambda'};
        [args_specified,args_location] = ismember(valid_argnames,varargin(1:2:end));
        
        if args_specified(1) == true
            % user defined
            LAMBDA = varargin{(args_location(1)*2-1)+1};
            
            %verify that user didnt mess up with the argument
            if strcmp(LAMBDA,'multiple')
            elseif strcmp(LAMBDA,'single')
            else
                error('Lamba argument is not in the list! Opts: 1. multiple; 2. single');
            end
        else
            % default
            LAMBDA = 'multiple';
        end
        if args_specified(2) == true
            LEVELS = varargin{(args_location(2)*2-1)+1};
            
            %verify that user didnt mess up with the argument
            if strcmp(LEVELS,'automatic')
            elseif (LEVELS > 0) && (LEVELS < 7)
            else
                error('Levels argument is not in the list! Opts: 1. "automatic"; 2. {1,2,...,7}');
            end
        else
            % default
            LEVELS = 'automatic';
        end
        if args_specified(3) == true
            FIXED_LAMBDA = varargin{(args_location(3)*2-1)+1};
        else
            FIXED_LAMBDA = [];
        end
                
    end

    
    % 1. LEVELS 
    if strcmp(LEVELS,'automatic')
        % 1. Automatic decomposition Level that maximizes sparsity                Eq. [6]&[10]
        [decomposition_level_r] = nD_wavelet_level_selector(real(x), WaveletName);
        [decomposition_level_i] = nD_wavelet_level_selector(imag(x), WaveletName);
        %
        L = min([decomposition_level_r , decomposition_level_i]);
    
    else
        % Arbitrary decompostion level
        L = LEVELS;
        
    end

    % 2. LAMBDA
    dimensions = numel(size(x));
    
    %%% NOTE: For now dimensions are expected to be between 1-3 but only spatial
    % ones. Version 2.0 needs to take into account multiple echo-times or
    % other variables from the dataset.
    
    if strcmp( LAMBDA , 'multiple') && isempty(FIXED_LAMBDA)
        %%%% NEED TO REMOVE THE TRIMMING PARAMETER T
        % analysis on each high-frequency level based on k-means
        %tic
        list_lambdas = selection_adaptive_lambda_kmeans(abs(x) , WaveletName, L , dimensions, LAMBDA); % Eq. [7]
        list_lambdas(1,:) = 0; % lambdas are only for detail coefficients
        %time_lambdas = toc; 
        %disp(time_lambdas);
    elseif strcmp( LAMBDA , 'single') && isempty(FIXED_LAMBDA)
        % analysis using all high-frequency levels
        list_lambdas = selection_adaptive_lambda_exp2model(abs(x) , WaveletName, L , dimensions, LAMBDA); % Eq. [7]
        list_lambdas(1,:) = 0; % lambdas are only for detail coefficients
    else
        num_dims = numel(size(abs(x)));
     
        switch num_dims
            case 1
            case 2
            case 3
                list_lambdas = ones((7*L)+1,1)*FIXED_LAMBDA;
        end
    end
    
    
    % 3. Create Proximal Wavelet L1-norm operator                   Eq. [8]
    prox_reg = create_proximal_operator_L1(WaveletName, L, list_lambdas , dimensions );
    %%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% OPTIMIZATION
    
    %%% This was implemented as a "practical" monotonic FISTA algorithm
    
    niter = 50; % We are going to set it to 50 iterations just to check convergence
                % but we are going to use the results from iter 15.
    alpha = 1;
    yn = x;
    cost_opt_past = 1000000;
    checker = round(15:15:niter);
    convergence_rate = zeros(niter,1);
    
    % cost function optimization
    % F = \| Ax - y \|_2^2 + \lambda \| Wx \|_1
    % F = f + g
    
    mae = @(a_) norm(vec(a_));
    cost_f = @(a_) mae(a_)^2;  
    
    %%% WaveCS parameters and built-in functions
    %%% from k-space to image
    
    crop_size = WaveCSparams.crop_size;
    pad_size = WaveCSparams.pad_size;
    PSF = WaveCSparams.PSF;
    %Nfactor = WaveCSparams.Nfactor;
    %spsf = size(PSF,1);
    %sx = size(C,1);
    %factor_cycles = WaveCSparams.factor_cycles;
    %oversampling_factor = WaveCSparams.oversampling_factor;
    number_channels = WaveCSparams.number_channels;
    clear WaveCSparams;
    
    fft1Dc = @(Ksp, Dim) fftshift( fft( ifftshift(Ksp , Dim) ,[],Dim) ,Dim) / sqrt(size(Ksp,Dim));
    ifft1Dc = @(Im, Dim) fftshift( ifft( ifftshift(Im , Dim) ,[],Dim) ,Dim) * sqrt(size(Im, Dim));

    iF1d = @(ksp_,dim_) ifft1Dc(ksp_,dim_);
    F1d = @(ksp_,dim_) fft1Dc(ksp_,dim_);

    % Forward
    Cs = @(sensitivity_coils , im_) padarray(repmat(im_,[1,1,1,number_channels]) .* sensitivity_coils, [pad_size,0,0,0]);
    PFx = @(Cs_ , wave_point_spread_function) bsxfun(@times, F1d(Cs_,1) , wave_point_spread_function);
    Fyz = @(im_ , S_) S_ .* ( F1d( F1d(im_, 2) ,3 ) ) ;


    %  Inverse
    iFyz = @(ksp_, S_) iF1d( iF1d(ksp_ .* S_, 3) , 2 );
    iPFx = @(Im_yz_kspace_x , wave_point_spread_function) ...
        iF1d( bsxfun(@times,Im_yz_kspace_x , wave_point_spread_function) , 1 );
    Ct = @(sensitivity_coils , im_) sum( ...
    conj(sensitivity_coils) .* im_(1+end/2-crop_size/2:end/2+crop_size/2,:,:,:), 4);


    %%% monotomic FISTA :: mFISTA
    % if cost function starts to diverge, it restarts the updating step...
    YnSets = zeros([size(yn),length(checker)-1]);
    count_checker = 1;
    iter_time = zeros(niter,1);
    
    disp(strcat('Iteration:',num2str(1),'/',num2str(niter)));
    for it = 1:niter
        tic
        if sum(checker == it) ~= 0
            disp(strcat('Iteration:',num2str(it),'/',num2str(niter)));
        end
        %%% change here for WAVECS 
        tempyn = PFx( Cs(C, yn) , PSF);
        tempyn = Fyz( tempyn , S );
        tempyn = iFyz( y - tempyn , S );
        %r = Ct( C , iPFx( tempyn , PSFt ) );
        r = Ct( C , iPFx( tempyn , conj(PSF) ) );
        
        % r = C' * (S' * (F' * (y - (F * (S * (C * yn))))));
        
        cost_opt = cost_f(r) + cost_g(yn , list_lambdas, WaveletName, L, dimensions);

        %if it > 50 
            if cost_opt < cost_opt_past
                cost_opt_past = cost_opt;
            else
                alpha = 1;
            end
        %end

        xn = prox_reg( yn + alpha * ( r ), alpha);
        temp_alpha = alpha;
        alpha = (1 + sqrt(1 + 4*alpha^2)) / (2); 
        yn = xn + ((temp_alpha - 1)/alpha)*(xn - yn);

        convergence_rate(it) = cost_opt;
        
        if sum(checker == it) ~= 0
            YnSets(:,:,:,count_checker) = yn;
            count_checker = count_checker + 1;
        end
        iter_time(it) = toc;
    end
    
    %%% END OPTIMIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% 

    Yn = yn; 
    
    OptKit = struct;
    OptKit.Convergence_rate = convergence_rate;
    OptKit.Iter_time = iter_time;
    OptKit.Wavelet_Level = L;
    OptKit.Lambda = list_lambdas;
    OptKit.Final_F = cost_f(r);
    OptKit.Final_G = cost_g_no_lambda(yn, WaveletName, L, dimensions);
    
end

%%% END FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%

%%%% ADDITIONAL FUNCTIONS
% there are here to avoid data transfer...and speed up performance,
% Furthermore, we need the lambdas and wavelet transform to compute cost
% function of Regularizer.

% Enumeration:
% 1. create_proximal_operator_L1
%   i.      helper1    
%   ii.     define_factors
%   iii.    helper2
%   iv.     helper3
% 2. cost_g
%

function [Operator] = create_proximal_operator_L1(wname, N, List_of_Lambdas, signal_dims)

    switch signal_dims
        case 1
            Operator = @(x_ , alpha_) helper1(x_ , wname, N, List_of_Lambdas, alpha_);
        case 2
            Operator = @(x_ , alpha_) helper2(x_ , wname, N, List_of_Lambdas, alpha_);
        case 3
            Operator = @(x_ , alpha_) helper3(x_ , wname, N, List_of_Lambdas, alpha_);
        otherwise
            disp('Unknown number of dimensions to define a wavelet proximal operator');
            Operator = [];
    end

end

function [signal_rec] = helper1(signal, wname, N, List_of_Lambdas, alpha )
    
    set_range = @(Book_,Extraction_Number_)(1 + sum(Book_(1:(Extraction_Number_-1))) : sum(Book_(1:Extraction_Number_)));

    [signal, rshift] = randshift(signal);
    dwtmode('ppd', 0);
    
    s = size(signal);
    signal_rec = zeros(s);
    
    n_dims = 1;
    
    for i = 1:prod(s(2:end))

        % wavelet coefficients
        [coeffs_real, book_real] = wavedec(real(signal(:,i)),N,wname);
        [coeffs_imag, book_imag] = wavedec(imag(signal(:,i)),N,wname);
        
        for rj = 2:N+1
                                       
                rj_data_real = nD_extract_coefficients_from_wavelets(coeffs_real, book_real, n_dims , rj , wname, []);
                rj_data_imag = nD_extract_coefficients_from_wavelets(coeffs_imag, book_imag, n_dims , rj , wname, []);


                rj_data_complex = complexSoftThresh(rj_data_real + 1i*rj_data_imag, List_of_Lambdas(rj)*alpha );
                

                coeffs_real( set_range( book_real , rj ) ) = real(rj_data_complex);
                coeffs_imag( set_range( book_imag , rj ) ) = imag(rj_data_complex);
        end

            
        % combination
        real_part = waverec(coeffs_real, book_real, wname);
        imag_part = waverec(coeffs_imag, book_imag, wname);
        signal_rec(:,i) = real_part + 1i * imag_part;
    end

    signal_rec = randunshift( signal_rec, rshift );

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
        %disp('Unknown behaviour for L1 proximal operator');
        Fact_3 = 0;
        Fact_4 = 0;  
    end

end

function signal_rec = helper2(signal, wname, N, List_of_Lambdas, alpha)
  
    [signal, rshift] = randshift(signal);


    dwtmode('ppd', 0);

    s = size(signal);
    signal_rec = zeros(s);
    
    n_dims = 2;
    dtype = {'H','V','D'};

    for i = 1:prod(s(3:end))

        % real part
        [coeffs_real, book_real] = wavedec2(real(signal(:,:,i)),N,wname);    
        [coeffs_imag, book_imag] = wavedec2(imag(signal(:,:,i)),N,wname);
        
        for rj = 2:N+1
            for rj2 = 1:3
                [approx_factor , fact_1 , fact_2, fact_3, fact_4, init_index] = define_factors(rj, dtype{rj2});

                set_ranges_2D = @(Book , Extraction_Number) ...
                    (1 + approx_factor*sum(prod(Book(1,:),2)) + fact_1*sum(prod(Book(init_index:(Extraction_Number-1),:),2)) + ...
                                           fact_3*sum(prod(Book(Extraction_Number,:),2))) : ...
                    ( approx_factor*sum(prod(Book(1,:),2)) + fact_2*sum(prod(Book(init_index:Extraction_Number-1,:),2)) + ...
                                           fact_4*sum(prod(Book(Extraction_Number,:),2)));
                
                %disp(length(set_ranges_2D(book_real,rj)));
                
                rj_data_real = nD_extract_coefficients_from_wavelets(coeffs_real, book_real, n_dims , rj, 'both', dtype{rj2});
                rj_data_imag = nD_extract_coefficients_from_wavelets(coeffs_imag, book_imag, n_dims , rj, 'both', dtype{rj2});


                rj_data_complex = complexSoftThresh(rj_data_real + 1i*rj_data_imag, List_of_Lambdas(rj) * alpha);
                

                coeffs_real( set_ranges_2D( book_real , rj ) ) = real(rj_data_complex);
                coeffs_imag( set_ranges_2D( book_imag , rj ) ) = imag(rj_data_complex);
                
            
            end
        end

        real_part = waverec2(coeffs_real, book_real, wname);
        imag_part = waverec2(coeffs_imag, book_imag, wname);
        signal_rec(:,:,i) = real_part + 1i * imag_part;

    end

    signal_rec = randunshift( signal_rec, rshift );


end

function signal_rec = helper3(signal, wname, N, List_of_Lambdas, alpha)
%%% built-in functions for the job!

    [signal, rshift] = randshift(signal);


    dwtmode('ppd', 0);

    s = size(signal);

    signal_rec = zeros(s);

    for i = 1:prod(s(4:end))

        
        real_stwave = wavedec3(real(signal(:,:,:,i)),N,wname);
        imag_stwave = wavedec3(imag(signal(:,:,:,i)),N,wname);

        Nw = length(real_stwave.dec);
        for rj = 2:Nw
            
            size_wav = size(real_stwave.dec{rj});
            
            rj_data_real = vec(real_stwave.dec{rj});
            rj_data_imag = vec(imag_stwave.dec{rj});

            rj_data_complex = complexSoftThresh(rj_data_real + 1i*rj_data_imag, List_of_Lambdas(rj) * alpha);
            
            real_stwave.dec{rj} = reshape( real(rj_data_complex) , size_wav );
            imag_stwave.dec{rj} = reshape( imag(rj_data_complex) , size_wav );
            
        end
        % combination
        signal_rec(:,:,:,i) = waverec3( real_stwave ) + 1i * waverec3( imag_stwave );        
    end

    signal_rec = randunshift( signal_rec, rshift );

end


%%%%%% cost function regularizer \lambda \| Wx \|_1
function [cost_function_regularizer] = cost_g(yn , List_of_Lambdas, wname, N, dimensions)
    
    cost_function_regularizer = 0;
    
    switch dimensions
        case 1
            s = size(yn);
            for i = 1 : prod( s(2:end) )    
                [coeffs_real, book_real] = wavedec(real(yn(:,i)),N,wname);
                [coeffs_imag, book_imag] = wavedec(imag(yn(:,i)),N,wname);
                
                for rj = 2:N+1

                        rj_data_real = nD_extract_coefficients_from_wavelets(coeffs_real, book_real, 1 , rj, 'both', 'all');
                        rj_data_imag = nD_extract_coefficients_from_wavelets(coeffs_imag, book_imag, 1 , rj, 'both', 'all');

                        
                        cost_function_regularizer = cost_function_regularizer + List_of_Lambdas(rj)*norm(vec(rj_data_real + 1i*rj_data_imag),1);

                end
                
            end
        case 2
            s = size(yn);
            for i = 1 : prod( s(3:end) )    
                [coeffs_real, book_real] = wavedec2(real(yn(:,:,i)),N,wname);
                [coeffs_imag, book_imag] = wavedec2(imag(yn(:,:,i)),N,wname);
                
                for rj = 2:N+1

                    rj_data_real = nD_extract_coefficients_from_wavelets(coeffs_real, book_real, 2 , rj, 'both', 'all');
                    rj_data_imag = nD_extract_coefficients_from_wavelets(coeffs_imag, book_imag, 2 , rj, 'both', 'all');


                    cost_function_regularizer = cost_function_regularizer + List_of_Lambdas(rj)*norm(vec(rj_data_real + 1i*rj_data_imag),1);

                end
                
            end
            
        case 3
            s = size(yn);
            for i = 1 : prod(s(4:end))
                real_stwave = wavedec3(real(yn(:,:,:,i)),N,wname);
                imag_stwave = wavedec3(imag(yn(:,:,:,i)),N,wname);

                Nw = length(real_stwave.dec);
                for rj = 2:Nw

                    rj_data_real = vec(real_stwave.dec{rj});
                    rj_data_imag = vec(imag_stwave.dec{rj});

                    cost_function_regularizer = cost_function_regularizer + List_of_Lambdas(rj)*norm(vec(rj_data_real + 1i*rj_data_imag),1);
                end
            end
    end
    
end

%%%%%%

%%%%%% cost function regularizer no lambda \| Wx \|_1
function [cost_function_regularizer] = cost_g_no_lambda(yn , wname, N, dimensions)
    
    cost_function_regularizer = 0;
    
    switch dimensions
        case 1
            s = size(yn);
            for i = 1 : prod( s(2:end) )    
                [coeffs_real, book_real] = wavedec(real(yn(:,i)),N,wname);
                [coeffs_imag, book_imag] = wavedec(imag(yn(:,i)),N,wname);
                
                for rj = 2:N+1

                        rj_data_real = nD_extract_coefficients_from_wavelets(coeffs_real, book_real, 1 , rj, 'both', 'all');
                        rj_data_imag = nD_extract_coefficients_from_wavelets(coeffs_imag, book_imag, 1 , rj, 'both', 'all');

                        
                        cost_function_regularizer = cost_function_regularizer + norm(vec(rj_data_real + 1i*rj_data_imag),1);

                end
                
            end
        case 2
            s = size(yn);
            for i = 1 : prod( s(3:end) )    
                [coeffs_real, book_real] = wavedec2(real(yn(:,:,i)),N,wname);
                [coeffs_imag, book_imag] = wavedec2(imag(yn(:,:,i)),N,wname);
                
                for rj = 2:N+1

                    rj_data_real = nD_extract_coefficients_from_wavelets(coeffs_real, book_real, 2 , rj, 'both', 'all');
                    rj_data_imag = nD_extract_coefficients_from_wavelets(coeffs_imag, book_imag, 2 , rj, 'both', 'all');


                    cost_function_regularizer = cost_function_regularizer + norm(vec(rj_data_real + 1i*rj_data_imag),1);

                end
                
            end
            
        case 3
            s = size(yn);
            for i = 1 : prod(s(4:end))
                real_stwave = wavedec3(real(yn(:,:,:,i)),N,wname);
                imag_stwave = wavedec3(imag(yn(:,:,:,i)),N,wname);

                Nw = length(real_stwave.dec);
                for rj = 2:Nw

                    rj_data_real = vec(real_stwave.dec{rj});
                    rj_data_imag = vec(imag_stwave.dec{rj});

                    cost_function_regularizer = cost_function_regularizer + norm(vec(rj_data_real + 1i*rj_data_imag),1);
                end
            end
    end
    
end
