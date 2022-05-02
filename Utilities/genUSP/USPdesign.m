function [UndersamplingPattern,UndersamplingPattern2D, FinalUSF , Dir2Save] = USPdesign(Method,Ry,Rz,Ncalib, FlagElliptical , Size_ksp)
%%% This function generates 2D undersampling patterns for ND data (N > 2)
%%% It assumes that the firt dimension, or the frequency is fully sampled
%%% and the undersampling is performed in the sencond and third dimensions.
%%% This function replicates the 3D undersampling to posterior dimensions.

%%% This function has a set of flavours but all of them start from
%%%     1. 'grappa' || 'cs'
%%% and then (independent options):
%%%     a. Include a calibration region (a.ka. Ncalib)
%%%     b. Remove corners from the pattern (a.k.a FlagElliptical)
%%%     c. more in the future...

%%%     Ultimately, the pattern is defined as USP = 1 -> a. || b. || c. ...

%%% Finally, Dir2Save is a suggestion for the directory to store results based on
%%% the generated undersampling, and this function provides a sense of the
%%% final acceleration factor.

% by gvm, 2021

    % Method of choice
    usfString = [];
    if strcmp(Method,'gr')
        corePattern = includeUniformSampling(Ry,Rz,Size_ksp);
        usfString = strcat('_gr_R',num2str(Ry),'x',num2str(Rz));
    elseif strcmp(Method,'cs')

        usf = Ry*Rz; % usf + Clreg makes the real USF

        randnum = @(mu_ , sigma_) mu_ + sigma_ .* (rand() - 0.5);
        %mask = wangs_VD_undersampling_pattern([size_ksp(2:3)],usf,randnum(3,.5),randnum(2.5,.6),'Calibration region',Ncalib);
        corePattern = wangs_VD_undersampling_pattern([Size_ksp(2:3)],usf,randnum(5,.5),randnum(1,.6));
        usfString = strcat('_cs_R',num2str(usf));
    else
        error('unknown undermsapling method. Options: {gr , cs}');
    end

    % Include calibration region 
    clString = [];
    if Ncalib > 0
       clRegion = includeCalibrationRegion(Ncalib, Size_ksp);
       mask = or(corePattern,clRegion);
       clString = strcat('_',num2str(Ncalib),'x',num2str(Ncalib),'CLreg');
    else
       mask = corePattern;
    end

    % last step: remove corners
    ellipString = []; 
    if FlagElliptical
        UndersamplingPattern = removeCorners(mask);
        ellipString = strcat('_Elliptical');
    else
        UndersamplingPattern = mask;
    end

    % final undersampling pattern
    FinalUSF = numel(UndersamplingPattern) / sum(UndersamplingPattern(:));
    disp(strcat('Total undersampling is:',num2str(FinalUSF),'x'));
    RtotalString = strcat('Rt',num2str(round(FinalUSF)));

    Dir2Save = strcat(RtotalString,usfString,clString,ellipString);
    % dir2save

    % finally, replicate
    UndersamplingPattern2D = squeeze(UndersamplingPattern);
    UndersamplingPattern = repmat(permute(UndersamplingPattern,[3,1,2]),[Size_ksp(1),1,1,Size_ksp(4:end)]);

end

function [mask] = includeUniformSampling(Ry,Rz,size_ksp)

    mask = zeros(size_ksp(2:3));
    mask(1:Ry:end,1:Rz:end) = 1;

    %S = repmat(permute(mask,[3,1,2]),[size_ksp(1),1,1,size_ksp(4:end)]);
end

function [mask]=includeCalibrationRegion(Ncalib, size_ksp)
    
    X = size_ksp(2);
    Y = size_ksp(3);

    mask = zeros(X,Y);
    CL = Ncalib;
    if mod(CL,2) == 0
        toaddCl = CL/2;
    else
        toaddCl = (CL-1)/2;
    end
    
    if mod(X,2) == 0
        midX = ceil((X-1)/2);
    else
        midX = ceil((X+1)/2);
    end
    
    if mod(Y,2) == 0
        midY = ceil((Y-1)/2);
    else
        midY = ceil((Y+1)/2);
    end
    
    if mod(CL,2) == 0
        x_range = (midX -  toaddCl) : ( midX + toaddCl -1) ;
        y_range = (midY  -  toaddCl) : ( midY + toaddCl -1) ;
    else
        x_range = (midX -  toaddCl) : ( midX + toaddCl) ;
        y_range = (midY  -  toaddCl) : ( midY + toaddCl) ;
    end

    mask(x_range, y_range) = 1;
end

function Mask = removeCorners(Mask)

%%% to remove corners we use the formula of a ellipse:
%%%          (x-h)^2            (y-k)^2
%%%         _________    +     _________      = 1
%%%      radius_in x ^2       radius_in y ^2

    [sx,sy] = size(Mask);
    
    % radius
    %R = min(sx,sy);
    %R = R/2;
    radx = sx/2;
    rady = sy/2;
    
    % centre
    if mod(sx,2) == 0
        midX = ceil((sx-1)/2);
    else
        midX = ceil((sx+1)/2);
    end
    
    if mod(sx,2) == 0
        midY = ceil((sy-1)/2);
    else
        midY = ceil((sy+1)/2);
    end
    
    % find positions
    [px,py] = find(Mask==1);
    
    % remove corners
    % Rpos = sqrt((px-midX).^2 + (py-midY).^2);
    Rpos = ( (px-midX).^2 /radx^2) + ( (py-midY).^2/rady^2 );
    for i = 1 : length(Rpos)
        %if Rpos(i) > R 
        if Rpos(i) > 1
            Mask(px(i),py(i)) = 0;
        end
    end
    
end