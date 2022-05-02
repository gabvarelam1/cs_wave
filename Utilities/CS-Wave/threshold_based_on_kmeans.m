function [threshold] = threshold_based_on_kmeans(waveletCoefficients)

    N = 512;
    data = (waveletCoefficients(:));
    
    
    %%% option 1:
    %%% directly from the wavelet coefficients
    if false
        map_of_groups = kmeans(data,2);
        [val,pos] = min(abs(data));

        noise_coefficients = data(map_of_groups == map_of_groups(pos));
        signal_coefficients = data(map_of_groups ~= map_of_groups(pos)); 

        threshold = max(noise_coefficients);
    end
    
    %%%% option 2:
    %%% based on histogram
    if true
        
        [count, bins] = genHistogram(data, N);
        [val , pos] = max([count]);

        map_of_groups = kmeans([count],2);
        bin_pos = find(map_of_groups == map_of_groups(pos),1,'last')+3;

        threshold = bins(bin_pos);% + .05*bins(bin_pos);
    end
end

function [Count, Bins] = genHistogram(Data, N)

    Maxd = max(Data);
    Mind = min(Data);
    Bins = linspace(Mind-.1*Mind, Maxd+.1*Maxd, N);
    
    Count = zeros(N,1);
    for iter = 2 : N
       Count(iter-1) = sum((Data >= Bins(iter-1)) .* (Data < Bins(iter)) ); 
    end
end