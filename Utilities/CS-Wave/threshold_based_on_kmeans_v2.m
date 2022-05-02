function threshold = threshold_based_on_kmeans_v2(waveletCoefficients)
    % Automatic lambda selection, similar to Varela-Mattatall G, Baron CA,
    % Menon RS. Magn. Reson. Med. 2021;86:1403–1419. Here we'll base the lambda
    % selection per level on histograms of wavelet xform of "zero filled"
    % image, which is more generally just the adjoint of the encoding matrix A.
    % We just use the location where the histogram of the magnitude image
    % is maximized.
    
    h = histogram(abs(waveletCoefficients(:)));
    vals = h.Values;
    vals(1) = []; % we remove the count from 0-valued coefficients 
                  % to make easier the splitting b
    [~,pos] = max(vals);
    kgroups = 2;
    map_of_groups = kmeans(vals',kgroups);
    bin_pos = find(map_of_groups == map_of_groups(pos),1,'last');
    
    % We divide by sqrt(2) b/c we're dealing with magnitude
    % while the bfista algorithm is using complex values for softthresholding.
    
    accelerator = 1;
    threshold = accelerator*(h.BinEdges(bin_pos) / sqrt(2)); % assuming a conservative value for the threshold
                                               % think of this as floor() instead of the average
                                               % between edges.
    close all;
end


