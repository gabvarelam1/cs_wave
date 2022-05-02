function GI = Gini_Index(SetCoefficients)

    s = sort(abs(SetCoefficients(:)),'ascend');
    N = length(s);
    s = reshape(s,[1,N]);
    v = 1:N;
    l1norm = norm(s,1);
   
    
    GI = 1 - 2*sum( (s / l1norm ) .* ( (N-v+0.5)/N )  );

return;