function Res_Qmetric = Qmetric(Multicoil_complex_data,Dim_for_coils) 

%%% Adaptation from gab 2022
%%% Eq [15] from review
%%% ""
%%% Q = 100 x | Sum_j M_j e^{i(phi_j corrected)} | / Sum_j Mj

%%% "When the phases of the individual signal vectors are in good agreement (i.e., they are matched),
%%% the length of the complex sum of the signals (numerator) is equal to the sum of the length
%%% of the individual signal vectors (denominator)...." (page 8)

%%% "...Q should only be evaluated in voxel containing signal..." (page 8)

mj = abs(Multicoil_complex_data);
aj = angle(Multicoil_complex_data);

Res_Qmetric = 100 * abs( sum( mj .* exp(1i*aj) , Dim_for_coils ) ) ./ ...
                        sum(mj , Dim_for_coils);

end