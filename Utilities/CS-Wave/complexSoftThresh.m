function [res] = complexSoftThresh(ComplexData, Lambda)

    res =  ( abs(ComplexData) > Lambda) .* (exp( 1i*angle(ComplexData) ) .* ( abs(ComplexData) - Lambda))  ;

end