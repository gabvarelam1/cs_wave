function [coefficients] = nD_extract_coefficients_from_wavelets(Coefficients_Set, Book, Number_Dimensions_Data , Extraction_Number , Coefficient_Types, Detail_Type)

    %%% function to extract coefficients from a matlab nD wavelet operator.
    % it changes according to the number of dimensions because the operator
    % is messy... :/
    
    %%% gab, 2019
    
    
%%%
% set_ranges_approx = @(book_ , iter_) (1 + sum(prod(book_(1:(iter_-1),:),2)) ) : sum(prod(book_(1:iter_,:),2));
% set_ranges_details_all = @(book_ , iter_) (1 + 3*sum(prod(book_(1:(iter_-1),:),2)) ) : 3*sum(prod(book_(1:iter_,:),2));

    
    switch Number_Dimensions_Data
        case 1
            coefficients = Coefficients_Set( (1 + sum(Book(1:(Extraction_Number-1))) : sum(Book(1:Extraction_Number))));
        case 2
            if strcmp(Coefficient_Types,'both')
                init_index = 2;
                switch Extraction_Number
                    case 1
                        % just approx coefficients
                        fact_1 = 1;
                        fact_2 = 1;
                        init_index = 1;
                        approx_factor = 0;
                    case 2
                        % boundary between approx and detail coefficients
                        fact_1 = 0;
                        fact_2 = 3;
                        approx_factor = 1;
                    otherwise
                        % detail coefficients
                        fact_1 = 3;
                        fact_2 = 3;
                        approx_factor = 1;
                end
            elseif strcmp(Coefficient_Types,'approx')
                fact_1 = 1;
                fact_2 = 1;
                approx_factor = 0;
                init_index = 1;
            elseif strcmp(Coefficient_Types,'detail')
                fact_1 = 3;
                fact_2 = 3;
                approx_factor = 0;
                init_index = 1;
            else
                fact_1 = [];
                fact_2 = [];
                approx_factor = [];
                init_index = [];
                coefficients = [];
                disp('Error in function nD_extract_coefficients_from_wavelets');
                return
            end
            
            %%% approx or all-details
            if isempty(Detail_Type) || strcmp(Detail_Type,'all')
                coefficients = ...
                    Coefficients_Set( (1 + approx_factor*sum(prod(Book(1,:),2)) + fact_1*sum(prod(Book(init_index:(Extraction_Number-1),:),2))  ) : ...
                                    ( approx_factor*sum(prod(Book(1,:),2)) + fact_2*sum(prod(Book(init_index:Extraction_Number,:),2)) ) );
            elseif strcmp(Detail_Type,'H')
                % set_ranges_details_H = @(book_ , iter_) (1 + 3*sum(prod(book_(1:(iter_-1),:),2)) ) : sum(prod(book_(1:iter_,:),2));
                fact_3 = 0;
                fact_4 = 1;
                coefficients = ...
                    Coefficients_Set( (1 + approx_factor*sum(prod(Book(1,:),2)) + fact_1*sum(prod(Book(init_index:(Extraction_Number-1),:),2)) + ...
                                       fact_3*sum(prod(Book(init_index:Extraction_Number,:),2))) : ...
                                       ( approx_factor*sum(prod(Book(1,:),2)) + fact_2*sum(prod(Book(init_index:Extraction_Number-1,:),2)) + ...
                                       fact_4*sum(prod(Book(Extraction_Number,:),2))) );
            elseif strcmp(Detail_Type,'V')
                % set_ranges_details_V = @(book_ , iter_) (1 + 3*sum(prod(book_(1:(iter_-1),:),2)) +   sum(prod(book_(1:(iter_-1),:),2) ) : ...
%                                             2*sum(prod(book_(1:iter_,:),2)))
                fact_3 = 1;
                fact_4 = 2;
                coefficients = ...
                    Coefficients_Set( (1 + approx_factor*sum(prod(Book(1,:),2)) + fact_1*sum(prod(Book(init_index:(Extraction_Number-1),:),2)) + ...
                                       fact_3*sum(prod(Book(Extraction_Number,:),2))) : ...
                                       ( approx_factor*sum(prod(Book(1,:),2)) + fact_2*sum(prod(Book(init_index:Extraction_Number-1,:),2)) + ...
                                       fact_4*sum(prod(Book(Extraction_Number,:),2))) );
            elseif strcmp(Detail_Type,'D')
                % set_ranges_details_D = @(book_ , iter_) (1 + 3*sum(prod(book_(1:(iter_-1),:),2)) + 2*sum(prod(book_(1:(iter_-1),:),2)) : ...
%                                             3*sum(prod(book_(1:iter_,:),2)));
                fact_3 = 2;
                fact_4 = 3;
                coefficients = ...
                    Coefficients_Set( (1 + approx_factor*sum(prod(Book(1,:),2)) + fact_1*sum(prod(Book(init_index:(Extraction_Number-1),:),2)) + ...
                                       fact_3*sum(prod(Book(Extraction_Number,:),2))) : ...
                                       ( approx_factor*sum(prod(Book(1,:),2)) + fact_2*sum(prod(Book(init_index:Extraction_Number-1,:),2)) + ...
                                       fact_4*sum(prod(Book(Extraction_Number,:),2))) );

            else
                coefficients = [];
            end
        case 3
    end

end