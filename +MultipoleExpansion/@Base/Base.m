classdef (Abstract) Base < handle
    
    properties (SetAccess = private)
        
        xyz;
        coeffs; % (maxOrder+1)^2 by 1, ordering [00 z y x 20 21s 21c ...]
        
    end
    
    methods
        
        function obj = InitializeFromGDMA(obj, matpsiGDMA, iSite)
            obj.xyz = matpsiGDMA.xyzSites(:, iSite);
            obj.coeffs = ...
                matpsiGDMA.multipoles(1:(matpsiGDMA.limit(iSite)+1)^2, iSite);
        end
        
        function interaction = InteractionWith(multiExpA, multiExpB)
            interaction ...
                = multiExpA.coeffs' ...
                * multiExpA.InteractionMatrixWith(multiExpB) ...
                * multiExpB.coeffs;
        end
        
        function matrix = InteractionMatrixWith(multiExpA, multiExpB)
            xyzA2B = multiExpB.xyz - multiExpA.xyz; % vector A to B
            
            matrix = multiExpA.MatrixWithMaxOrderPlus1( ...
                xyzA2B([3 2 1]) ./ norm(xyzA2B), ... % use normalized zyx
                multiExpB.MaxOrderPlus1());
            
            exponents_matrix ... % -N
                = -repmat((0:multiExpA.MaxOrderPlus1()-1)', 1, multiExpB.MaxOrderPlus1()) ...
                - repmat((1:multiExpB.MaxOrderPlus1()), multiExpA.MaxOrderPlus1(), 1);
            
            bigR_to_minusN_matrix = norm(xyzA2B) .^ exponents_matrix; % R^(-N)
            matrix = matrix .* bigR_to_minusN_matrix( ...
                ceil(sqrt(1:length(multiExpA.coeffs))), ... % mappingA
                ceil(sqrt(1:length(multiExpB.coeffs))) );   % mappingB
        end
        
    end
    
    methods (Access = protected)
        
        res = Block00(~, ~);
        res = Block10(~, vecZYX);
        res = Block11(~, vecZYX);
        res = Block20(~, vecZYX);
        res = Block21(~, vecZYX);
        res = Block22(~, vecZYX);
        res = Block30(~, vecZYX);
        res = Block31(~, vecZYX);
        res = Block32(~, vecZYX);
        res = Block33(~, vecZYX);
        res = Block40(~, vecZYX);
        res = Block41(~, vecZYX);
        res = Block42(~, vecZYX);
        res = Block43(~, vecZYX);
        res = Block44(~, vecZYX);
        res = Block50(~, vecZYX);
        res = Block51(~, vecZYX);
        res = Block52(~, vecZYX);
        res = Block53(~, vecZYX);
        res = Block54(~, vecZYX);
        res = Block55(~, vecZYX);
        
    end
    
    methods (Abstract, Access = protected)
        
        res = BlockWithMaxOrder0(obj, vecAB);
        res = BlockWithMaxOrder1(obj, vecAB);
        res = BlockWithMaxOrder2(obj, vecAB);
        res = BlockWithMaxOrder3(obj, vecAB);
        res = BlockWithMaxOrder4(obj, vecAB);
        res = BlockWithMaxOrder5(obj, vecAB);
        
    end
    
    methods (Access = private)
        
        function maxOrderPlus1 = MaxOrderPlus1(obj)
            maxOrderPlus1 = sqrt(length(obj.coeffs));
        end
        
        res = MatrixWithMaxOrderPlus1(obj, vecAB, maxOrderPlus1);
        
    end
    
end