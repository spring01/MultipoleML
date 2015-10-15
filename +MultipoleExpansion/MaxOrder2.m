classdef MaxOrder2 < MultipoleExpansion.Base
    
    methods (Access = protected)
        
        function res = BlockWithMaxOrder0(obj, vecAB)
            res = [ ...
                obj.Block00(vecAB); ...     % 00
                obj.Block10(vecAB); ...     % 01
                obj.Block20(vecAB)];        % 02
        end
        
        function res = BlockWithMaxOrder1(obj, vecAB)
            res = [ ...
                obj.Block10(-vecAB)'; ...   % 01
                obj.Block11(vecAB); ...     % 11
                obj.Block21(vecAB)];        % 21
        end
        
        function res = BlockWithMaxOrder2(obj, vecAB)
            res = [ ...
                obj.Block20(-vecAB)'; ...   % 02
                obj.Block21(-vecAB)'; ...   % 12
                obj.Block22(vecAB)];        % 22
        end
        
        function res = BlockWithMaxOrder3(obj, vecAB)
            res = [ ...
                obj.Block30(-vecAB)'; ...   % 03
                obj.Block31(-vecAB)'; ...   % 13
                obj.Block32(-vecAB)'];      % 23
        end
        
        function res = BlockWithMaxOrder4(obj, vecAB)
            res = [ ...
                obj.Block40(-vecAB)'; ...   % 04
                obj.Block41(-vecAB)'; ...   % 14
                obj.Block42(-vecAB)'];      % 24
        end
        
        function res = BlockWithMaxOrder5(obj, vecAB)
            res = [ ...
                obj.Block50(-vecAB)'; ...   % 05
                obj.Block51(-vecAB)'; ...   % 15
                obj.Block52(-vecAB)'];      % 25
        end
        
    end
    
end
