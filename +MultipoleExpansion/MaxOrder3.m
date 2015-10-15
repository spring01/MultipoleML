classdef MaxOrder3 < MultipoleExpansion.Base
    
    methods (Access = protected)
        
        function res = BlockWithMaxOrder0(obj, vecAB)
            res = [ ...
                obj.Block00(vecAB); ...     % 00
                obj.Block10(vecAB); ...     % 10
                obj.Block20(vecAB); ...     % 20
                obj.Block30(vecAB)];        % 30
        end
        
        function res = BlockWithMaxOrder1(obj, vecAB)
            res = [ ...
                obj.Block10(-vecAB)'; ...   % 01
                obj.Block11(vecAB); ...     % 11
                obj.Block21(vecAB); ...     % 21
                obj.Block31(vecAB)];        % 31
        end
        
        function res = BlockWithMaxOrder2(obj, vecAB)
            res = [ ...
                obj.Block20(-vecAB)'; ...   % 02
                obj.Block21(-vecAB)'; ...   % 12
                obj.Block22(vecAB); ...     % 22
                obj.Block32(vecAB)];        % 32
        end
        
        function res = BlockWithMaxOrder3(obj, vecAB)
            res = [ ...
                obj.Block30(-vecAB)'; ...   % 03
                obj.Block31(-vecAB)'; ...   % 13
                obj.Block32(-vecAB)'; ...   % 23
                obj.Block33(vecAB)];        % 33
        end
        
        function res = BlockWithMaxOrder4(obj, vecAB)
            res = [ ...
                obj.Block40(-vecAB)'; ...   % 04
                obj.Block41(-vecAB)'; ...   % 14
                obj.Block42(-vecAB)'; ...   % 24
                obj.Block43(-vecAB)'];      % 34
        end
        
        function res = BlockWithMaxOrder5(obj, vecAB)
            res = [ ...
                obj.Block50(-vecAB)'; ...   % 05
                obj.Block51(-vecAB)'; ...   % 15
                obj.Block52(-vecAB)'; ...   % 25
                obj.Block53(-vecAB)'];      % 35
        end
        
    end
    
end
