classdef MaxOrder5 < MultipoleExpansion.Base
    
    methods (Access = protected)
        
        function res = BlockWithMaxOrder0(obj, vecAB)
            res = [ ...
                obj.Block00(vecAB); ...     % 00
                obj.Block10(vecAB); ...     % 10
                obj.Block20(vecAB); ...     % 20
                obj.Block30(vecAB); ...     % 30
                obj.Block40(vecAB); ...     % 40
                obj.Block50(vecAB)];        % 50
        end
        
        function res = BlockWithMaxOrder1(obj, vecAB)
            res = [ ...
                obj.Block10(-vecAB)'; ...   % 01
                obj.Block11(vecAB); ...     % 11
                obj.Block21(vecAB); ...     % 21
                obj.Block31(vecAB); ...     % 31
                obj.Block41(vecAB); ...     % 41
                obj.Block51(vecAB)];        % 51
        end
        
        function res = BlockWithMaxOrder2(obj, vecAB)
            res = [ ...
                obj.Block20(-vecAB)'; ...   % 02
                obj.Block21(-vecAB)'; ...   % 12
                obj.Block22(vecAB); ...     % 22
                obj.Block32(vecAB); ...     % 32
                obj.Block42(vecAB); ...     % 42
                obj.Block52(vecAB)];        % 52
        end
        
        function res = BlockWithMaxOrder3(obj, vecAB)
            res = [ ...
                obj.Block30(-vecAB)'; ...   % 03
                obj.Block31(-vecAB)'; ...   % 13
                obj.Block32(-vecAB)'; ...   % 23
                obj.Block33(vecAB); ...     % 33
                obj.Block43(vecAB); ...     % 43
                obj.Block53(vecAB)];        % 53
        end
        
        function res = BlockWithMaxOrder4(obj, vecAB)
            res = [ ...
                obj.Block40(-vecAB)'; ...   % 04
                obj.Block41(-vecAB)'; ...   % 14
                obj.Block42(-vecAB)'; ...   % 24
                obj.Block43(-vecAB)'; ...   % 34
                obj.Block44(vecAB); ...     % 44
                obj.Block54(vecAB)];        % 54
        end
        
        function res = BlockWithMaxOrder5(obj, vecAB)
            res = [ ...
                obj.Block50(-vecAB)'; ...   % 05
                obj.Block51(-vecAB)'; ...   % 15
                obj.Block52(-vecAB)'; ...   % 25
                obj.Block53(-vecAB)'; ...   % 35
                obj.Block54(-vecAB)'; ...   % 45
                obj.Block55(vecAB)];        % 55
        end
        
    end
    
end
