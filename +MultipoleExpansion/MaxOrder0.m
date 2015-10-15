classdef MaxOrder0 < MultipoleExpansion.Base
    
    methods (Access = protected)
        
        function res = BlockWithMaxOrder0(obj, vecAB)
            res = obj.Block00(vecAB);
        end
        
        function res = BlockWithMaxOrder1(obj, vecAB)
            res = obj.Block10(-vecAB)';
        end
        
        function res = BlockWithMaxOrder2(obj, vecAB)
            res = obj.Block20(-vecAB)';
        end
        
        function res = BlockWithMaxOrder3(obj, vecAB)
            res = obj.Block30(-vecAB)';
        end
        
        function res = BlockWithMaxOrder4(obj, vecAB)
            res = obj.Block40(-vecAB)';
        end
        
        function res = BlockWithMaxOrder5(obj, vecAB)
            res = obj.Block50(-vecAB)';
        end
        
    end
    
end
