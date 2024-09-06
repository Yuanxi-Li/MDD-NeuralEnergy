classdef IonChannel
    %create a class for ion channels
    %including ion types, current and g(if there is)
    properties
        Neurontype
        Iontype
        current
        g
    end
    methods
        function obj = IonChannel(Neurontype, Iontype, current, g)
            obj.Neurontype = Neurontype;
            obj.Iontype = Iontype;
            obj.current = current;
            obj.g = g;
        end
    end
end
            
            
            
            
            