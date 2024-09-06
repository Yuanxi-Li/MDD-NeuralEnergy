classdef MSNComponent
% Create Class for components:
% Type: Soma, Proximal, Middle, Distal
% Length, Diameter: micron, 10^-6 meter
% Capacitance: micron-F/(cm^2)'''
    properties
        Type
        Length
        Diameter
        Capacitance
    end
    
    methods
        function obj = MSNComponent(Type,Length,Diameter,Capacitance)
            obj.Type = Type;
            obj.Length = Length;
            obj.Diameter = Diameter;
            obj.Capacitance = Capacitance;
        end
    end
end
            