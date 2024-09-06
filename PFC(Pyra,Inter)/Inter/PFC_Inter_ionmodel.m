classdef PFC_Inter_ionmodel
% Create Class for ion in PFC Pyra:
% Name, gmax, phi, rev
    properties
        Name
        gmax
        phi
        rev
    end
    
    methods
        function obj = PFC_Inter_ionmodel(Name, gmax, phi, rev)
            obj.Name = Name;
            obj.gmax = gmax;
            obj.phi = phi;
            obj.rev = rev;
        end
    end
end
            