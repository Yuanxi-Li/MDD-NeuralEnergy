classdef SynapseChannel
    %create a class for ion channels
    %including ion types, current and g(if there is)
    properties
        Synapsetype
        Presynapse
        Postsynapse
        gmax
        alpha
        beta
        rev
    end
    methods
        function obj = SynapseChannel(Synapsetype, Presynapse, Postsynapse, gmax,...
                alpha, beta, rev)
            obj.Synapsetype = Synapsetype;
            obj.Presynapse = Presynapse;
            obj.Postsynapse = Postsynapse;
            obj.gmax = gmax;
            obj.alpha = alpha;
            obj.beta = beta;
            obj.rev = rev;
        end
    end
end
            
            
            
            
            