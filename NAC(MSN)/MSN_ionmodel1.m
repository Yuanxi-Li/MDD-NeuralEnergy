classdef MSN_ionmodel1
% Create Class for ion except Ca2+:
% NaF, NaP, KAf, KAs, KIR, KRP, Leak, BKKCa, SKKCa
% Name, gmax, a, mVhalf, hvhalf, mk, hk, rev
    properties
        Name
        gmax
        a
        mVhalf
        hVhalf
        mk
        hk
        rev
    end
    
    methods
        function obj = MSN_ionmodel1(Name, gmax, a, mVhalf, hVhalf, mk, hk, rev)
            obj.Name = Name;
            obj.gmax = gmax;
            obj.a = a;
            obj.mVhalf = mVhalf;
            obj.hVhalf = hVhalf;
            obj.mk = mk;
            obj.hk = hk;
            obj.rev = rev;
        end
    end
end
            