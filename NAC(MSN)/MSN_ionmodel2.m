classdef MSN_ionmodel2
% Create Class for ion including Ca2+:
% CaL1_2, CaL1_3, CaN, CaQ, CaR, CaT
% Name, Pbar(cm/s), a, mVhalf, hvhalf, mk, hk
    properties
        Name
        Pbar
        a
        mVhalf
        hVhalf
        mk
        hk
    end
    
    methods
        function obj = MSN_ionmodel2(Name, Pbar, a, mVhalf, hVhalf, mk, hk)
            obj.Name = Name;
            obj.Pbar = Pbar;
            obj.a = a;
            obj.mVhalf = mVhalf;
            obj.hVhalf = hVhalf;
            obj.mk = mk;
            obj.hk = hk;
        end
    end
end
            