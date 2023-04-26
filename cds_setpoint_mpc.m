function [sys, x0, str, ts] = cds_setpoint_mpc(t,x,u,flag)
    switch flag
    case 0
        [sys, x0, str,ts] = mdlInitializeSizes;
    case 1
        sys = mdlDerivatives(t,x,u);
    case 3
        sys = mdlOutputs(t,x,u);
    case {2,4,9}
        sys = [];
    otherwise
        error(['Unhandled flag = ', num2str(flag)]);
    end
end
function [sys,x0,str,ts]=mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates = 0;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 12;
    sizes.NumInputs = 0 ;
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0 = [];
    str = [];
    ts = [0 0];
end


function sys = mdlOutputs(t,x,u)
%  l1  = 0.0008*t^3-0.0125*t^2+0.072*t+0.5;
%  l2  = 0.0008*t^3-0.0125*t^2+0.072*t+0.5;
%  l3 =  -0.0008*t^3+0.0175*t^2-0.11*t+0.7;
%  gamma = -0.02*t^3+0.1625*t^2-0.1*t;
%  dl1 = 3*0.0008*t^2-2*0.0125*t+0.072;
%  dl2 = 3*0.0008*t^2-2*0.0125*t+0.072;
%  dl3 =  -3*0.0008*t^2+2*0.0175*t-0.11;
%  dgamma = -3*0.02*t^2+2*0.1625*t-0.1;
%  ddl1 = 6*0.0008*t-2*0.0125;
%  ddl2 = 6*0.0008*t-2*0.0125;
%  ddl3 =  -6*0.0008*t+2*0.0175;
%  ddgamma = -6*0.02*t+2*0.1625;

 l1  = -0.0045*t^3+0.034*t^2-0.025*t+0.45;
 l2  =  -0.0008*t^3+0.0175*t^2-0.11*t+0.7;
 l3 =  0.0008*t^3-0.0125*t^2+0.072*t+0.5;
 gamma = 0.004*t^3-0.0625*t^2+0.4*t;
 dl1  = -3*0.0045*t^2+2*0.034*t-0.025;
 dl2  =  -3*0.0008*t^2+2*0.0175*t-0.11;
 dl3 =  3*0.0008*t^2-2*0.0125*t+0.072;
 dgamma = 3*0.004*t^2-2*0.0625*t+0.4;
 ddl1  = -6*0.0045*t+2*0.034;
 ddl2  =  -6*0.0008*t+2*0.0175;
 ddl3 =  6*0.0008*t-2*0.0125;
 ddgamma = 6*0.004*t-2*0.0625;
 
sys(1) = l1;
sys(2) = dl1;
sys(3) = l2;
sys(4) = dl2;

sys(5) = l3;
sys(6) = dl3;

sys(7) = gamma;
sys(8) = dgamma;
sys(9) = ddl1;
sys(10) = ddl2;

sys(11) = ddl3;
sys(12) = ddgamma;

end

