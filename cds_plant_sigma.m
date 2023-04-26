function [sys, x0, str, ts] = cds_plant_sigma(t,x,u,flag)
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
    sizes.NumContStates = 8;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 8;
    sizes.NumInputs =10 ;
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [0.45 0.7 0.5 0 0 0 0 0];
    str = [];
    ts = [];
end

function sys = mdlDerivatives(t,x,u) 
a = 0.5 ;mp = 15; m1 =0.5; m2 = 3; mdc =3; Ipy = 1; Ipx = 1; Igamma = 1;
m11 = m2 + mdc + mp/9 + Ipy/(4*a^2); m21 = mp/9;
m22 = m2+mdc+mp/9+ Ipx/(12*a^2); m33 = m2 +mdc +mp/9; m44 = Igamma;
M = [m11 m21 m21 0; m21 m22 m21 0; m21 m21 m33 0; 0 0 0 m44];
invM = inv(M);
c11 = 15*Ipy/(4*a^4) ; c22 = 5*Ipx/(12*a^4);
C = zeros(4,4) ; C(1,1) = c11*x(5) ; C(2,2) =c22*x(6);
ef = u(10);
D = ((m2+mp/3)*9.81+ef)*[1;1;1;0];
q = [x(1);x(2);x(3);x(4)]; dq = [x(5);x(6);x(7);x(8)];
F = [u(1);u(2);u(3);u(4)]; noise = [u(9);u(9);u(9);0];
ddq = invM*(F - C*dq -D - noise);
sys(1) = x(5);
sys(2) = x(6);
sys(3) = x(7);
sys(4)= x(8);
sys(5) = ddq(1);
sys(6) = ddq(2);
sys(7) = ddq(3);
sys(8) = ddq(4);

end

function sys = mdlOutputs(t,x,u)
sys(1)=x(1);
sys(2)=x(2);
sys(3)=x(3);
sys(4)=x(4);
sys(5)=x(5);
sys(6)=x(6);
sys(7)=x(7);
sys(8)=x(8);

end

