function [sys, x0, str, ts] = cds_main_bsmc_sigma(t,x,u,flag)
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
    sizes.NumOutputs = 8;
    sizes.NumInputs =24 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [];
    str = [];
    ts = [];
end


function sys = mdlOutputs(t,x,u) 
c1 = 1; c2 = 1; k1 =5;k2=5;k3=5;
a = 0.5 ;mp = 15; m1 =0.5; m2 = 3; mdc =3; Ipy = 1; Ipx = 1; Igamma = 1;
q = [u(1);u(3);u(5);u(7)];
dq = [u(2);u(4);u(6);u(8)];

qr = [u(9);u(11);u(13);u(15)];
dqr = [u(10);u(12);u(14);u(16)];
ddqr = [u(17);u(18);u(19);u(20)];
sigma = [u(21);u(22);u(23);u(24)];

e1 = c1*(q - qr); de1 = c1*(dq - dqr);
alpha = dqr - k1*e1;
dalpha = ddqr - k1*de1;
e2 = c2*(dq- alpha);

m11 = m2 + mdc + mp/9 + Ipy/(4*a^2); m21 = mp/9;
m22 = m2+mdc+mp/9+ Ipx/(12*a^2); m33 = m2 +mdc +mp/9; m44 = Igamma;
M = [m11 m21 m21 0; m21 m22 m21 0; m21 m21 m33 0; 0 0 0 m44];

c11 = 15*Ipy/(4*a^4) ; c22 = 5*Ipx/(12*a^4);
C = zeros(4,4) ; C(1,1) = c11*dq(1) ; C(2,2) =c22*dq(2);
D = (m2+mp/3)*9.81*[1;1;1;0];


deltaD = sigma;
ep = [1;1;1;0.1];
se2 = e2;
for i=1:1:4
    if abs(e2(i)/ep(i)) >1
        se2(i)  = sign(e2(i));
    end
end

Feq = C*dq + D - M*e1 + M*dalpha + deltaD;
Fsw = -M*(k2*se2 + k3*e2);
F = Feq + Fsw;
F_max= 250;
for i=1:1:4
    if abs(F(i)) > F_max
        F(i) = sign(F(i))*F_max;
    end
end

sys(1) = F(1);
sys(2) = F(2);
sys(3) = F(3);
sys(4) = F(4);
sys(5) = deltaD(1);
sys(6) = deltaD(2);
sys(7) = deltaD(3);
sys(8) = deltaD(4);
end

