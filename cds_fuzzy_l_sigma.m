function [sys, x0, str, ts] = cds_fuzzy_l_sigma(t,x,u,flag)
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
    sizes.NumContStates = 4000;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 4;
    sizes.NumInputs =20 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [0*ones(3993,1); 0*ones(7,1)];
    str = [];
    ts = [];
end

function sys = mdlDerivatives(t,x,u)
k1 =15;  c1=1;  c2=1; T1 =800; T2 = 2;
q = [u(1);u(3);u(5);u(7)];
dq = [u(2);u(4);u(6);u(8)];
qr = [u(9);u(11);u(13);u(15)];
dqr = [u(10);u(12);u(14);u(16)];

e1 = c1*(q - qr); de1 = c1*(dq - dqr);
alpha = dqr - k1*e1; 
e2 = c2*(dq- alpha);

l1 = e1(1) ; l2 = e1(2); l3 = e1(3); gamma= e1(4); 

% l_params = [-0.3 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.3 0.02123]; 
% gamma_params = [-pi/6, -pi/9, -pi/18,  0, pi/18 ,pi/9 , pi/6, pi/36];

l_params = [-0.125 -0.1 -0.075 -0.05 -0.025 0 0.025 0.05 0.075 0.1 0.125 0.0125]; 
gamma_params = [-pi/6, -pi/9, -pi/18,  0, pi/18 ,pi/9 , pi/6, pi/36];

ll = length(l_params) -1 ; gl = length(gamma_params) -1 ;
suml = ll^3;
sumg = gl;

for i = 1:1:suml
  thetal(i,1) = x(i); 
  thetal(i+suml,1) = x(i+suml);
  thetal(i+suml*2,1) = x(i+suml*2);
end
for i =1:1:sumg
  thetag(i,1) = x(i+suml*3);
end

FS1 = 0 ;

for i1 = 1:1:ll
    ul1 = gaussmf(l1,[l_params(end),l_params(i1)]);
    for i2 = 1:1:ll
        ul2 = gaussmf(l2,[l_params(end),l_params(i2)]);
        for i3 = 1:1:ll
            ul3 = gaussmf(l3,[l_params(end),l_params(i3)]);                                             
            idx = i3+ ll*(i2-1)+ ll^2*(i1-1);
            FS2(idx) = ul1*ul2*ul3;
            FS1 = FS1 + ul1*ul2*ul3 ;                    
        end
    end
end
  
FSl = FS2/FS1;
clear FS2;
FS1 = 0;
for i = 1:1:gl
    ugamma = gaussmf(gamma,[gamma_params(end),gamma_params(i)]);
    FS2(i) = ugamma;
    FS1 = FS1 + ugamma;
  
end
FSg = FS2/FS1;

z_matrix = zeros(1,suml);

Wl = [FSl,z_matrix,z_matrix;
     z_matrix,FSl,z_matrix;
     z_matrix,z_matrix,FSl];
Wg = FSg;

Wl_dot = T1*Wl'*e2(1:3);
Wg_dot = T2*Wg'*e2(end);
clear Wl;
for i = 1:1:suml
  sys(i) = Wl_dot(i); 
  sys(i+suml) = Wl_dot(i+suml);
  sys(i+suml*2) = Wl_dot(i+suml*2);
end
for i =1:1:sumg
    sys(i+suml*3) = Wg_dot(i);
end
end
%-------------------------------------------------------------------
%---------------------------------------------------------------------
function sys = mdlOutputs(t,x,u)
c1=1;  
mp = 15; m1 =0.5; m2 = 3; mdc =3; Ipy = 1; Ipx = 1; Igamma = 1; a = 0.5;
m11 = m2 + mdc + mp/9 + Ipy/(4*a^2); m21 = mp/9;
m22 = m2+mdc+mp/9+ Ipx/(12*a^2); m33 = m2 +mdc +mp/9; m44 = Igamma;
M = [m11 m21 m21 0; m21 m22 m21 0; m21 m21 m33 0; 0 0 0 m44];

q = [u(1);u(3);u(5);u(7)];
dq = [u(2);u(4);u(6);u(8)];
qr = [u(9);u(11);u(13);u(15)];
dqr = [u(10);u(12);u(14);u(16)];


e1 = c1*(q - qr); de1 = c1*(dq - dqr);


l1 = e1(1) ; l2 = e1(2); l3 = e1(3); gamma= e1(4); 

l_params = [-0.125 -0.1 -0.075 -0.05 -0.025 0 0.025 0.05 0.075 0.1 0.125 0.0125]; 
gamma_params = [-pi/6, -pi/9, -pi/18,  0, pi/18 ,pi/9 , pi/6, pi/36];

% l_params = [-0.3 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.3 0.02123]; 
% gamma_params = [-pi/6, -pi/9, -pi/18,  0, pi/18 ,pi/9 , pi/6, pi/36];

ll = length(l_params) -1 ; gl = length(gamma_params) -1 ;

suml = ll^3;
sumg = gl;

for i = 1:1:suml
  thetal(i,1) = x(i); 
  thetal(i+suml,1) = x(i+suml);
  thetal(i+suml*2,1) = x(i+suml*2);
end
for i =1:1:sumg
  thetag(i,1) = x(i+suml*3);
end

FS1 = 0 ;

for i1 = 1:1:ll
    ul1 = gaussmf(l1,[l_params(end),l_params(i1)]);
    for i2 = 1:1:ll
        ul2 = gaussmf(l2,[l_params(end),l_params(i2)]);
        for i3 = 1:1:ll
            ul3 = gaussmf(l3,[l_params(end),l_params(i3)]);                                             
            idx = i3+ ll*(i2-1)+ ll^2*(i1-1);
            FS2(idx) = ul1*ul2*ul3;
            FS1 = FS1 + ul1*ul2*ul3 ;                    
        end
    end
end
  
FSl = FS2/FS1;
clear FS2;
FS1 = 0;
for i = 1:1:gl
    ugamma = gaussmf(gamma,[gamma_params(end),gamma_params(i)]);
    FS2(i) = ugamma;
    FS1 = FS1 + ugamma;
  
end
FSg = FS2/FS1;

z_matrix = zeros(1,suml);

Wl = [FSl,z_matrix,z_matrix;
     z_matrix,FSl,z_matrix;
     z_matrix,z_matrix,FSl];
Wg = FSg;

sigma = Wl*thetal;
sigma(4) = Wg*thetag;
deltaD = -M*sigma;


sys(1) = deltaD(1);
sys(2) = deltaD(2);
sys(3) = deltaD(3);
sys(4) = deltaD(4);
end

