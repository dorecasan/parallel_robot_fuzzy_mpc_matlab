function cineq = myIneqConFunction(X,U,e,data)
p = data.PredictionHorizon;
Ts = data.Ts;
U1 = U(1:p,data.MVIndex(1));
U2 = U(1:p,data.MVIndex(2));
U3 = U(1:p,data.MVIndex(3));
U4 = U(1:p,data.MVIndex(4));
X1 = X(2:p+1,1);X5 = X(2:p+1,5);
X2 = X(2:p+1,2);X6 = X(2:p+1,6);
X3 = X(2:p+1,3);X7 = X(2:p+1,7);
X4 = X(2:p+1,4);X8 = X(2:p+1,8);
%------------ For BSMC -------------------------------------------------
k1 = 1; k2 =1; k3 =1; k4 = 1 ; k5 =1;                               
a = 0.5 ;    ep = [1 1 1 0.1]  ;                     % Some static parameters
mp = 15; m1 =0.5; m2 = 3; mdc =3; Ipy = 1; Ipx = 1; Igamma = 1;
m11 = m2 + mdc + mp/9 + Ipy/(4*a^2); m21 = mp/9;
m22 = m2+mdc+mp/9+ Ipx/(12*a^2); m33 = m2 +mdc +mp/9; m44 = Igamma;
M = [m11 m21 m21 0; m21 m22 m21 0; m21 m21 m33 0; 0 0 0 m44];
invM = inv(M);
c1 = 15*Ipy/(4*a^4) ; c2 = 5*Ipx/(12*a^4);
external_force = [U(1,5);U(1,6);U(1,7);0];
D = ((m2+mp/3)*9.81)*[1;1;1;0]+external_force;

Xref = data.References;
X10 = X(1,1); X20 = X(1,2); X30 = X(1,3); X40 =X(1,4);
dX10 = X(1,5); dX20 = X(1,6); dX30 = X(1,7); dX40 =X(1,8);
q0 = [X10;X20;X30;X40];
dq0 = [dX10;dX20;dX30;dX40];
qd0 = [Xref(1,1);Xref(1,2);Xref(1,3);Xref(1,4)];
dqd0 = [Xref(1,5);Xref(1,6);Xref(1,7);Xref(1,8)];
dqd1 = [Xref(2,5);Xref(2,6);Xref(2,7);Xref(2,8)];
ddqd0 = (dqd1 - dqd0)/Ts;
C0 = [c1*dX10*dX10;c2*dX20*dX20; 0;0];
e1 = k4*(q0 -qd0); alpha = dqd0-k1*e1; e2 = k5*(dq0 - alpha); dalpha  = ddqd0-k1*(dq0-dqd0);
se2 = sign(e2);
VP = -e1'*e2 -k2*e2'*se2 - k3*e2'*e2 +e2'*(invM*(C0+D)+dalpha);
eM = e2'*invM;
U10 = U1(1) ; U20 = U2(1) ; U30 = U3(1) ; U40 = U4(1);

%-----------------------------------------------------------------------
Smax = 0.7; Smin = 0.3;
Umax =  160; Tmax = 20;
cineq = [U1-Umax;U2-Umax;U3-Umax;U4-Tmax;-U1-Umax;-U2-Umax;-U3-Umax;-U4-Tmax;
         X1-Smax;X2-Smax;X3-Smax;X4-pi/2;-X1+Smin;-X2+Smin;-X3+Smin;-X4-pi/2;
         eM(1)*U10+eM(2)*U20+eM(3)*U30+eM(4)*U40 - VP]; % contraction condition after simplifying
end