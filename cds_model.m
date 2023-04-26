function z_out = cds_model(x,u)

a = 0.5 ;                            % Some static parameters
mp = 15; m1 =0.5; m2 = 3; mdc =3; Ipy = 1; Ipx = 1; Igamma = 1;
m11 = m2 + mdc + mp/9 + Ipy/(4*a^2); m21 = mp/9;
m22 = m2+mdc+mp/9+ Ipx/(12*a^2); m33 = m2 +mdc +mp/9; m44 = Igamma;
M = [m11 m21 m21 0; m21 m22 m21 0; m21 m21 m33 0; 0 0 0 m44];
invM = inv(M);
c1 = 15*Ipy/(4*a^4) ; c2 = 5*Ipx/(12*a^4);
external_force = [u(5);u(6);u(7);u(8)];
D = ((m2+mp/3)*9.81)*[1;1;1;0]+external_force;

dl1 = x(5);
dl2 = x(6);
dl3 = x(7);
dgamma = x(8);
ddl1 =invM(1,1)*(u(1)-D(1)-c1*dl1*dl1)+invM(1,2)*(u(2) - D(2)-c2*dl2*dl2)+invM(1,3)*(u(3)-D(3))+invM(1,4)*(u(4)-D(4));
ddl2 =invM(2,1)*(u(1)-D(1)-c1*dl1*dl1)+invM(2,2)*(u(2) - D(2)-c2*dl2*dl2)+invM(2,3)*(u(3)-D(3))+invM(2,4)*(u(4)-D(4));
ddl3 =invM(3,1)*(u(1)-D(1)-c1*dl1*dl1)+invM(3,2)*(u(2) - D(2)-c2*dl2*dl2)+invM(3,3)*(u(3)-D(3))+invM(3,4)*(u(4)-D(4));
ddgamma =invM(4,4)*(u(4)-D(4));

z_out = [dl1;dl2;dl3;dgamma;ddl1;ddl2;ddl3;ddgamma];
end