function y = cal_q(u)

%% input

p = u(1:4); 
    Pz      = p(1);
    alpha   = p(2);
    beta    = p(3);
    gamma   = u(4);

a=0.5;az=0.05;

% Pz = 500;
% alpha = 0;
% beta = 0;
% gamma = 0;



A1 = [a*sin(pi/6) a*cos(pi/6) 0]';
A2 = [-a 0 0]';
A3 = [a*sin(pi/6) -a*cos(pi/6) 0]';

B1 = [a*sin(pi/6) a*cos(pi/6) 0]';
B2 = [-a 0 0]';
B3 = [a*sin(pi/6) -a*cos(pi/6) 0]';


Rp = [cos(beta)*cos(gamma) cos(gamma)*sin(alpha)*sin(beta)-cos(alpha)*sin(gamma) sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)*sin(beta);
      cos(beta)*sin(gamma) cos(alpha)*cos(gamma)+sin(alpha)*sin(beta)*sin(gamma) cos(alpha)*sin(beta)*sin(gamma)-cos(gamma)*sin(alpha);
      -sin(beta)                    cos(beta)*sin(alpha)                                        cos(alpha)*cos(beta)                   ];
Ra = [cos(gamma) -sin(gamma) 0;
      sin(gamma) cos(gamma)  0; 
           0         0       1];
Oa = [0 0 az]';
P  = [0 0 Pz]';

A1_r = Ra*A1 + Oa;  
A2_r = Ra*A2 + Oa; 
A3_r = Ra*A3 + Oa;        

B1_r = Rp*B1 + P;
B2_r = Rp*B2 + P;
B3_r = Rp*B3 + P;

AB1 = B1_r - A1_r;
AB2 = B2_r - A2_r;
AB3 = B3_r - A3_r;

l1 = abs(sqrt(AB1'*AB1));
l2 = abs(sqrt(AB2'*AB2));
l3 = abs(sqrt(AB3'*AB3));



 
 y = [l1;l2;l3;gamma];
%y = AB3;

end

