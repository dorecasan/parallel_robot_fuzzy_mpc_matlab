function y = cal_inv(u)

%% input

q = u(1:4);

l1 = q(1);
l2 = q(2);
l3 = q(3);
gamma = q(4);

% l1 = 50;
% l2 = 50;
% l3 = 50;
% gamma = 0;

% if l1 > 300
%     l1 = 300;
% end
% if l1 < 10
%     l1 = 10;  
% end
% 
% if l2 > 300
%     l2 = 300;
% end
% if l2 < 10
%     l2 = 10;  
% end
% 
% if l2 > 300
%     l2 = 300;
% end
% if l3 < 10
%     l3 = 10;  
% end

global X az a

A1 = [a*sin(pi/6) a*cos(pi/6) 0]';
A2 = [-a 0 0]';
A3 = [a*sin(pi/6) -a*cos(pi/6) 0]';


B1 = [a*sin(pi/6) a*cos(pi/6) 0]';
B2 = [-a 0 0]';
B3 = [a*sin(pi/6) -a*cos(pi/6) 0]';



X = [150;0;0];

for i = 1:100

    Pz = X(1);    
    alpha = X(2);
    beta = X(3);

Rp = [cos(beta)*cos(gamma) cos(gamma)*sin(alpha)*sin(beta)-cos(alpha)*sin(gamma) sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)*sin(beta);
      cos(beta)*sin(gamma) cos(alpha)*cos(gamma)+sin(alpha)*sin(beta)*sin(gamma) cos(alpha)*sin(beta)*sin(gamma)-cos(gamma)*sin(alpha);
      -sin(beta)                    cos(beta)*sin(alpha)                                        cos(alpha)*cos(beta)                   ];
Ra = [cos(gamma) -sin(gamma) 0;
      sin(gamma) cos(gamma)  0; 
           0         0       1];
       


Oa = [0 0 az]';
P  = [0 0 Pz]';

f1 = (Rp*B1 + P - Ra*A1 - Oa)'*(Rp*B1 + P - Ra*A1 - Oa) - l1^2  ;
f2 = (Rp*B2 + P - Ra*A2 - Oa)'*(Rp*B2 + P - Ra*A2 - Oa) - l2^2  ;
f3 = (Rp*B3 + P - Ra*A3 - Oa)'*(Rp*B3 + P - Ra*A3 - Oa) - l3^2  ;

df1z = Pz - az + conj(Pz) - conj(az) - (a*sin(conj(beta)))/2 - (conj(a)*sin(beta))/2 + (3^(1/2)*conj(a)*cos(beta)*sin(alpha))/2 + (3^(1/2)*a*cos(conj(beta))*sin(conj(alpha)))/2;
 
 
df1_alpha = (3^(1/2)*a*(cos(conj(gamma))*sin(conj(alpha)) - cos(conj(alpha))*sin(conj(beta))*sin(conj(gamma)))*((conj(a)*sin(gamma))/2 + (3^(1/2)*conj(a)*cos(gamma))/2 - (conj(a)*cos(beta)*cos(gamma))/2 - (3^(1/2)*conj(a)*(cos(alpha)*cos(gamma) + sin(alpha)*sin(beta)*sin(gamma)))/2))/2 - (3^(1/2)*a*(sin(conj(alpha))*sin(conj(gamma)) + cos(conj(alpha))*cos(conj(gamma))*sin(conj(beta)))*((conj(a)*cos(gamma))/2 - (3^(1/2)*conj(a)*sin(gamma))/2 - (conj(a)*cos(beta)*cos(gamma))/2 + (3^(1/2)*conj(a)*(cos(alpha)*sin(gamma) - cos(gamma)*sin(alpha)*sin(beta)))/2))/2 - (3^(1/2)*conj(a)*(sin(alpha)*sin(gamma) + cos(alpha)*cos(gamma)*sin(beta))*((a*cos(conj(gamma)))/2 - (a*cos(conj(beta))*cos(conj(gamma)))/2 + (3^(1/2)*a*(cos(conj(alpha))*sin(conj(gamma)) - cos(conj(gamma))*sin(conj(alpha))*sin(conj(beta))))/2 - (3^(1/2)*a*sin(conj(gamma)))/2))/2 + (3^(1/2)*conj(a)*(cos(gamma)*sin(alpha) - cos(alpha)*sin(beta)*sin(gamma))*((a*sin(conj(gamma)))/2 - (a*cos(conj(beta))*cos(conj(gamma)))/2 - (3^(1/2)*a*(cos(conj(alpha))*cos(conj(gamma)) + sin(conj(alpha))*sin(conj(beta))*sin(conj(gamma))))/2 + (3^(1/2)*a*cos(conj(gamma)))/2))/2 + (3^(1/2)*a*cos(conj(alpha))*cos(conj(beta))*(conj(Pz) - conj(az) - (conj(a)*sin(beta))/2 + (3^(1/2)*conj(a)*cos(beta)*sin(alpha))/2))/2 + (3^(1/2)*conj(a)*cos(alpha)*cos(beta)*(Pz - az - (a*sin(conj(beta)))/2 + (3^(1/2)*a*cos(conj(beta))*sin(conj(alpha)))/2))/2;
 
 
df1_beta = ((conj(a)*cos(gamma)*sin(beta))/2 - (3^(1/2)*conj(a)*cos(beta)*cos(gamma)*sin(alpha))/2)*((a*cos(conj(gamma)))/2 - (a*cos(conj(beta))*cos(conj(gamma)))/2 + (3^(1/2)*a*(cos(conj(alpha))*sin(conj(gamma)) - cos(conj(gamma))*sin(conj(alpha))*sin(conj(beta))))/2 - (3^(1/2)*a*sin(conj(gamma)))/2) + ((conj(a)*cos(gamma)*sin(beta))/2 - (3^(1/2)*conj(a)*cos(beta)*sin(alpha)*sin(gamma))/2)*((a*sin(conj(gamma)))/2 - (a*cos(conj(beta))*cos(conj(gamma)))/2 - (3^(1/2)*a*(cos(conj(alpha))*cos(conj(gamma)) + sin(conj(alpha))*sin(conj(beta))*sin(conj(gamma))))/2 + (3^(1/2)*a*cos(conj(gamma)))/2) - ((a*cos(conj(beta)))/2 + (3^(1/2)*a*sin(conj(alpha))*sin(conj(beta)))/2)*(conj(Pz) - conj(az) - (conj(a)*sin(beta))/2 + (3^(1/2)*conj(a)*cos(beta)*sin(alpha))/2) + ((a*cos(conj(gamma))*sin(conj(beta)))/2 - (3^(1/2)*a*cos(conj(beta))*cos(conj(gamma))*sin(conj(alpha)))/2)*((conj(a)*cos(gamma))/2 - (3^(1/2)*conj(a)*sin(gamma))/2 - (conj(a)*cos(beta)*cos(gamma))/2 + (3^(1/2)*conj(a)*(cos(alpha)*sin(gamma) - cos(gamma)*sin(alpha)*sin(beta)))/2) + ((a*cos(conj(gamma))*sin(conj(beta)))/2 - (3^(1/2)*a*cos(conj(beta))*sin(conj(alpha))*sin(conj(gamma)))/2)*((conj(a)*sin(gamma))/2 + (3^(1/2)*conj(a)*cos(gamma))/2 - (conj(a)*cos(beta)*cos(gamma))/2 - (3^(1/2)*conj(a)*(cos(alpha)*cos(gamma) + sin(alpha)*sin(beta)*sin(gamma)))/2) - ((conj(a)*cos(beta))/2 + (3^(1/2)*conj(a)*sin(alpha)*sin(beta))/2)*(Pz - az - (a*sin(conj(beta)))/2 + (3^(1/2)*a*cos(conj(beta))*sin(conj(alpha)))/2);
 
 
df2z = Pz - az + conj(Pz) - conj(az) + a*sin(conj(beta)) + conj(a)*sin(beta);
 
 
df2_alpha = 0;
 
 
df2_beta = a*cos(conj(beta))*(conj(Pz) - conj(az) + conj(a)*sin(beta)) + conj(a)*cos(beta)*(Pz - az + a*sin(conj(beta))) + a*cos(conj(gamma))*sin(conj(beta))*(conj(a)*cos(gamma) - conj(a)*cos(beta)*cos(gamma)) + conj(a)*cos(gamma)*sin(beta)*(a*cos(conj(gamma)) - a*cos(conj(beta))*cos(conj(gamma))) + a*cos(conj(gamma))*sin(conj(beta))*(conj(a)*sin(gamma) - conj(a)*cos(beta)*cos(gamma)) + conj(a)*cos(gamma)*sin(beta)*(a*sin(conj(gamma)) - a*cos(conj(beta))*cos(conj(gamma)));
 
 
df3z = Pz - az + conj(Pz) - conj(az) - (a*sin(conj(beta)))/2 - (conj(a)*sin(beta))/2 - (3^(1/2)*conj(a)*cos(beta)*sin(alpha))/2 - (3^(1/2)*a*cos(conj(beta))*sin(conj(alpha)))/2;
 
 
df3_alpha = (3^(1/2)*a*(sin(conj(alpha))*sin(conj(gamma)) + cos(conj(alpha))*cos(conj(gamma))*sin(conj(beta)))*((conj(a)*cos(gamma))/2 + (3^(1/2)*conj(a)*sin(gamma))/2 - (conj(a)*cos(beta)*cos(gamma))/2 - (3^(1/2)*conj(a)*(cos(alpha)*sin(gamma) - cos(gamma)*sin(alpha)*sin(beta)))/2))/2 - (3^(1/2)*a*(cos(conj(gamma))*sin(conj(alpha)) - cos(conj(alpha))*sin(conj(beta))*sin(conj(gamma)))*((conj(a)*sin(gamma))/2 - (3^(1/2)*conj(a)*cos(gamma))/2 - (conj(a)*cos(beta)*cos(gamma))/2 + (3^(1/2)*conj(a)*(cos(alpha)*cos(gamma) + sin(alpha)*sin(beta)*sin(gamma)))/2))/2 + (3^(1/2)*conj(a)*(sin(alpha)*sin(gamma) + cos(alpha)*cos(gamma)*sin(beta))*((a*cos(conj(gamma)))/2 - (a*cos(conj(beta))*cos(conj(gamma)))/2 - (3^(1/2)*a*(cos(conj(alpha))*sin(conj(gamma)) - cos(conj(gamma))*sin(conj(alpha))*sin(conj(beta))))/2 + (3^(1/2)*a*sin(conj(gamma)))/2))/2 - (3^(1/2)*conj(a)*(cos(gamma)*sin(alpha) - cos(alpha)*sin(beta)*sin(gamma))*((a*sin(conj(gamma)))/2 - (a*cos(conj(beta))*cos(conj(gamma)))/2 + (3^(1/2)*a*(cos(conj(alpha))*cos(conj(gamma)) + sin(conj(alpha))*sin(conj(beta))*sin(conj(gamma))))/2 - (3^(1/2)*a*cos(conj(gamma)))/2))/2 + (3^(1/2)*a*cos(conj(alpha))*cos(conj(beta))*(conj(az) - conj(Pz) + (conj(a)*sin(beta))/2 + (3^(1/2)*conj(a)*cos(beta)*sin(alpha))/2))/2 + (3^(1/2)*conj(a)*cos(alpha)*cos(beta)*(az - Pz + (a*sin(conj(beta)))/2 + (3^(1/2)*a*cos(conj(beta))*sin(conj(alpha)))/2))/2;
 
 
df3_beta = ((conj(a)*cos(gamma)*sin(beta))/2 + (3^(1/2)*conj(a)*cos(beta)*cos(gamma)*sin(alpha))/2)*((a*cos(conj(gamma)))/2 - (a*cos(conj(beta))*cos(conj(gamma)))/2 - (3^(1/2)*a*(cos(conj(alpha))*sin(conj(gamma)) - cos(conj(gamma))*sin(conj(alpha))*sin(conj(beta))))/2 + (3^(1/2)*a*sin(conj(gamma)))/2) + ((conj(a)*cos(gamma)*sin(beta))/2 + (3^(1/2)*conj(a)*cos(beta)*sin(alpha)*sin(gamma))/2)*((a*sin(conj(gamma)))/2 - (a*cos(conj(beta))*cos(conj(gamma)))/2 + (3^(1/2)*a*(cos(conj(alpha))*cos(conj(gamma)) + sin(conj(alpha))*sin(conj(beta))*sin(conj(gamma))))/2 - (3^(1/2)*a*cos(conj(gamma)))/2) + ((a*cos(conj(beta)))/2 - (3^(1/2)*a*sin(conj(alpha))*sin(conj(beta)))/2)*(conj(az) - conj(Pz) + (conj(a)*sin(beta))/2 + (3^(1/2)*conj(a)*cos(beta)*sin(alpha))/2) + ((a*cos(conj(gamma))*sin(conj(beta)))/2 + (3^(1/2)*a*cos(conj(beta))*cos(conj(gamma))*sin(conj(alpha)))/2)*((conj(a)*cos(gamma))/2 + (3^(1/2)*conj(a)*sin(gamma))/2 - (conj(a)*cos(beta)*cos(gamma))/2 - (3^(1/2)*conj(a)*(cos(alpha)*sin(gamma) - cos(gamma)*sin(alpha)*sin(beta)))/2) + ((a*cos(conj(gamma))*sin(conj(beta)))/2 + (3^(1/2)*a*cos(conj(beta))*sin(conj(alpha))*sin(conj(gamma)))/2)*((conj(a)*sin(gamma))/2 - (3^(1/2)*conj(a)*cos(gamma))/2 - (conj(a)*cos(beta)*cos(gamma))/2 + (3^(1/2)*conj(a)*(cos(alpha)*cos(gamma) + sin(alpha)*sin(beta)*sin(gamma)))/2) + ((conj(a)*cos(beta))/2 - (3^(1/2)*conj(a)*sin(alpha)*sin(beta))/2)*(az - Pz + (a*sin(conj(beta)))/2 + (3^(1/2)*a*cos(conj(beta))*sin(conj(alpha)))/2);

f = [f1 f2 f3]';

df = [df1z df1_alpha df1_beta;
      df2z df2_alpha df2_beta;   
      df3z df3_alpha df3_beta];

  
X = X - (df^(-1))*f;

end
 
X(1) = abs(X(1));
y = X;

end

