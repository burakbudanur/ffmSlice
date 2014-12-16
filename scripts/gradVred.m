function Ared = gradVred(x)
%Nazmi Burak Budanur
%Kuramoto-Sivashinsky reduced stability matrix

A = gradV(x);
T = Lg();
xp = slicep();
tp = T*xp;
tx = T*x;
vx = vel(x);

Ared = A - (tx*((A'*(tx'*tp) - T'*(vx'*tp))*tp)')/(tx'*tp)^2 - ((vx'*tp)/(tx'*tp))*T;
