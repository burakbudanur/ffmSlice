clear all
clc

tw1 = [4.35443487e-01,   2.55718139e-16,   3.55754017e-01, 3.69921902e-02, ...
2.65780484e-01,   1.23087691e-01, 8.96316601e-02, 1.99746439e-01, ...
-5.05277891e-02, 9.44414104e-02, -3.97381177e-02, 1.50718227e-02, ...
-1.69870437e-02, -4.78697400e-03, -4.29874538e-03, -6.15070596e-03, ...
3.91009306e-04, -3.00081892e-03, 8.69737447e-04,  -7.44165447e-04, ...
4.23685427e-04, 1.94447859e-07, 1.19333832e-04, 1.01415612e-04, ... 
8.76041620e-06,   5.66807136e-05,  -1.10380807e-05, 1.74490878e-05,  ...
-7.02070045e-06,   2.11126404e-06]';

ftw = @(tw) [vksred(tw); SliceCond(tw)];  % Constrained equation that 
										  % must be satisfied by a 
										  % travelling wave
[tw1, fval, info] = fsolve(ftw, tw1);  % Get it precise

A =  gradVred(tw1);

[V, lambda] = eig(A);
lambda=eig(A);
lambdarealsorted = sort(real(lambda), 'descend');
i1 = find(lambdarealsorted(1)==real(lambda)); i1 = i1(1);
i2 = find(lambdarealsorted(6)==real(lambda)); i2 = i2(1);

v1 = real(V(:, i1));
v2 = imag(V(:, i1));
v3 = real(V(:, i2));

%Gram-Schmidt orthogonalization:
v2o = v2/norm(v2);

v1o = v1 - dot(v1,v2o)*v2o;
v1o = v1o/norm(v1o);

v3o = v3 - dot(v3,v1o)*v1o - dot(v3,v2o)*v2o;
v3o = v3o/norm(v3o);

%Statespace eigenbasis:
e1 = v1o; e2 = v2o; e3 = v3o; ee = [e1, e2, e3];

n = 5;
tf = 115;
i = 1;
xsol=0; uu=0;
hold on
N = 20;

mu = real(lambda(i1));
nu = imag(lambda(i1));

for phi = 0:2*(mu/nu)*pi/N:2*(mu/nu)*pi-2*(mu/nu)*pi/N
    x0 = tw1 + 1e-6*exp(phi)*e1;
    
    [xsol, uu, tt, taut] = ksETDRK4red(x0, tf);
        
    xsolrel = xsol - repmat(tw1', size(xsol, 1), 1);
    xsolssp = ee'*xsolrel';
    xsolssp = xsolssp';
    
    plot3(xsolssp(:,1),xsolssp(:,2),xsolssp(:,3));
    i = i+1;
end

%Mark the other relative equilibrium in this projection:
tw2= [2.32900179177042e-17, 0.419185107640453, 0.177530724935696, ...
0.146282501219190, -0.421516319672415, 0.287443694911592, ...
0.0933649938222870, -0.320489299533689, -0.00580747467566936, ...
0.0162361128812925, 0.0493485521038562, 0.0269378819483284, ...
-0.0298815966184568, 0.00764867193425785, 0.00381954953196330, ...
-0.00509574672330880, -0.000346492806473284, -0.00191533222584679, ...
0.00126630256353163, 0.00123563778819157, -0.000513489812992672, ...
-4.76928288759806e-05, 2.22135119056738e-07, 9.94332807607478e-06, ...
1.31720259229251e-05, -6.81917084273278e-05, 1.31699038813821e-05, ...
2.57228901352395e-05, -3.27286075134879e-06, -1.76547598708647e-06]';
tw2 = LieEl(-pi/2)*tw2;

[tw2, fval, info] = fsolve(ftw, tw2);
tw2rel = tw2 - tw1;

%Project this reqv onto the ssp basis:
tw2ssp = ee'*tw2rel;

rpo0 = [1.92725947347017e-19 3.75553820268048e-01 1.16318112762086e-01 ...
2.75138291542745e-01 -2.61084549544603e-01 -4.98014744888981e-01 ...
2.63249305002884e-02 5.47881940614809e-02 7.77434443829294e-02 ...
6.79067693226824e-02 -4.61139194205068e-02 -3.40247311597783e-02 ...
1.23863730919249e-03 2.80026309458154e-03 7.32296341258052e-03 ...
1.65527969386703e-03 -2.77682791610979e-03 -5.23679491319270e-04 ...
-8.63573372755386e-05 1.44474346812294e-04 3.78079166279577e-04 ...
-9.75826098379176e-05 -1.10511584069682e-04 2.64149399587028e-05 ... 
-6.47700765741695e-06 1.26640502043697e-05 1.26600135032373e-05 ...
-1.08301101019014e-05 -2.79469850392735e-06 1.32808096113333e-06]';

Trpo = 33.5010341720007; %rpo period

rpo0 = LieEl(-pi/2)*rpo0; %Bring the rpo onto the slice

[xsolred, uured, ttred, tautred] = ksETDRK4red(rpo0, 2*Trpo);
    
xsolrel = xsolred - repmat(tw1', size(xsolred, 1), 1);
xsolssp = ee'*xsolrel';
xsolssp = xsolssp';
plot3(xsolssp(:,1),xsolssp(:,2),xsolssp(:,3), 'r', 'LineWidth', 2);

plot3(0, 0, 0, '.m', 'MarkerSize', 10)
text(0,0,0, 'TW_1', 'FontSize', 25, 'Color', 'm')
plot3(tw2ssp(1), tw2ssp(2), tw2ssp(3), '.g', 'MarkerSize', 10)
text(tw2ssp(1), tw2ssp(2), tw2ssp(3), 'TW_2', 'FontSize', 25, 'Color', 'g')
    
xlabel('v_1')
ylabel('v_2')
zlabel('v_3')
