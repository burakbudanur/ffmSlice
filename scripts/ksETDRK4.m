function [xsol, uu, tt] = ksETDRK4(x0, tf)
%Adaptation of kursiv.m from Kassam ^ Trefethen (2005)

%Spatial grid:
d = 22;
N = 32;
Nh = N/2;
x = d*(1:N)'/N;

%Initial condition:
v = [0; x0(1:2:end-1)+1i*x0(2:2:end); 0; x0(end-1:-2:1)-1i*x0(end:-2:2)];

%Template:
slicep = zeros(N,1);
slicep(2) = 1;
slicep(N) = 1;

%ETDRK4 scalars:
h = 0.1;                        % time step
k = (2.*pi./d).*[0:N/2-1 0 -N/2+1:-1]';  % wave numbers
T = 1i * [0:N/2-1 0 -N/2+1:-1]'; %U(1) generator

%Linear term:  
L = k.^2 - k.^4;

E = exp(h*L); E2 = exp(h*L/2);
M = 16;                         % no. of points for complex means
r = exp(1i*pi*((1:M)-.5)/M);    % roots of unity
LR = h*L(:,ones(M,1)) + r(ones(N,1), :);
Q  = h*real(mean(           (exp(LR/2) - 1)./LR               ,2));
f1 = h*real(mean(   (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3   ,2));
f2 = h*real(mean(       (2+LR+exp(LR).*(-2+LR))./LR.^3        ,2));
f3 = h*real(mean(   (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3    ,2));

%Nonlinear term:
g = 0.5i*k*N;
N = @(x) g.*fft(real(ifft(x)).^2);

t = 0;

%Time-stepping loop:
taut=0; tt = 0; vv = v; u=real(ifft(v)); uu = u;
tmax = tf; 

while t < tf
    
    t = t+h;
    
    Nv = N(v);
    a = E2.*v + Q.*Nv;
    Na = N(a);
    b = E2.*v + Q.*Na;
    Nb = N(b);
    c = E2.*a + Q.*(2*Nb-Nv);
    Nc = N(c);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    
    vv = [vv, v];
    
    u=real(ifft(v));
    uu = [uu, u];
   
    tt = [tt, t]; 
    
end

xsol = zeros(Nh*2-2, size(vv,2));
xsol(1:2:end-1, :) = real(vv(2:Nh,:));  
xsol(2:2:end, :) = imag(vv(2:Nh,:));
xsol = xsol';
uu = uu';
