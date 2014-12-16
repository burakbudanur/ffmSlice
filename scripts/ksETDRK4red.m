function [xsol, uu, tt, taut] = ksETDRK4red(x0, tf)
%
%Adaptation of kursiv.m from Kassam ^ Trefethen (2005)
%
%Outputs:
%xsol = solution in the state space
%uu = solution in the configuration space
%tt = lab time of the solution
%taut = slice time of the solution
%
%Inputs:
%x0 = initial condition (must be in the first mode slice)
%tf = simulation duration

%Spatial grid:
d = 22;
N = 32;
Nh = N/2;
x = d*(1:N)'/N;

v = [0; x0(1:2:end-1)+1i*x0(2:2:end); 0; x0(end-1:-2:1)-1i*x0(end:-2:2)];

%Template:
slicep = zeros(N,1);
slicep(2) = 1;
slicep(N) = 1;

%ETDRK4 scalars:
h0 = 0.1;                        % time step
k = (2.*pi./d).*[0:N/2-1 0 -N/2+1:-1]';  % wave numbers
T = 1i * [0:N/2-1 0 -N/2+1:-1]'; %U(1) generator
%Linear term:  
Lold = k.^2 - k.^4;
L = Lold;

E0 = exp(h0*L); E20 = exp(h0*L/2);
M = 16;                         % no. of points for complex means
r = exp(1i*pi*((1:M)-.5)/M);    % roots of unity
LR0 = h0*L(:,ones(M,1)) + r(ones(N,1), :);
Q0  = h0*real(mean(           (exp(LR0/2) - 1)./LR0                 ,2));
f10 = h0*real(mean(   (-4-LR0+exp(LR0).*(4-3*LR0+LR0.^2))./LR0.^3   ,2));
f20 = h0*real(mean(       (2+LR0+exp(LR0).*(-2+LR0))./LR0.^3        ,2));
f30 = h0*real(mean(   (-4-3*LR0-LR0.^2+exp(LR0).*(4-LR0))./LR0.^3   ,2));

h = h0;
E = E0; E2 = E20;

LR = LR0;
Q  = Q0;
f1 = f10;
f2 = f20;
f3 = f30;

%Nonlinear term:
g = 0.5i*k*N;
Nold = @(x) g.*fft(real(ifft(x)).^2); fft(real(ifft(v)).^2);
vold = @(x) Lold.*x + Nold(x);

%Awful workaround for Matab compatibility:
pick2ndel = zeros(N, 1);
pick2ndel(2) = 1;
vold2 = @(x) (dot(pick2ndel, (Lold.*x + Nold(x))));

%Reduced velocity:
tp = T.*slicep;
vhat = @(x) vold(x) - ((vold(x)'*tp)/((T.*x)'*tp)) * (T.*x);

N = @(x) Nold(x) + (real(x(2))-1)*vold(x) - imag(vold2(x))*(T.*x);

vnew = @(x) L.*x + N(x);

t = 0;
tau = 0;

%Time-stepping loop:
taut=0; tt = 0; vv = v; u=real(ifft(v)); uu = u;
tmax = tf; 

while t < tf
    
    tau = tau+h;
    
    Nv = N(v);
    a = E2.*v + Q.*Nv;
    Na = N(a);
    b = E2.*v + Q.*Na;
    Nb = N(b);
    c = E2.*a + Q.*(2*Nb-Nv);
    Nc = N(c);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    
    taut = [taut tau];
    vv = [vv, v];

    u=real(ifft(v));
    uu = [uu, u];
    
    t = t + trapz([0 h], real(vv(2,end-1:end)));
    tt = [tt, t]; 
    
end

xsol = zeros(Nh*2-2, size(vv,2));
xsol(1:2:end-1, :) = real(vv(2:Nh,:));  
xsol(2:2:end, :) = imag(vv(2:Nh,:));
xsol = xsol';
uu = uu';
