function vel = vksred(x)

d = 22;
N = 32;
Nh = N/2;

T = 1i * [0:N/2-1 0 -N/2+1:-1]'; %U(1) generator

v = [0; x(1:2:end-1)+1i*x(2:2:end); 0; x(end-1:-2:1)-1i*x(end:-2:2)];

k = (2.*pi./d).*[0:N/2-1 0 -N/2+1:-1]';  % wave numbers

Lold = k.^2 - k.^4;
L = Lold;
g = 0.5i*k*N;


Nold = @(x) g.*fft(real(ifft(x)).^2); 
vold = @(x) Lold.*x + Nold(x);

pick2ndel = zeros(1, N);
pick2ndel(2) = 1;

%N = @(x) Nold(x) + (real(x(2))-1)*vold(x) - imag(vold(x)(2))*(T.*x);
N = @(x) Nold(x) + (real(x(2))-1)*vold(x) - imag(dot(pick2ndel, vold(x)))*(T.*x);

vred = @(x) L.*x + N(x);

velv = vred(v);

vel = zeros(Nh*2-2,1);

vel(1:2:end-1) = real(velv(2:Nh));  
vel(2:2:end, :) = imag(velv(2:Nh));
