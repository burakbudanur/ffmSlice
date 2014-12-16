function vx = vel(x)

d = 22;
N = 32;
Nh = N/2;

T = 1i * [0:N/2-1 0 -N/2+1:-1]'; %U(1) generator

v = [0; x(1:2:end-1)+1i*x(2:2:end); 0; x(end-1:-2:1)-1i*x(end:-2:2)];

k = (2.*pi./d).*[0:N/2-1 0 -N/2+1:-1]';  % wave numbers

L = k.^2 - k.^4;

g = 0.5i*k*N;

N = @(x) g.*fft(real(ifft(x)).^2); 
vc = @(x) L.*x + N(x);

vcomplex = vc(v);

vx = zeros(Nh*2-2,1);

vx(1:2:end-1) = real(vcomplex(2:Nh));  
vx(2:2:end, :) = imag(vcomplex(2:Nh));
