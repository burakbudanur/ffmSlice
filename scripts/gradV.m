function A = gradV(x)
%Nazmi Burak Budanur
%Kuramoto-Sivashinsky stability matrix function
%based on Xiong Ding's stability.m
d = 22;
N = 32;
Nh = N/2;

x=[zeros(1,size(x,2));x(1:2:end-1,:)+1i*x(2:2:end,:); ...
   zeros(1,size(x,2));x(end-1:-2:1,:)-1i*x(end:-2:2,:)];

tmpa=real(ifft(x));

w=zeros(N,N-2); w(2:2*N+1:end) = 1;  w(N+2:2*N+1:end) = 1i;

v=zeros(N,N-2); v(2:2*N+1:end) = 1;  v(N+2:2*N+1:end) = 1i;
v(N:2*N-1:end) = 1;  v(2*N:2*N-1:end) = -1i; 

k = (2.*pi./d).*[0:Nh-1 0 -Nh+1:-1]'; L =repmat( k.^2 - k.^4,1,N-2);
g = repmat(1i*k*N,1,N-2);

A=[];
Ac=L.*w+g.*fft(repmat(tmpa,1,N-2).*real(ifft(v)));

for ii=1:N-2
	tA=Ac(2:Nh,ii); tA=[real(tA), imag(tA)]'; A=[A,tA(:)];
	%tA(:) will rearrange tA collomn wise.
end
