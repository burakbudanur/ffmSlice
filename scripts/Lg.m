function T = Lg()

N = 16;

T = zeros(2*N-2,2*N-2);
    
for i=1:1:N-1

	T(2*i-1, 2*i) = -i;
	T(2*i, 2*i-1) = i;

end
