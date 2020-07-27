function [U] = computevelocity1D(lambdahar,dx,P,N)
U = zeros(1,N+1);

for i = 2:N
    U(i)=-lambdahar(i)*(P(i)-P(i-1))/dx;
end

end