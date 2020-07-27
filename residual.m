function R = residual(q,phi,dt,dx,Sw,Swn,fwn,U,N);
R = zeros(N,1);
for i = 2:N
    R(i) = (dt/phi)*q(i)-(Swn(i)-Sw(i))-dt/(dx*phi)*(fwn(i)*U(i+1)-fwn(i-1)*U(i));
end

return 