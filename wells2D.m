function [A,q] = wells2D(pi,pw,lambda,A,q,grid)

for i = 1:length(pw)
    q(grid(i)) = q(grid(i))+ pi(i).*lambda(grid(i))*pw(i);
    A(grid(i),grid(i)) = A(grid(i),grid(i))+pi(i).*lambda(grid(i));
end

end