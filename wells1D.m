function [A,q] = wells1D(pi,pw,lambda,A,q,cell)

for i = 1:length(pw)
    q(cell(i)) = q(cell(i))+ pi(i).*lambda(cell(i))*pw(i);
    A(cell(i),cell(i)) = A(cell(i),cell(i))+pi(i).*lambda(cell(i));
end

% i = 1; % left boundary
%     % T(1) *(P(1)-PL)
%     A(i,i) = A(i,i)+ T(i);
%     q(i) = q(i)+T(i)*PL;
% 
%     i = N; % left boundary
%     % T(N+1) *(P(N)-PR)
%     A(i,i) = A(i,i)+ T(i+1);
%     q(i) = q(i)+T(i+1)*PR;

end