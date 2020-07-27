function [UX,UY] = computevelocity2D(lambdaharx,lambdahary,dx,dy,P,NX,NY,i,j)
UX = zeros(NY,NX+1);
UY = zeros(NY+1,NX);
[I,IN,IE,IS,IW] = findindex2D(j,i,NX);
% P = reshape(P,NY,NX);
P = vec2mat(P,NX);

for m = 1:NY
    for n = 2:NX
    UX(m,n)=-lambdaharx(m,n)*(P(m,n)-P(m,n-1))/dx;
    end
end

for m = 2:NY
    for n = 1:NX
    UY(m,n)=-lambdahary(m,n)*(P(m,n)-P(m-1,n))/dy;
    end
end

% for m = 2:NY+1
%     for n = 2:NX+1
%     UX(m,n)=-lambdaharx(m,n)*(P(m,n)-P(m,n-1))/dx;
%     UY(m,n)=-lambdahary(m,n)*(P(m,n)-P(m-1,n))/dy;
% %     U(i)=-lambdahar(i)*(P(i)-P(i-1))/dx;
%     end
% end

end