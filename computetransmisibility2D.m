function [TX,TY,lambdaharx,lambdahary] = computetransmisibility2D(lambda,dx,NX,dy,NY)

lambdaharx(:,1) = lambda(:,1);
lambdaharx(:,2:NX) = 2*lambda(:,2:NX).*lambda(:,1:NX-1)./(lambda(:,2:NX)+lambda(:,1:NX-1));
lambdaharx(:,NX+1) = lambda(:,NX);

lambdahary(1,:) = lambda(1,:);
lambdahary(2:NY,:) = 2*lambda(2:NY,:).*lambda(1:NY-1,:)./(lambda(2:NY,:)+lambda(1:NY-1,:));
lambdahary(NY+1,:) = lambda(NY,:);

TX(:,1)=lambdaharx(:,1)./((dx^2)/2); 
TX(:,2:NX)=lambdaharx(:,2:NX)./(dx^2);
TX(:,NX+1)=lambdaharx(:,NX+1)./((dx^2)/2);

TY(1,:)=lambdahary(1,:)./((dy^2)/2);
TY(2:NY,:)=lambdahary(2:NY,:)./(dy^2);
TY(NY+1,:)=lambdahary(NY+1,:)./((dy^2)/2);

% figure(1)
% hold on
% plot(x,lambda,'linewidth',1)
% plot(xi,lambdahar,'linewidth',1)
% xlabel('x location')
% ylabel('lambda / lambdahar')
% legend('lambda','lambdahar')
% title('lambda in incompressible flow')
% hold off
% 
% figure
% plot(xi,T)
% xlabel('xi')
% ylabel('T')
% figure
% plot(xi,lambdahar)
% xlabel('x')
% ylabel('lambdahar')

end