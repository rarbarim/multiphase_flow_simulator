function R = residual(q,phi,dt,dx,dy,Sw,Swn,fwn,UX,UY,NX,NY);
R = zeros(NX*NY,1);

for i = 1:NX
    for j = 1:NY
        [I,IN,IE,IS,IW] = findindex2D(j,i,NX);
        gridindex(j,i) = (j-1)*NX+i;
        
        if i > 1 && j > 1
            R(I) = (dt/phi)*q(I)-(Swn(I)-Sw(I))-dt/(dx*phi)*(fwn(I)*UX(j,i+1)-fwn(IW)*UX(j,i)) - dt/(dy*phi)*(fwn(I)*UY(j+1,i)-fwn(IN)*UY(j,i));
%         end
        
        elseif i == 1 && j ~= 1
            R(I) = (dt/phi)*q(I)-(Swn(I)-Sw(I))-dt/(dx*phi)*(fwn(I)*UX(j,i+1)-0*UX(j,i)) - dt/(dy*phi)*(fwn(I)*UY(j+1,i)-fwn(IN)*UY(j,i));
%         end
        
        elseif j == 1 && i ~= 1
            R(I) = (dt/phi)*q(I)-(Swn(I)-Sw(I))-dt/(dx*phi)*(fwn(I)*UX(j,i+1)-fwn(IW)*UX(j,i)) - dt/(dy*phi)*(fwn(I)*UY(j+1,i)-0*UY(j,i));
        end
    end
end

return

%%
% for i = 1:NX
%     for j = 1:NY
%         [I,IN,IE,IS,IW] = findindex2D(j,i,NX);
%         gridindex(j,i) = (j-1)*NX+i;
%         
%         if i > 1 && j > 1
%             R(I) = (dt/phi)*q(I)-(Swn(I)-Sw(I))-dt/(dx*phi)*(fwn(I)*UX(j,i+1)-fwn(IW)*UX(j,i)) - dt/(dy*phi)*(fwn(I)*UY(j+1,i)-fwn(IN)*UY(j,i));
%         end
%         
%         if i == 1 && j ~= 1
%             R(I) = (dt/phi)*q(I)-(Swn(I)-Sw(I))-dt/(dx*phi)*(fwn(I)*UX(j,i+1)-0*UX(j,i)) - dt/(dy*phi)*(fwn(I)*UY(j+1,i)-fwn(IN)*UY(j,i));
%         end
%         
%         if j == 1 && i ~= 1
%             R(I) = (dt/phi)*q(I)-(Swn(I)-Sw(I))-dt/(dx*phi)*(fwn(I)*UX(j,i+1)-fwn(IW)*UX(j,i)) - dt/(dy*phi)*(fwn(I)*UY(j+1,i)-0*UY(j,i));
%         end
%     end
% end
% 
% return