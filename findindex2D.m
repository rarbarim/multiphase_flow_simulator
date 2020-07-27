function [I,IN,IE,IS,IW] = findindex2D(j,i,NX)
I = (j-1)*NX+i;
IE = (j-1)*NX+i+1;
IW = (j-1)*NX+i-1;
IS = (j)*NX+i;
IN = (j-2)*NX+i;

return