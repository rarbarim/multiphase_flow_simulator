function [fw,fo,dfwds] = fractionalflow(lambdaw,lambdao,dlambdaw,dlambdao);
       
    fw = (lambdaw)./(lambdaw+lambdao);
    fo = 1-fw;
    
    dfwds = ((dlambdaw.*(lambdao+lambdaw))-((dlambdao+dlambdaw).*lambdaw))./((lambdaw+lambdao).^2);
    
return