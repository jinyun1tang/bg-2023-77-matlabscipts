function hmat=hessianmat(fun,x00)
%finite difference approximation of the hessian matrix of 
%the multivariate function fun at x
%Jinyun Tang (jinyuntang@gmail.com)
%History: Mar 5th, 2023
nz=numel(x00);
hmat=zeros(nz);
x0=x00+eps*randn(size(x00));
for ii = 1 : nz
    for jj = 1 : nz        
        x1=x0;x2=x0;x3=x0;x4=x0;
        x1(ii)=x1(ii)+0.01*x0(ii);x1(jj)=x1(jj)+0.01*x0(jj);
        x2(ii)=x2(ii)+0.01*x0(ii);x2(jj)=x2(jj)-0.01*x0(jj);
        x3(ii)=x3(ii)-0.01*x0(ii);x3(jj)=x3(jj)+0.01*x0(ii);
        x4(ii)=x4(ii)-0.01*x0(ii);x4(jj)=x4(jj)-0.01*x0(jj);
        f1=feval(fun,x1);
        f2=feval(fun,x2);
        f3=feval(fun,x3);
        f4=feval(fun,x4);
        hmat(ii,jj)=(f1-f2-f3+f4)./(4*0.01^2*x0(ii)*x0(jj));
    end
end

end