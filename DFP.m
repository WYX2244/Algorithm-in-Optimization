function [y,min,k]=DFP(iteration,tolerance,x0,f)
[~,n]=size(x0);
syms('x',[1,n]);
g=@(x)(f);
grad_f=gradient(g,x);
hess=jacobian(grad_f,x);
B=double(subs(hess,x,x0));
H=inv(B);
y=x0;
for i=1:iteration
    grad=double(subs(grad_f,x,y));
    d=-H*grad;
    y=double(y+d');
    temp=double(subs(grad_f,x,y))-grad;
    H=H-(H*temp*temp'*H')/(temp'*H*temp)+(d*d')/(temp'*d);
    %以下为另一种H的迭代算法：pho=1/(temp'*d);
    %H=(eye(n)-pho*temp*d')'*H*(eye(n)-pho*temp*d')+pho*d*d';
    if norm(grad)<tolerance
        break;
    end
end
k=i;
min=double(subs(f,x,y));
end
