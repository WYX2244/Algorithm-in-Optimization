function [y,min,k]=BFGS(iteration,tolerance,x0,f)
[~,n]=size(x0);
syms('x',[1,n]);
g=@(x)(f);
grad_f=gradient(g,x);
hess=jacobian(grad_f,x);
B=double(subs(hess,x,x0));
y=x0;
for i=1:iteration
    grad=double(subs(grad_f,x,y));
    d=-inv(B)*grad;
    y=double(y+d');
    temp=double(subs(grad_f,x,y))-grad;
    B=B+(temp*temp')/(d'*temp)-(B*d*d'*B')/(d'*B*d);
    if norm(grad)<tolerance
        break;
    end
end
k=i;
min=double(subs(f,x,y));
end
