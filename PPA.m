function [y,min,i]=PPA(A,b,x0,mu,iteration,tolerance)
[n,~]=size(x0);
syms('x',[n,1]);
f=0.5*(norm(A*x-b))^2;
g=mu*norm(x,1);
grad_f=A'*(A*x-b);
[~,D]=eig(A'*A);
eigenvalue=diag(D);
lamda_max=max(eigenvalue); %l取A'A的最大特征值即可满足L-smooth条件
t=1/(lamda_max+36); %此处取t<1/L即可
y=x0;
v=x0;
i=0;
flag=1;
while flag
    f_temp=subs(f,x,y)+subs(g,x,y);
    f_temp=double(f_temp);
    gamma=2/(i+2);
    z=(1-gamma)*y+gamma*v;
    y_temp=y;
    grad=double(subs(grad_f,x,z));
    y=prox(z-t*grad,t*mu);
    v=y_temp+(y-y_temp)/gamma;
    i=i+1;
    f_value=double(subs(f,x,y));
    %使用近似点并用FISTA加速
    if i<=iteration&&abs(f_value-f_temp)>tolerance&&norm(y-y_temp)>tolerance
        flag=1;
    else
        flag=0;
    end %判断是否需要结束循环
end
min=double(subs(g,x,y));
min=min+f_value;
end

function y = prox(x, mu)%此为求mu|x|_1的proximal point
y = max(abs(x) - mu, 0);
y = sign(x) .* y;
end