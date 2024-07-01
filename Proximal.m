%本程序minimize 0.5||Ax-b||_2^2+mu||x||_1
function [y,min,i]=Proximal(A,b,x0,mu,iteration,tolerance)
[n,~]=size(x0);
syms('x',[n,1]);
f=0.5*(norm(A*x-b))^2;
g=mu*norm(x,1);
grad_f=A'*(A*x-b);
% 求解L—连续的L数值
[~,D]=eig(A'*A);
eigenvalue=diag(D);
lamda_max=max(eigenvalue); %l取A'A的最大特征值即可满足L-smooth条件
t=1/(lamda_max+36); %此处取t<1/L即可
i=0;
y=x0;
z=zeros(size(x0));
while i<=iteration && norm(y-z)>tolerance %此处使用梯度变化来判断
    temp=double(subs(grad_f,x,y));
    temp=y-t*temp;
    z=y;
    y=sign(temp).*max(abs(temp)-t*mu,0);
    i=i+1;
end
min=subs(f,x,y)+subs(g,x,y);
min=double(min);
end
