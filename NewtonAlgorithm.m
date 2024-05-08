function [y,min,i]=NewtonAlgorithm(iteration,tolerance,x0)
[~,n]=size(x0);
syms('x',[1,n]);
f(x(1),x(2))=x(1)^2+3*x(2)^3-3*x(1)*x(2);
grad_f(x(1),x(2))=gradient(f,x);
hessian_f(x(1),x(2))=jacobian(grad_f,x);
history_x=zeros(2,iteration);
history_f=zeros(1,iteration);
for i=1:iteration
    grad=double(grad_f(x0(1),x0(2)));
    hess=double(hessian_f(x0(1),x0(2)));
    d=-inv(hess)*grad;
    x0=double(x0+d');
    history_x(1,i)=x0(1);
    history_x(2,i)=x0(2);
    temp=f(x0(1),x0(2));
    history_f(1,i)=temp;
    if norm(grad)<tolerance
        break;
    end
end
y=x0;
min=double(f(y(1),y(2)));

% 以下为可视化
figure;
scatter3(history_x(1,1:i),history_x(2,1:i),history_f(1,1:i),'o');
hold on;
for j=1:i
  c=num2str(j);
  text(history_x(1,j),history_x(2,j),history_f(1,j),c);
end
title('迭代数值');
end
