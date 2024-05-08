function [y,min,k]=TrustRegion(iteration,tolerance,x0,f)
[~,n]=size(x0);
syms('x',[1,n]);
h=@(x)(f);
grad_f=gradient(h,x);
hessian_f=jacobian(grad_f,x);
B=subs(hessian_f,x,x0);
B=double(B);
g=double(subs(grad_f,x,x0));
y=x0;
p1=0.25;
p2=0.75;
gamma1=0.25;
gamma2=2; %这四排为参数设置 可根据具体情况更改
Delta=3; %给定初始半径
Delta_max=6; %给定最大半径
for i=1:iteration
    if norm(g)<tolerance
        break;
    end
    %使用随机方法找到下降方向
    flag=1;
    while flag
        p=rand(n,1); %生成（0，1）^n的随机向量
        p=2*p-ones(n,1); %变为（-1，1）^n的随机向量
        if p'*g<0
            if p'*B*p<=0
                p=(Delta/norm(p))*p;
            else
                if -g'*p*norm(p)/(p'*B*p)<=Delta
                    p=-g'*p*norm(p)/(p'*B*p)*p;
                else
                    p=(Delta/norm(p))*p;
                end
            end
            flag=0;
        end
    end
    rho=subs(f,x,y+p')-subs(f,x,y); %下为调试半径
    rho=double(rho/(g'*p+0.5*p'*B*p));
    if rho>p2
        if gamma2*Delta>Delta_max
            Delta=Delta_max;
        else
            Delta=gamma2*Delta;
        end
    else
        Delta=gamma1*Delta;
    end
    y=y+p';
    g=double(subs(grad_f,x,y));
    B=double(subs(hessian_f,x,y));
end
k=i;
min=double(subs(f,x,y));
end


