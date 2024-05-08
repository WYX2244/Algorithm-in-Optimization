function [z,min,k]=GradDe(x0,f,e)
[m,n]=size(x0);
syms('x',[1 n]);
if ischar(f)
    f=sym(f);
end
grad=gradient(sym(f),x);
d=-grad';
flag=1;
k=1;
z=x0;
while(flag)
    d_temp=subs(d,x,z);
    nor=norm(d_temp);
    if(nor>=e)
        z=z+0.1*d_temp;
        k=k+1;
    else
        flag=0;
    end
end
k=k-1;
z=vpa(z);
min=vpa(subs(f,x,z));
end
