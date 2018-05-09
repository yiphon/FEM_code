function fun = quad2_assembly_mat( k,kk,node_number )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   生成整体刚度矩阵。
a=node_number(1);
b=node_number(2);
c=node_number(3);
d=node_number(4);
e=node_number(5);
f=node_number(6);
g=node_number(7);
h=node_number(8);
dof=[2*a-1,2*a,2*b-1,2*b,2*c-1,2*c,2*d-1,2*d,2*e-1,2*e,2*f-1,2*f,2*g-1,2*g,2*h-1,2*h];
for m=1:1:16
    for n=1:1:16
        k(dof(m),dof(n))=k(dof(m),dof(n))+kk(m,n);
    end
end
fun=k;


end

