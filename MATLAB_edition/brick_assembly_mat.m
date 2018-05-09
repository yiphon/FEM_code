function y = brick_assembly_mat( k,kk,node_number )
%BRICK_ASSEMBLY_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成总体刚度矩阵。
a=node_number(1);
b=node_number(2);
c=node_number(3);
d=node_number(4);
e=node_number(5);
f=node_number(6);
g=node_number(7);
h=node_number(8);
dof=[3*a-2,3*a-1,3*a,3*b-2,3*b-1,3*b,3*c-2,3*c-1,3*c,3*d-2,3*d-1,3*d,3*e-2,3*e-1,3*e,3*f-2,3*f-1,3*f,3*g-2,3*g-1,3*g,3*h-2,3*h-1,3*h];
for i=1:1:24
    for j=1:1:24
        k(dof(i),dof(j))=k(dof(i),dof(j))+kk(i,j);
    end
end
y=k;
end

