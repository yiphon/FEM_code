function y = tetra_assembly_mat( k,kk,node_number )
%TETRA_ASSEMBLY_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成总体刚度矩阵。
a=node_number(1);
b=node_number(2);
c=node_number(3);
d=node_number(4);
dof=[3*a-2,3*a-1,3*a,3*b-2,3*b-1,3*b,3*c-2,3*c-1,3*c,3*d-2,3*d-1,3*d];
for i=1:1:12
    for j=1:1:12
        k(dof(i),dof(j))=k(dof(i),dof(j))+kk(i,j);
    end
end
y=k;
end

