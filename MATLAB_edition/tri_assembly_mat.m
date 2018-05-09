function y = tri_assembly_mat( k,kk,node_number)
%TRI_ASSEMBLY_MAT Summary of this function goes here
%   Detailed explanation goes here
%   完成总体刚度矩阵，kk为单元刚度矩阵，node_number=[l;m;n]。
l=node_number(1);
m=node_number(2);
n=node_number(3);
dof=[2*l-1,2*l,2*m-1,2*m,2*n-1,2*n];
for i=1:1:6
    for j=1:1:6
        k(dof(i),dof(j))=k(dof(i),dof(j))+kk(i,j);
    end
end
y=k;
end


