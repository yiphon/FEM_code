function y = truss2d_ele_mat_assembly( k,kk,node_number )
%TRUSS2D_ELE_MAT_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here
%   生成总体刚度矩阵。
a=node_number(1);
b=node_number(2);
dof=[2*a-1,2*a,2*b-1,2*b];
for i=1:1:4
    for j=1:1:4
        k(dof(i),dof(j))=k(dof(i),dof(j))+kk(i,j);
    end
end
y=k;
end

