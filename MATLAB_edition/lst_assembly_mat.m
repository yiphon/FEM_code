function y = lst_assembly_mat( k,kk,node_number )
%LST_ASSEMBLY_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成总体刚度矩阵。
l=node_number(1);
m=node_number(2);
n=node_number(3);
o=node_number(4);
p=node_number(5);
q=node_number(6);
dof=[2*l-1,2*l,2*m-1,2*m,2*n-1,2*n,2*o-1,2*o,2*p-1,2*p,2*q-1,2*q];
for i=1:1:12
    for j=1:1:12
        k(dof(i),dof(j))=k(dof(i),dof(j))+kk(i,j);
    end
end
y=k;


end

