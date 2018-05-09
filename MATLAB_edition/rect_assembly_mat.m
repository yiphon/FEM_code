function y = rect_assembly_mat( K,kk,node_number )
%RECT_ASSEMBLY_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成整体刚度矩阵，node_number=[i,j,k,l]'。
i=node_number(1);
j=node_number(2);
k=node_number(3);
l=node_number(4);
dof=[2*i-1,2*i,2*j-1,2*j,2*k-1,2*k,2*l-1,2*l];
for m=1:1:8
    for n=1:1:8
        K(dof(m),dof(n))=K(dof(m),dof(n))+kk(m,n);
    end
end
y=K;
end

