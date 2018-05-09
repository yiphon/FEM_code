function y = frame3d_assembly_mat( k,kk,node_number )
%FRAME3D_ASSEMBLY_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成整体刚度矩阵，kk为对应单元的刚度矩阵，node_number=[i;j]。
i=node_number(1);
j=node_number(2);
dof=[6*i-5;6*i-4;6*i-3;6*i-2;6*i-1;6*i;6*j-5;6*j-4;6*j-3;6*j-2;6*j-1;6*j];
for m=1:1:12
    for n=1:1:12
        k(dof(m),dof(n))=k(dof(m),dof(n))+kk(m,n);
    end
end
y=k;
end


