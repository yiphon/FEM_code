function y = frame2d_assembly_mat( k,kk,node_number )
%FRAME2D_ASSEMBLY_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成整体刚度矩阵，kk为对应单元的刚度矩阵，node_number=[i;j]。
i=node_number(1);
j=node_number(2);
dof=[3*i-2;3*i-1;3*i;3*j-2;3*j-1;3*j];
for m=1:1:6
    for n=1:1:6
        k(dof(m),dof(n))=k(dof(m),dof(n))+kk(m,n);
    end
end
y=k;
end

