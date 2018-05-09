function y = truss3d_ele_mat_assemble( k,kk,node_number )
%ELE_MAT_ASSEMBLE_3DTRUSS Summary of this function goes here
%   Detailed explanation goes here
%   完成单元刚度矩阵，i、j为单元起点终点的节点编号。
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

