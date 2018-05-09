function y = tri_ele_stiff_mat( a,t,b,d )
%TRI_ELE_STIFF_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成单元刚度矩阵，b为单元的B矩阵。
kk=a*t*b'*d*b;
y=kk;
end

