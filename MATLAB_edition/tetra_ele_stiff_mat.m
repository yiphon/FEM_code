function y = tetra_ele_stiff_mat( b,d,v )
%TETRA_ELE_STIFF_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元刚度矩阵。
y=v*b'*d*b;

end

