function y = frame2d_ele_stiff_mat( r,k_ele )
%FRAME2D_ELE_STIFF_MAT Summary of this function goes here
%   Detailed explanation goes here
%   完成整体坐标系下单元刚度矩阵，r与k是对应单元的转移矩阵和局部坐标单元刚度矩阵。
a=r'*k_ele*r;
y=a;
end

