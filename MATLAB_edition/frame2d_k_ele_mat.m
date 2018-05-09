function y = frame2d_k_ele_mat( ele_info )
%FRAME2D_K_ELE_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成局部坐标系下的单元刚度矩阵，ele_info形式为[e;a;iz;l]。
e=ele_info(1);
a=ele_info(2);
iz=ele_info(3);
l=ele_info(4);
x1=e*a/l;
x2=12*e*iz/l^3;
x3=6*e*iz/l^2;
x4=4*e*iz/l;
x5=2*e*iz/l;
aa=[x1,0,0,-x1,0,0;0,x2,x3,0,-x2,x3;0,x3,x4,0,-x3,x5;-x1,0,0,x1,0,0;0,-x2,-x3,0,x2,-x3;0,x3,x5,0,-x3,x4];
y=aa;
end

