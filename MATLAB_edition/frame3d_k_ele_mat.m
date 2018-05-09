function y = frame3d_k_ele_mat( ele_info )
%FRAME3D_K_ELE_MAT Summary of this function goes here
%   Detailed explanation goes here
%   完成局部坐标系下的单位刚度矩阵，ele_info=[e;a;iy;iz;j;l;;;;g]。
e=ele_info(1);
a=ele_info(2);
iy=ele_info(3);
iz=ele_info(4);
j=ele_info(5);
l=ele_info(6);
g=ele_info(10);
x1=e*a/l;
x2=12*e*iz/l^3;
x3=6*e*iz/l^2;
x4=12*e*iy/l^3;
x5=6*e*iy/l^2;
x6=g*j/l;
x7=4*e*iy/l;
x8=2*e*iy/l;
x9=4*e*iz/l;
x10=2*e*iz/l;
k=[x1,0,0,0,0,0,-x1,0,0,0,0,0;
    0,x2,0,0,0,x3,0,-x2,0,0,0,x3;
    0,0,x4,0,-x5,0,0,0,-x4,0,-x5,0;
    0,0,0,x6,0,0,0,0,0,-x6,0,0;
    0,0,-x5,0,x7,0,0,0,x5,0,x8,0;
    0,x3,0,0,0,x9,0,-x3,0,0,0,x10;
    -x1,0,0,0,0,0,x1,0,0,0,0,0;
    0,-x2,0,0,0,-x3,0,x2,0,0,0,-x3;
    0,0,-x4,0,x5,0,0,0,x4,0,x5,0;
    0,0,0,-x6,0,0,0,0,0,x6,0,0;
    0,0,-x5,0,x8,0,0,0,x5,0,x7,0;
    0,x3,0,0,0,x10,0,-x3,0,0,0,x9];
y=k;

end

