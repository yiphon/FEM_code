function y = axisy_ele_stiff_mat_1( b,d,r,z )
%AXISY_ELE_STIFF_MAT_1 Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元刚度矩阵，r\z是s\t的函数。
syms s t

j=det(jacobian([r,z],[s,t]));
kk_int0=2*pi*r*j*b'*d*b;
kk_int1=int(kk_int0,s,-1,1);
kk_int=int(kk_int1,t,-1,1);
y=kk_int;


end

