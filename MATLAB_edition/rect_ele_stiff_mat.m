function y = rect_ele_stiff_mat( t,length_info,b,d )
%RECT_ELE_STIFF_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成单位刚度矩阵。
syms x y
x0=length_info(1);
y0=length_info(2);
kk_int0=t*b'*d*b;
kk_int1=int(kk_int0,x,-x0,x0);
kk_int=int(kk_int1,y,-y0,y0);
y=kk_int;
end

