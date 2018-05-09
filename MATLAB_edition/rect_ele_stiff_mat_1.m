function y = rect_ele_stiff_mat_1(t_thickness,b,d,x,y )
%RECT_ELE_STIFF_MAT_1 Summary of this function goes here
%   Detailed explanation goes here
%   生成单位刚度矩阵。x\y为s\t的函数。
syms s t

j=det(jacobian([x,y],[s,t]));
kk_int0=t_thickness*j*b'*d*b;
kk_int1=int(kk_int0,s,-1,1);
kk_int=int(kk_int1,t,-1,1);
y=kk_int;

end

