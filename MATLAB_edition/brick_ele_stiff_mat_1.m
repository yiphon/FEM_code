function y = brick_ele_stiff_mat_1( b,d,x,y,z )
%BRICK_ELE_STIFF_MAT_1 Summary of this function goes here
%   Detailed explanation goes here
%   生成单位刚度矩阵。x\y为s\t\zi的函数。
syms s t zi

j=det(jacobian([x,y,z],[s,t,zi]));
kk_int0=j*b'*d*b;
kk_int1=int(kk_int0,s,-1,1);
kk_int2=int(kk_int1,t,-1,1);
kk_int=int(kk_int2,zi,-1,1);
y=kk_int;

end

