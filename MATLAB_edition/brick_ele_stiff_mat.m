function aa = brick_ele_stiff_mat( b,d,x,y,z )
%BRICK_ELE_STIFF_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元刚度矩阵。
syms s t zi
j=det(jacobian([x,y,z]));
kk=b'*d*b;
kk_ints=int(j*kk,s,-1,1);
kk_intst=int(kk_ints,t,-1,1);
kk_intstzi=int(kk_intst,zi,-1,1);
aa=kk_intstzi;
end

