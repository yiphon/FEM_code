function fun = lst_ele_stiff_mat( n_c,t,b,d )
%LST_ELE_STIFF_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成单元刚度矩阵，此时非常值积分，需要描述三角形单元的范围区间。
syms x y
x1=n_c(1,1);
y1=n_c(2,1);
x2=n_c(1,2);
y2=n_c(2,2);
x3=n_c(1,3);
y3=n_c(2,3);
l1=y1+(y2-y1)*(x-x1)/(x2-x1);
l2=y1+(y3-y1)*(x-x1)/(x3-x1);
l3=y2+(y3-y2)*(x-x2)/(x3-x2);
kk=t*b'*d*b;
kk1_int1=int(kk,y,l1,l2);
kk1_int2=int(kk1_int1,x,x1,x3);
kk2_int1=int(kk,y,l1,l3);
kk2_int2=int(kk2_int1,x,x3,x2);
fun=double(kk1_int2+kk2_int2);

end

