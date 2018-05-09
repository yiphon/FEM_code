function fun = lst_stress_center( n_c,stress )
%LST_STRESS_CENTER Summary of this function goes here
%   Detailed explanation goes here
%   求各单元中心位置的应力值。
syms x y
x1=n_c(1,1);
y1=n_c(2,1);
x2=n_c(1,2);
y2=n_c(2,2);
x3=n_c(1,3);
y3=n_c(2,3);
x_center=(x1+x2+x3)/3;
y_center=(y1+y2+y3)/3;
fun=double(subs(stress,{x,y},{x_center,y_center}));

end

