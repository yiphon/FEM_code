function y = frame3d_max_bending_sigma( e,l,u )
%FRAME3D_MAX_BENDING_SIGMA Summary of this function goes here
%   Detailed explanation goes here
%   求单元最大弯曲应力值。
syms x y z
n=[1-3*x^2/l^2+2*x^3/l^3;x-2*x^2/l+x^3/l^2;3*x^2/l^2-2*x^3/l^3;x^3/l^2-x^2/l];
n=diff(n,x);
uy=u([3,5,9,11]);   %绕y轴转需要考虑的位移。
uz=u([2,6,8,12]);
func=-e*n'*uy-e*n'*uz
func=matlabFunction(func)   %此处为方便求算取了负值，其最小值为最大值。
 x=0:0.01:l;
 plot(x,func(x));
 hold on
[x0,ymax]=fminbnd(func,0,l);
ymax=-ymax;
y_left=-func(0);
y_right=-func(l);
candidate_y=[y_left,ymax,y_right]
candidate_x=[0,x0,l]
[y_max,index]=max(candidate_y);
y=[y_max;candidate_x(index)];
end

