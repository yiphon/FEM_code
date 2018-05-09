function y = beam_strain( E,y0,x0,u,l)
%BEAM_STRAIN Summary of this function goes here
%   Detailed explanation goes here
%   计算某位置的应力应变。x0、y0为坐标，u是该单元的节点位移矩阵，四行一列，l为单元长度。
syms x;
n=[1-3*x^2/l^2+2*x^3/l^3;x-2*x^2/l+x^3/l^2;3*x^2/l^2-2*x^3/l^3;x^3/l^2-x^2/l]';
v=n*u;
strain=-y0*diff(v,2);
x=x0;
y=double([subs(strain);E*subs(strain)]);
end

