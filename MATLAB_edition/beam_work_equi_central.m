function y = beam_work_equi_central( x0,l,f )
%BEAM_WORK_EQUI_CENTRAL Summary of this function goes here
%   Detailed explanation goes here
%   集中载荷的等效节点力，x0为作用位置，l为单元总长，f为大小。
syms x;
n=[1-3*x^2/l^2+2*x^3/l^3;x-2*x^2/l+x^3/l^2;3*x^2/l^2-2*x^3/l^3;x^3/l^2-x^2/l];
x=x0;
n=f*subs(n);
y=n;
end

