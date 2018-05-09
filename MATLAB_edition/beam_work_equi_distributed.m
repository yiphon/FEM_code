function y = beam_work_equi_distributed( l,q )
%BEAM_WORK_EQUIVALENCE Summary of this function goes here
%   Detailed explanation goes here
%   分布力载荷在节点产生的等效力,此处特殊考虑均匀横向分布力，大小为q。
syms x;
n=[1-3*x^2/l^2+2*x^3/l^3;x-2*x^2/l+x^3/l^2;3*x^2/l^2-2*x^3/l^3;x^3/l^2-x^2/l];
f=zeros(4,1);
for i=1:1:4
    f(i)=int(n(i)*q,0,l);
end
y=f;
end

