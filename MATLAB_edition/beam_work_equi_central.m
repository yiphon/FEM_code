function y = beam_work_equi_central( x0,l,f )
%BEAM_WORK_EQUI_CENTRAL Summary of this function goes here
%   Detailed explanation goes here
%   �����غɵĵ�Ч�ڵ�����x0Ϊ����λ�ã�lΪ��Ԫ�ܳ���fΪ��С��
syms x;
n=[1-3*x^2/l^2+2*x^3/l^3;x-2*x^2/l+x^3/l^2;3*x^2/l^2-2*x^3/l^3;x^3/l^2-x^2/l];
x=x0;
n=f*subs(n);
y=n;
end

