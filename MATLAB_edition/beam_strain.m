function y = beam_strain( E,y0,x0,u,l)
%BEAM_STRAIN Summary of this function goes here
%   Detailed explanation goes here
%   ����ĳλ�õ�Ӧ��Ӧ�䡣x0��y0Ϊ���꣬u�Ǹõ�Ԫ�Ľڵ�λ�ƾ�������һ�У�lΪ��Ԫ���ȡ�
syms x;
n=[1-3*x^2/l^2+2*x^3/l^3;x-2*x^2/l+x^3/l^2;3*x^2/l^2-2*x^3/l^3;x^3/l^2-x^2/l]';
v=n*u;
strain=-y0*diff(v,2);
x=x0;
y=double([subs(strain);E*subs(strain)]);
end

