function y = beam_work_equi_distributed( l,q )
%BEAM_WORK_EQUIVALENCE Summary of this function goes here
%   Detailed explanation goes here
%   �ֲ����غ��ڽڵ�����ĵ�Ч��,�˴����⿼�Ǿ��Ⱥ���ֲ�������СΪq��
syms x;
n=[1-3*x^2/l^2+2*x^3/l^3;x-2*x^2/l+x^3/l^2;3*x^2/l^2-2*x^3/l^3;x^3/l^2-x^2/l];
f=zeros(4,1);
for i=1:1:4
    f(i)=int(n(i)*q,0,l);
end
y=f;
end

