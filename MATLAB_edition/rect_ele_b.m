function y = rect_ele_b( jacob )
%RECT_ELE_B Summary of this function goes here
%   Detailed explanation goes here
%   ���ɵ�Ԫ��B����jacobΪ�õ�Ԫ��ֵ�������ſɱȾ���
b=sym(zeros(3,8));
for i=1:1:4
    b(:,2*i-1:2*i)=[jacob(i,1),0;0,jacob(i,2);jacob(i,2),jacob(i,1)];
end
y=b;

end

