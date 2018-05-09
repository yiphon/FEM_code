function y = rect_ele_b_1( Ni,x,y )
%RECT_ELE_B_1 Summary of this function goes here
%   Detailed explanation goes here
%   ���ɵ�Ԫ��B����Ni��x��y��Ϊs\t�ĺ�����
syms s t 
b=sym(zeros(3,2));
j=det(jacobian([x,y],[s,t]));
ni_x=det(jacobian([Ni,y],[s,t]))/j;
ni_y=det(jacobian([x,Ni],[s,t]))/j;
b=[ni_x,0;0,ni_y;ni_y,ni_x];
y=b;

end

