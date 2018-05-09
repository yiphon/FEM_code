function y = frame2d_max_sigma( e,l,u,f,a)
%FRAME2D_MAX_SIGMA Summary of this function goes here
%   Detailed explanation goes here
%   ���ɵ�Ԫ����Ӧ��ֵ��x���ꡣuΪ��Ԫ�������ڵ�λ������
syms x
n=[1-3*x^2/l^2+2*x^3/l^3;x-2*x^2/l+x^3/l^2;3*x^2/l^2-2*x^3/l^3;x^3/l^2-x^2/l];
n=diff(n,x);
u=u([2,3,5,6]);
func=-e*n'*u;
func=matlabFunction(func);   %�˴�Ϊ��������ȡ�˸�ֵ������СֵΪ���ֵ��
f_axial=f(4);
 x=0:0.01:45;
 plot(x,func(x));
 hold on
s1=f_axial/a;
[x0,ymax]=fminbnd(func,0,l);
ymax=-ymax+s1;
y_left=-func(0)+s1;
y_right=-func(l)+s1;
candidate_y=[y_left,ymax,y_right];
candidate_x=[0,x0,l];
[y_max,index]=max(candidate_y);
y=[y_max;candidate_x(index)];


end

