function y = frame2d_r_mat( ele_info )
%FRAME2D_R_MAT Summary of this function goes here
%   Detailed explanation goes here
%   ����ÿ����Ԫ��ת�ƾ���,ele_info��[cos;sin]��ʽ��
c=ele_info(1);
s=ele_info(2);
rr=[c,s,0;-s,c,0;0,0,1];
r=[rr,zeros(3,3);zeros(3,3),rr];
y=r;
end

