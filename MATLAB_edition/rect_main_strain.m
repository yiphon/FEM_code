function y = rect_main_strain( strain )
%RECT_MAIN_STRAIN Summary of this function goes here
%   Detailed explanation goes here
%   �������Ӧ���������򣬶Ծ���Ԫ����ԭ��ģ�strainΪ[sx;sy;rxy]����
sx=strain(1);
sy=strain(2);
rxy=strain(3)/2;
a=zeros(3,1);
a(1)=(sx+sy)/2+sqrt(((sx-sy)/2)^2+rxy^2);
a(2)=(sx+sy)/2-sqrt(((sx-sy)/2)^2+rxy^2);
a(3)=0.5*atan2(2*rxy,sx-sy)*180/pi;
y=a;
end

