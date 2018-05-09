function y = rect_main_strain( strain )
%RECT_MAIN_STRAIN Summary of this function goes here
%   Detailed explanation goes here
%   求最大主应变与主方向，对矩形元是求原点的，strain为[sx;sy;rxy]。。
sx=strain(1);
sy=strain(2);
rxy=strain(3)/2;
a=zeros(3,1);
a(1)=(sx+sy)/2+sqrt(((sx-sy)/2)^2+rxy^2);
a(2)=(sx+sy)/2-sqrt(((sx-sy)/2)^2+rxy^2);
a(3)=0.5*atan2(2*rxy,sx-sy)*180/pi;
y=a;
end

