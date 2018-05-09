function y = tri_main_sigma( sigma )
%TRI_MAIN_SIGMA Summary of this function goes here
%   Detailed explanation goes here
%   求各单元的主应力和最大主应力的主方向。
sx=sigma(1);
sy=sigma(2);
txy=sigma(3);
a=zeros(3,1);
a(1)=(sx+sy)/2+sqrt(((sx-sy)/2)^2+txy^2);
a(2)=(sx+sy)/2-sqrt(((sx-sy)/2)^2+txy^2);
a(3)=0.5*atan2(2*txy,sx-sy)*180/pi;
y=a;
end

