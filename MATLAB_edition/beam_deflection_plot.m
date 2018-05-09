function y = beam_deflection_plot( l,u,eles,ele )
%BEAM_DEFLECTION_PLOT Summary of this function goes here
%   Detailed explanation goes here
%   绘出各单元的挠度曲线，u为四行一列。
syms x
deflection=1000*[1-3*x.^2/l^2+2*x.^3/l^3,x-2*x.^2/l+x.^3/l^2,3*x.^2/l^2-2*x.^3/l^3,x.^3/l^2-x.^2/l]*u;
slope=diff(deflection)/1000*180/pi;
x=0:0.01:l;
b=x;
c=x;
b=subs(deflection);
c=subs(slope);
hold on;
subplot(eles,2,2*ele-1);
plot(x,b);
ylabel('Deflection/mm');
subplot(eles,2,2*ele);
plot(x,c);
ylabel('slope/degree');

end

