%矩形元算例。
%一般四边形四点等参单元理论上也可以用以下的方法处理。

source=xlsread('rect_eg_input.xlsx');
source(all(isnan(source),2),:)=[];
eles=source(1,1);
nodes=source(1,2);
e=source(2,1);
nu=source(2,2);
t_thickness=source(2,3);
node_number=source(3:2+eles,1:4)';   %每个单元对应节点的编号。
node_coordinate=source(3+eles:2+eles+nodes,1:2)';
node_coordinate=rect_ele_coordinate(node_coordinate,node_number);   %每个单元对应节点的坐标。
unknown_u_index=source(3+eles+nodes,:);
unknown_u_index(:,all(isnan(unknown_u_index),1))=[];
unknown_u_index=unknown_u_index';       %unknown_u_index是U列向量中未知位移位置的索引。
known_f=source(4+eles+nodes,1:size(unknown_u_index,1))';     %f_known为已知的节点力，包含等效节点载荷。
known_u=source(5+eles+nodes,1:2*nodes-size(unknown_u_index,1))';    %u_known为已知的节点位移。


%求解插值函数的基本思路：
%   首先利用各节点的自然坐标值构建矩阵A，每行是[1,si,ti,si*ti]。 
%   解方程组A*[a1;a2;a3;a4]=[u1;u2;u3;u4],ui列阵是各节点虚设的u位移，只是一个为计算方便引入的符号。
%   场变量u位移表示为[a1,a2,a3,a4]*[1,s,t,s*t]'，从中再整理符号量u1\u2\u3\u4的系数，他们的系数就是Ni。
%   全局坐标x用[Ni]*[x1;x2;x3;x4]表示，xi是已知的坐标数值，x\y就用s\t表示了。
%   以上只是一个单元的情况。这个例子中取的矩形元，各单元尺寸完全一样，各单元对应的[Ni]是一样的。
%   不规则情况，各单元的节点值只要还是和矩形元一样自然坐标按四个±1取，[Ni]也是一样的。


syms s t u1 u2 u3 u4 v1 v2 v3 v4
x_expression=sym(zeros(eles,1));     %x符号列阵用来存储各单元的场变量x以s\t表示的表达式。y符号矩阵对应。
y_expression=sym(zeros(eles,1));


st_coordinate=sym([1,1;-1,1;-1,-1;1,-1]');       %st_coordinate符号列阵表示在自然坐标下，每个单元对应的节点的坐标。
u=sym([u1;u2;u3;u4]);     %u以符号表示单元四个节点的u位移，只作为符号在计算中用来求插值函数，没有实际对应的已知值。v符号列阵同。
v=sym([v1;v2;v3;v4]);

st_ele_mat=sym(zeros(4,4));         %st_ele_mat是在A*[a1;a2;a3;a4]=[u1;u2;u3;u4]方程组中的A，每行是[1,si,ti,si*ti]。
for i=1:4
    si=st_coordinate(1,i);
    ti=st_coordinate(2,i);
    st_ele_mat(i,:)=[1,si,ti,si*ti];
end

 
ai=st_ele_mat\u;        %ai为[a1;a2;a3;a4]矩阵。
u_expression=ai'*sym([1;s;t;s*t]);       %u_expression是以s\t表示的场变量u函数。
Ni=rect_coeffs_ui(u_expression,u)  %生成Ni，依靠对ui的系数的提取。

x_ele=sym(zeros(4,eles));
y_ele=sym(zeros(4,eles));   %x_ele是各单元节点的总体坐标x的集合，是已知的数值。
for i=1:eles
    x_ele(:,i)=node_coordinate([1,3,5,7],i);
    y_ele(:,i)=node_coordinate([2,4,6,8],i);
end


for i=1:eles
    x_expression(i)=Ni'*x_ele(:,i);
    y_expression(i)=Ni'*y_ele(:,i);
end     %x、y的s\t表达式已经得到。


b=sym(zeros(3,8,eles));     %b为各单元B矩阵的集合。
for i=1:eles
    for j=1:4
        b(:,2*j-1:2*j,i)=rect_ele_b_1(Ni(j),x_expression(i),y_expression(i));
    end
end
b
 
 d=e/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
 kk=sym(zeros(8,8,eles));    %kk是各单元刚度矩阵的集合。
 for i=1:1:eles
      kk(:,:,i)=rect_ele_stiff_mat_1(t_thickness,b(:,:,i),d,x_expression(i),y_expression(i));
 end
kk=double(kk);
 
k=zeros(2*nodes,2*nodes);   %k为总体刚度矩阵。
for i=1:1:eles
    k=rect_assembly_mat(k,kk(:,:,i),node_number(:,i));
end


unknown_var=tri_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %利用分块矩阵计算未知的变量，力和位移。
ua = unknown_var{1,1};
ua=double(ua);
fc = unknown_var{1,2};
fc=double(fc);
U = unknown_var{1,3};
F = unknown_var{1,4};
kcc = unknown_var{1,5};
kca = unknown_var{1,6};
kac = unknown_var{1,7};
kaa = unknown_var{1,8};

u=zeros(8,eles);    %u为各单元对应节点的位移列阵集合。
for i=1:1:eles
      u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2);U(node_number(3,i)*2-1:node_number(3,i)*2);
          U(node_number(4,i)*2-1:node_number(4,i)*2)];
end

stress=sym(zeros(3,eles));   %sigma为单元应力矩阵，此时算出的是关于坐标的变量，不再是常量，strain为应变。
strain=sym(zeros(3,eles));
for i=1:1:eles
    stress(:,i)=d*b(:,:,i)*u(:,i);
    strain(:,i)=b(:,:,i)*u(:,i);
end
num_stress=vpa(stress,6);
num_strain=vpa(strain,6);

stress_zero_point=double(subs(stress,{s,t},{0,0}));
strain_zero_point=double(subs(strain,{s,t},{0,0}));
principal_stress_set=zeros(3,eles);    %求解主应力与主方向。
for i=1:1:eles
    principal_stress_set(:,i)=tri_main_sigma(stress_zero_point(:,i));
end
principal_stress=principal_stress_set(1:2,:);
principal_direction=principal_stress_set(3,:);


%以下写文件。
f=fopen('rect_eg_output.xls','w');

fprintf(f, 'element stiffness matrices:\r\n\r\n');
for i=1:eles
    for j=1:8
        for m=1:8
            fprintf(f,'%8.4f\t',kk(j,m,i));
        end
        fprintf(f, '\r\n');
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\nglobal stiffness matrix:\r\n\r\n');
for j=1:2*nodes
    for m=1:2*nodes
        fprintf(f,'%8.4f\t',k(j,m));
    end
    fprintf(f, '\r\n');
end

fprintf(f,'\r\n\r\nK_cc= \r\n\r\n');
for i=1:size(kcc,1)
    for j=1:size(kcc,2)
        fprintf(f,'%8.4f\t',kcc(i,j));
    end
    fprintf(f,'\r\n');
end
fprintf(f,'\r\n\r\nK_ca= \r\n\r\n');
for i=1:size(kca,1)
    for j=1:size(kca,2)
        fprintf(f,'%8.4f\t',kca(i,j));
    end
    fprintf(f,'\r\n');
end
fprintf(f,'\r\n\r\nK_ac= \r\n\r\n');
for i=1:size(kac,1)
    for j=1:size(kac,2)
        fprintf(f,'%8.4f\t',kac(i,j));
    end
    fprintf(f,'\r\n');
end
fprintf(f,'\r\n\r\nK_aa= \r\n\r\n');
for i=1:size(kaa,1)
    for j=1:size(kaa,2)
        fprintf(f,'%8.4f\t',kaa(i,j));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nU_a= \r\n\r\n');
for i=1:size(ua,1)
    fprintf(f,'%8.4f\t',ua(i));
end
fprintf(f,'\r\n\r\nF_c= \r\n\r\n');
for i=1:size(fc,1)
    fprintf(f,'%8.4f\t',fc(i));
end
fprintf(f,'\r\n\r\nU= \r\n\r\n');
for i=1:size(U,1)
    fprintf(f,'%8.4f\t',U(i));
end
fprintf(f,'\r\n\r\nF= \r\n\r\n');
for i=1:size(F,1)
    fprintf(f,'%8.4f\t',F(i));
end

fprintf(f,'\r\n\r\nstress= \r\n\r\n');
for i=1:size(num_stress,1)
    for j=1:size(num_stress,2)
        fprintf(f,'%s\t\t',char(num_stress(i,j)));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nstrain= \r\n\r\n');
for i=1:size(num_strain,1)
    for j=1:size(num_strain,2)
        fprintf(f,'%s\t\t',char(num_strain(i,j)));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nstress_zero_point= \r\n\r\n');
for i=1:size(stress_zero_point,1)
    for j=1:size(stress_zero_point,2)
        fprintf(f,'%8.4f\t\t',stress_zero_point(i,j));
    end
    fprintf(f,'\r\n');
end


fprintf(f,'\r\n\r\nstrain_zero_point= \r\n\r\n');
for i=1:size(strain_zero_point,1)
    for j=1:size(strain_zero_point,2)
        fprintf(f,'%8.4f\t\t',strain_zero_point(i,j));
    end
    fprintf(f,'\r\n');
end


fprintf(f,'\r\n\r\nprincipal_stress= \r\n\r\n');
for i=1:size(principal_stress,1)
    for j=1:size(principal_stress,2)
        fprintf(f,'%8.4f\t\t',principal_stress(i,j));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nprincipal_direction= \r\n\r\n');
for i=1:size(principal_direction,1)
    for j=1:size(principal_direction,2)
        fprintf(f,'%8.4f\t\t',principal_direction(i,j));
    end
    fprintf(f,'\r\n');
end

fclose(f);









% e=210e9;
% nu=0.3;
% t=0.025;
% eles=3;
% nodes=8;
% node_number=[8,7,4,5;
%     5,4,1,2;
%     6,5,2,3]';  %每个单元对应节点的编号，按逆时针一象限到四象限排列。
% length_info=zeros(2,eles);  %记录每个单元的半长度、半宽度。
% for i=1:1:eles
%     length_info(1,i)=0.125;
%     length_info(2,i)=0.125;
% end     %不同的分割单元方式，每个单元对应的值是可能不一样的。
% syms x y
% xi_nodes=[1;-1;-1;1];   %xi_nodes针对局部坐标系下一个单元各节点的横坐标，每个单元都一样,是常值。eta_nodes对应纵坐标。
% eta_nodes=[1;1;-1;-1];
% n=sym(zeros(4,eles));   %n为各单元插值函数的集合，sym类型内容，4行eles列。
% xi_eles=sym(zeros(eles,1)); %xi_eles针对各个单元中的无量纲横坐标x/a，每个单元因为a可能不一样，各元素也不一定一样。sym类型内容。eta_eles对应纵坐标。
% eta_eles=sym(zeros(eles,1));
% for i=1:1:eles
%     xi_eles(i)=x/length_info(1,i);
%     eta_eles(i)=y/length_info(2,i);
% end
% for i=1:1:4
%     for j=1:1:eles
%     n(i,j)=0.25*(1+xi_nodes(i)*xi_eles(j))*(1+eta_nodes(i)*eta_eles(j));
%     end
% end
% jacob_eles=sym(zeros(4,2*eles));    %生成各单元对应的插值函数组的雅可比矩阵集合。
% for i=1:1:eles
%     jacob_eles(:,2*i-1:2*i)=jacobian(n(:,i));
% end
% b=sym(zeros(3*eles,8));     %b为各单元B矩阵的集合。
% for i=1:1:eles
%     b(3*i-2:3*i,:)=rect_ele_b(jacob_eles(:,2*i-1:2*i));
% end     %生成B矩阵集合。
% d=e/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
% kk=sym(zeros(8*eles,8));    %kk是各单元刚度矩阵的集合。
% for i=1:1:eles
%      kk(8*i-7:8*i,:)=rect_ele_stiff_mat(t,length_info(:,i),b(3*i-2:3*i,:),d);
% end
% kk=double(kk);
% k=zeros(2*nodes,2*nodes);   %k为总体刚度矩阵。
% for i=1:1:eles
%     k=rect_assembly_mat(k,kk(8*i-7:8*i,:),node_number(:,i));
% end
% unknown_u_index=[3,4,5,6,9,10,11,12,15,16]';   %unknown_u_index是未知位移序号的索引。
%  known_f=[zeros(5,1);-12.5e3;zeros(4,1)];
%  known_u=zeros(6,1);
%  unknown_var=tri_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %利用分块矩阵计算未知的变量，力和位移。
%  unknown_u=unknown_var{1,1}
%  unknown_f=unknown_var{1,2}
% U=unknown_var{1,3}  %U为总位移列阵。
% F=unknown_var{1,4}   %F为总受力列阵。
% u=zeros(8,eles);    %u为各单元对应节点的位移列阵集合。
% for i=1:1:eles
%       u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2);U(node_number(3,i)*2-1:node_number(3,i)*2);
%           U(node_number(4,i)*2-1:node_number(4,i)*2)];
% end
% stress=sym(zeros(3,eles));   %sigma为单元应力矩阵，此时算出的是关于坐标的变量，不再是常量，strain为应变。
% strain=sym(zeros(3,eles));
% for i=1:1:eles
%     stress(:,i)=d*b(3*i-2:3*i,:)*u(:,i);
%     strain(:,i)=b(3*i-2:3*i,:)*u(:,i);
% end
% num_stress=vpa(stress,6)
% num_strain=vpa(strain,6)
% stress_zero_point=double(subs(stress,{x,y},{0,0}))
% strain_zero_point=double(subs(strain,{x,y},{0,0}))
% principal_stress_set=zeros(3,eles);    %求解主应力与主方向。
%   for i=1:1:eles
%       principal_stress_set(:,i)=tri_main_sigma(stress_zero_point(:,i));
%   end
%  principal_stress=principal_stress_set(1:2,:)
%  principal_direction=principal_stress_set(3,:)
%  principal_strain_set=zeros(3,eles);    %求解主应变与主方向。
%   for i=1:1:eles
%       principal_strain_set(:,i)=rect_main_strain(strain_zero_point(:,i));
%   end
%  principal_strain=principal_strain_set(1:2,:)
%  principal_direction=principal_strain_set(3,:)


 
