%轴对称矩形元算例。
source=xlsread('axisy_eg_input.xlsx');
source(all(isnan(source),2),:)=[];
eles=source(1,1);
nodes=source(1,2);
e=source(2,1);
nu=source(2,2);
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
%   全局坐标r用[Ni]*[r1;r2;r3;r4]表示，xi是已知的坐标数值，r\z就用s\t表示了。
%   以上只是一个单元的情况。这个例子中取的矩形元，各单元尺寸完全一样，各单元对应的[Ni]是一样的。
%   不规则情况，各单元的节点值只要还是和矩形元一样自然坐标按四个±1取，[Ni]也是一样的。


syms s t u1 u2 u3 u4
r_expression=sym(zeros(eles,1));     %x符号列阵用来存储各单元的场变量x以s\t表示的表达式。y符号矩阵对应。
z_expression=sym(zeros(eles,1));


st_coordinate=sym([1,1;-1,1;-1,-1;1,-1]');       %st_coordinate符号列阵表示在自然坐标下，每个单元对应的节点的坐标。
u=sym([u1;u2;u3;u4]);     %u以符号表示单元四个节点的u位移，只作为符号在计算中用来求插值函数，没有实际对应的已知值。

st_ele_mat=sym(zeros(4,4));         %st_ele_mat是在A*[a1;a2;a3;a4]=[u1;u2;u3;u4]方程组中的A，每行是[1,si,ti,si*ti]。
for i=1:4
    si=st_coordinate(1,i);
    ti=st_coordinate(2,i);
    st_ele_mat(i,:)=[1,si,ti,si*ti];
end

 
ai=st_ele_mat\u;        %ai为[a1;a2;a3;a4]矩阵。
u_expression=ai'*sym([1;s;t;s*t]);       %u_expression是以s\t表示的场变量u函数。
Ni=rect_coeffs_ui(u_expression,u);  %生成Ni，依靠对ui的系数的提取。

r_ele=sym(zeros(4,eles));
z_ele=sym(zeros(4,eles));   %r_ele是各单元节点的总体坐标r的集合，是已知的数值。
for i=1:eles
    r_ele(:,i)=node_coordinate([1,3,5,7],i);
    z_ele(:,i)=node_coordinate([2,4,6,8],i);
end


for i=1:eles
    r_expression(i)=Ni'*r_ele(:,i);
    z_expression(i)=Ni'*z_ele(:,i);
end     %r、z的s\t表达式已经得到。


b=sym(zeros(4,8,eles));     %b为各单元B矩阵的集合。
for i=1:eles
    for j=1:4
        b(:,2*j-1:2*j,i)=axisy_ele_b_1(Ni(j),r_expression(i),z_expression(i));
    end
end

 d=e*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),nu/(1-nu),0;
                               nu/(1-nu),1,nu/(1-nu),0;
                               nu/(1-nu),nu/(1-nu),1,0;
                               0,0,0,(1-2*nu)/(2-2*nu)];
 d=sym(d);
 
 kk=sym(zeros(8,8,eles));    %kk是各单元刚度矩阵的集合。
 for i=1:1:eles
      kk(:,:,i)=axisy_ele_stiff_mat_1(b(:,:,i),d,r_expression(i),z_expression(i));
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
 
 stress=sym(zeros(4,eles));   %sigma为单元应力矩阵，此时算出的是关于坐标的变量，不再是常量，strain为应变。
 strain=sym(zeros(4,eles));
 for i=1:1:eles
     stress(:,i)=d*b(:,:,i)*u(:,i);
     strain(:,i)=b(:,:,i)*u(:,i);
 end
 num_stress=vpa(stress,6);
 num_strain=vpa(strain,6);
 
stress_zero_point=double(subs(stress,{s,t},{0,0}));
strain_zero_point=double(subs(strain,{s,t},{0,0}));

%以下写文件。
f=fopen('axisy_eg_output.xls','w');

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


fclose(f);

  








% e=200e9;
% nu=0.3;
% eles=4;
% nodes=9;
% node_number=[6,1,5,7;
%     7,5,2,8;
%     9,7,8,3;
%     4,6,7,9]';  %每个单元对应节点的编号，按逆时针、一到四象限排列。
% node_coordinate=0.001*[40,10;40,0;60,0;60,10;
%     40,5;50,10;50,5;50,0;60,5]';    %每个节点的坐标。
% 
% syms u1 v1 u2 v2 u3 v3 u4 v4 u5 v5 u6 v6 u7 v7 u8 v8 u9 v9 r z
% u_sym=sym([u1;u2;u3;u4;u5;u6;u7;u8;u9]);
% v_sym=sym([v1;v2;v3;v4;v5;v6;v7;v8;v9]);
% 
% node_info_for_n=zeros(nodes,4);   %node_info_for_n是各节点[1;ri;zi;ri*zi]的集合，用以求形状函数矩阵N时用。
% for i=1:1:nodes
%     node_info_for_n(i,:)=axisy_node_info_for_n(node_coordinate(:,i));
% end
% 
% ele_info_for_n=zeros(4*eles,4);     %ele_info_for_n是各单元对应节点[1;ri;zi;ri*zi]的集合，用以求形状函数矩阵N时用。
% for i=1:1:eles
%     ele_info_for_n(4*i-3:4*i,:)=axisy_ele_info_for_n(node_info_for_n,node_number(:,i));
% end
% 
% ele_u_sym=sym(zeros(4,eles));  %ele_u_sym为各单元为求解N中函数系数设置的各节点坐标，sym类，视作已知量。
% ele_v_sym=sym(zeros(4,eles));
% for i=1:1:eles
%     ele_u_sym(:,i)=axisy_ele_uv_sym(u_sym,node_number(:,i));
%     ele_v_sym(:,i)=axisy_ele_uv_sym(v_sym,node_number(:,i));
% end
% 
% n_coeffs_for_rz=sym(zeros(8,eles));     %n_coeffs_for_rz是各单元求解出来的位移插值形式以r\z为整理得到的系数，前四行对应u位移，后四行对应r位移。
% for i=1:1:eles
%     n_coeffs_for_rz(1:4,i)=axisy_n_coeffs_for_rz(ele_u_sym(:,i),ele_info_for_n(4*i-3:4*i,:));
%     n_coeffs_for_rz(5:8,i)=axisy_n_coeffs_for_rz(ele_v_sym(:,i),ele_info_for_n(4*i-3:4*i,:));
% end
% 
% uv_functions=sym(zeros(eles,1));  %uv_functions是最终得到的单元内任意一点位移的的函数。
% for i=1:1:eles
%     uv_functions(i)=sym([1,r,z,r*z])*n_coeffs_for_rz(1:4,i);
% end
% uv_functions
% 
% n=sym(zeros(eles,4));   %n为各单元对应的位移函数矩阵的集合，即uv_functions以ui\vi整理得到的系数，即插值函数。
% for i=1:1:eles
%     p=ele_u_sym(:,i);
%     uv_functions(i)=collect(uv_functions(i),p);
%     p=p(end:-1:1);
%     n(i,:)=coeffs(uv_functions(i),p);
% end
% n
% 
% jacob_eles=sym(zeros(4,2*eles));    %生成各单元对应的插值函数组的雅可比矩阵集合。
%  for i=1:1:eles
%      jacob_eles(:,2*i-1:2*i)=jacobian(n(i,:));
%  end
%  jacob_eles
%  
%  b=sym(zeros(4*eles,8));     %b为各单元B矩阵的集合。
%  for i=1:1:eles
%      b(4*i-3:4*i,:)=axisy_ele_b(jacob_eles(:,2*i-1:2*i),n(i,:));
%  end     %生成B矩阵集合。
%  b
%  
%  d=e*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),nu/(1-nu),0;
%                                nu/(1-nu),1,nu/(1-nu),0;
%                                nu/(1-nu),nu/(1-nu),1,0;
%                                0,0,0,(1-2*nu)/(2-2*nu)];
%  d=sym(d);
%  
%  kk=sym(zeros(8*eles,8));    %kk是各单元刚度矩阵的集合。
%  for i=1:1:eles
%       kk(8*i-7:8*i,:)=axisy_ele_stiff_mat(node_number(:,i),node_coordinate,b(4*i-3:4*i,:),d);
%  end
%  kk
%  kk=double(kk)
%  
%  k=zeros(2*nodes,2*nodes);   %k为总体刚度矩阵。
%  for i=1:1:eles
%      k=rect_assembly_mat(k,kk(8*i-7:8*i,:),node_number(:,i));
%  end
% k
% 
%  unknown_u_index=[1;3;9;11;13;15];   %unknown_u_index是未知位移序号的索引。
%  known_f=[2513/2;2513/2;2513;0;0;0];
%  known_u=zeros(12,1);
%  unknown_var=tri_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %利用分块矩阵计算未知的变量，力和位移。
%  unknown_u=unknown_var{1,1}
%  unknown_f=unknown_var{1,2}
%  U=unknown_var{1,3}  %U为总位移列阵。
%  F=unknown_var{1,4}   %F为总受力列阵。
%   
%  u=zeros(8,eles);    %u为各单元对应节点的位移列阵集合。
%  for i=1:1:eles
%        u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2);U(node_number(3,i)*2-1:node_number(3,i)*2);
%            U(node_number(4,i)*2-1:node_number(4,i)*2)];
%  end
%  
%  stress=sym(zeros(4,eles));   %sigma为单元应力矩阵，此时算出的是关于坐标的变量，不再是常量，strain为应变。
%  strain=sym(zeros(4,eles));
%  for i=1:1:eles
%      stress(:,i)=d*b(4*i-3:4*i,:)*u(:,i);
%      strain(:,i)=b(4*i-3:4*i,:)*u(:,i);
%  end
%  num_stress=vpa(stress,6)
%  num_strain=vpa(strain,6)
%  
%  eles_center=zeros(2,eles);  %确定每个单元中心的坐标。
%   for i=1:1:eles
%      eles_center(:,i)=[(node_coordinate(1,node_number(1,i))+node_coordinate(1,node_number(3,i)))/2;
%          (node_coordinate(2,node_number(1,i))+node_coordinate(2,node_number(3,i)))/2];
%   end
%   eles_center
%   
%   stress_zero_point=sym(zeros(4,eles));
%   strain_zero_point=sym(zeros(4,eles)); %用以求单元中心点的应力应变。
%   for i=1:1:eles
%       for j=1:1:4
%       if ~isempty(findstr('r',findsym(stress(j,i)))) & ~isempty(findstr('z',findsym(stress(j,i))))
%           stress_zero_point(j,i)=subs(stress(j,i),{r,z},{eles_center(1,i),eles_center(2,i)});
%       elseif isempty(findstr('r',findsym(stress(j,i)))) & isempty(findstr('z',findsym(stress(j,i))))
%           stress_zero_point(j,i)=stress(j,i);
%       elseif isempty(findstr('z',findsym(stress(j,i))))
%           stress_zero_point(j,i)=subs(stress(j,i),r,eles_center(1,i));
%       elseif isempty(findstr('r',findsym(stress(j,i))))
%           stress_zero_point(j,i)=subs(stress(j,i),z,eles_center(2,i));
%       end
%       end
%   end
%   for i=1:1:eles
%       for j=1:1:4
%       if ~isempty(findstr('r',findsym(strain(j,i)))) & ~isempty(findstr('z',findsym(strain(j,i))))
%           strain_zero_point(j,i)=subs(strain(j,i),{r,z},{eles_center(1,i),eles_center(2,i)});
%       elseif isempty(findstr('r',findsym(strain(j,i)))) & isempty(findstr('z',findsym(strain(j,i))))
%           strain_zero_point(j,i)=strain(j,i);
%       elseif isempty(findstr('z',findsym(strain(j,i))))
%           strain_zero_point(j,i)=subs(strain(j,i),r,eles_center(1,i));
%       elseif isempty(findstr('r',findsym(strain(j,i))))
%           strain_zero_point(j,i)=subs(strain(j,i),z,eles_center(2,i));
%       end
%       end
%   end
%   
%    stress_zero_point=double(stress_zero_point)       
%    strain_zero_point=double(strain_zero_point)
