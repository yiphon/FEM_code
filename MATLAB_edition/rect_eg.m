%����Ԫ������
%һ���ı����ĵ�Ȳε�Ԫ������Ҳ���������µķ�������

source=xlsread('rect_eg_input.xlsx');
source(all(isnan(source),2),:)=[];
eles=source(1,1);
nodes=source(1,2);
e=source(2,1);
nu=source(2,2);
t_thickness=source(2,3);
node_number=source(3:2+eles,1:4)';   %ÿ����Ԫ��Ӧ�ڵ�ı�š�
node_coordinate=source(3+eles:2+eles+nodes,1:2)';
node_coordinate=rect_ele_coordinate(node_coordinate,node_number);   %ÿ����Ԫ��Ӧ�ڵ�����ꡣ
unknown_u_index=source(3+eles+nodes,:);
unknown_u_index(:,all(isnan(unknown_u_index),1))=[];
unknown_u_index=unknown_u_index';       %unknown_u_index��U��������δ֪λ��λ�õ�������
known_f=source(4+eles+nodes,1:size(unknown_u_index,1))';     %f_knownΪ��֪�Ľڵ�����������Ч�ڵ��غɡ�
known_u=source(5+eles+nodes,1:2*nodes-size(unknown_u_index,1))';    %u_knownΪ��֪�Ľڵ�λ�ơ�


%����ֵ�����Ļ���˼·��
%   �������ø��ڵ����Ȼ����ֵ��������A��ÿ����[1,si,ti,si*ti]�� 
%   �ⷽ����A*[a1;a2;a3;a4]=[u1;u2;u3;u4],ui�����Ǹ��ڵ������uλ�ƣ�ֻ��һ��Ϊ���㷽������ķ��š�
%   ������uλ�Ʊ�ʾΪ[a1,a2,a3,a4]*[1,s,t,s*t]'�����������������u1\u2\u3\u4��ϵ�������ǵ�ϵ������Ni��
%   ȫ������x��[Ni]*[x1;x2;x3;x4]��ʾ��xi����֪��������ֵ��x\y����s\t��ʾ�ˡ�
%   ����ֻ��һ����Ԫ����������������ȡ�ľ���Ԫ������Ԫ�ߴ���ȫһ��������Ԫ��Ӧ��[Ni]��һ���ġ�
%   ���������������Ԫ�Ľڵ�ֵֻҪ���Ǻ;���Ԫһ����Ȼ���갴�ĸ���1ȡ��[Ni]Ҳ��һ���ġ�


syms s t u1 u2 u3 u4 v1 v2 v3 v4
x_expression=sym(zeros(eles,1));     %x�������������洢����Ԫ�ĳ�����x��s\t��ʾ�ı��ʽ��y���ž����Ӧ��
y_expression=sym(zeros(eles,1));


st_coordinate=sym([1,1;-1,1;-1,-1;1,-1]');       %st_coordinate���������ʾ����Ȼ�����£�ÿ����Ԫ��Ӧ�Ľڵ�����ꡣ
u=sym([u1;u2;u3;u4]);     %u�Է��ű�ʾ��Ԫ�ĸ��ڵ��uλ�ƣ�ֻ��Ϊ�����ڼ������������ֵ������û��ʵ�ʶ�Ӧ����ֵ֪��v��������ͬ��
v=sym([v1;v2;v3;v4]);

st_ele_mat=sym(zeros(4,4));         %st_ele_mat����A*[a1;a2;a3;a4]=[u1;u2;u3;u4]�������е�A��ÿ����[1,si,ti,si*ti]��
for i=1:4
    si=st_coordinate(1,i);
    ti=st_coordinate(2,i);
    st_ele_mat(i,:)=[1,si,ti,si*ti];
end

 
ai=st_ele_mat\u;        %aiΪ[a1;a2;a3;a4]����
u_expression=ai'*sym([1;s;t;s*t]);       %u_expression����s\t��ʾ�ĳ�����u������
Ni=rect_coeffs_ui(u_expression,u)  %����Ni��������ui��ϵ������ȡ��

x_ele=sym(zeros(4,eles));
y_ele=sym(zeros(4,eles));   %x_ele�Ǹ���Ԫ�ڵ����������x�ļ��ϣ�����֪����ֵ��
for i=1:eles
    x_ele(:,i)=node_coordinate([1,3,5,7],i);
    y_ele(:,i)=node_coordinate([2,4,6,8],i);
end


for i=1:eles
    x_expression(i)=Ni'*x_ele(:,i);
    y_expression(i)=Ni'*y_ele(:,i);
end     %x��y��s\t���ʽ�Ѿ��õ���


b=sym(zeros(3,8,eles));     %bΪ����ԪB����ļ��ϡ�
for i=1:eles
    for j=1:4
        b(:,2*j-1:2*j,i)=rect_ele_b_1(Ni(j),x_expression(i),y_expression(i));
    end
end
b
 
 d=e/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
 kk=sym(zeros(8,8,eles));    %kk�Ǹ���Ԫ�նȾ���ļ��ϡ�
 for i=1:1:eles
      kk(:,:,i)=rect_ele_stiff_mat_1(t_thickness,b(:,:,i),d,x_expression(i),y_expression(i));
 end
kk=double(kk);
 
k=zeros(2*nodes,2*nodes);   %kΪ����նȾ���
for i=1:1:eles
    k=rect_assembly_mat(k,kk(:,:,i),node_number(:,i));
end


unknown_var=tri_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
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

u=zeros(8,eles);    %uΪ����Ԫ��Ӧ�ڵ��λ�����󼯺ϡ�
for i=1:1:eles
      u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2);U(node_number(3,i)*2-1:node_number(3,i)*2);
          U(node_number(4,i)*2-1:node_number(4,i)*2)];
end

stress=sym(zeros(3,eles));   %sigmaΪ��ԪӦ�����󣬴�ʱ������ǹ�������ı����������ǳ�����strainΪӦ�䡣
strain=sym(zeros(3,eles));
for i=1:1:eles
    stress(:,i)=d*b(:,:,i)*u(:,i);
    strain(:,i)=b(:,:,i)*u(:,i);
end
num_stress=vpa(stress,6);
num_strain=vpa(strain,6);

stress_zero_point=double(subs(stress,{s,t},{0,0}));
strain_zero_point=double(subs(strain,{s,t},{0,0}));
principal_stress_set=zeros(3,eles);    %�����Ӧ����������
for i=1:1:eles
    principal_stress_set(:,i)=tri_main_sigma(stress_zero_point(:,i));
end
principal_stress=principal_stress_set(1:2,:);
principal_direction=principal_stress_set(3,:);


%����д�ļ���
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
%     6,5,2,3]';  %ÿ����Ԫ��Ӧ�ڵ�ı�ţ�����ʱ��һ���޵����������С�
% length_info=zeros(2,eles);  %��¼ÿ����Ԫ�İ볤�ȡ����ȡ�
% for i=1:1:eles
%     length_info(1,i)=0.125;
%     length_info(2,i)=0.125;
% end     %��ͬ�ķָԪ��ʽ��ÿ����Ԫ��Ӧ��ֵ�ǿ��ܲ�һ���ġ�
% syms x y
% xi_nodes=[1;-1;-1;1];   %xi_nodes��Ծֲ�����ϵ��һ����Ԫ���ڵ�ĺ����꣬ÿ����Ԫ��һ��,�ǳ�ֵ��eta_nodes��Ӧ�����ꡣ
% eta_nodes=[1;1;-1;-1];
% n=sym(zeros(4,eles));   %nΪ����Ԫ��ֵ�����ļ��ϣ�sym�������ݣ�4��eles�С�
% xi_eles=sym(zeros(eles,1)); %xi_eles��Ը�����Ԫ�е������ٺ�����x/a��ÿ����Ԫ��Ϊa���ܲ�һ������Ԫ��Ҳ��һ��һ����sym�������ݡ�eta_eles��Ӧ�����ꡣ
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
% jacob_eles=sym(zeros(4,2*eles));    %���ɸ���Ԫ��Ӧ�Ĳ�ֵ��������ſɱȾ��󼯺ϡ�
% for i=1:1:eles
%     jacob_eles(:,2*i-1:2*i)=jacobian(n(:,i));
% end
% b=sym(zeros(3*eles,8));     %bΪ����ԪB����ļ��ϡ�
% for i=1:1:eles
%     b(3*i-2:3*i,:)=rect_ele_b(jacob_eles(:,2*i-1:2*i));
% end     %����B���󼯺ϡ�
% d=e/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
% kk=sym(zeros(8*eles,8));    %kk�Ǹ���Ԫ�նȾ���ļ��ϡ�
% for i=1:1:eles
%      kk(8*i-7:8*i,:)=rect_ele_stiff_mat(t,length_info(:,i),b(3*i-2:3*i,:),d);
% end
% kk=double(kk);
% k=zeros(2*nodes,2*nodes);   %kΪ����նȾ���
% for i=1:1:eles
%     k=rect_assembly_mat(k,kk(8*i-7:8*i,:),node_number(:,i));
% end
% unknown_u_index=[3,4,5,6,9,10,11,12,15,16]';   %unknown_u_index��δ֪λ����ŵ�������
%  known_f=[zeros(5,1);-12.5e3;zeros(4,1)];
%  known_u=zeros(6,1);
%  unknown_var=tri_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
%  unknown_u=unknown_var{1,1}
%  unknown_f=unknown_var{1,2}
% U=unknown_var{1,3}  %UΪ��λ������
% F=unknown_var{1,4}   %FΪ����������
% u=zeros(8,eles);    %uΪ����Ԫ��Ӧ�ڵ��λ�����󼯺ϡ�
% for i=1:1:eles
%       u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2);U(node_number(3,i)*2-1:node_number(3,i)*2);
%           U(node_number(4,i)*2-1:node_number(4,i)*2)];
% end
% stress=sym(zeros(3,eles));   %sigmaΪ��ԪӦ�����󣬴�ʱ������ǹ�������ı����������ǳ�����strainΪӦ�䡣
% strain=sym(zeros(3,eles));
% for i=1:1:eles
%     stress(:,i)=d*b(3*i-2:3*i,:)*u(:,i);
%     strain(:,i)=b(3*i-2:3*i,:)*u(:,i);
% end
% num_stress=vpa(stress,6)
% num_strain=vpa(strain,6)
% stress_zero_point=double(subs(stress,{x,y},{0,0}))
% strain_zero_point=double(subs(strain,{x,y},{0,0}))
% principal_stress_set=zeros(3,eles);    %�����Ӧ����������
%   for i=1:1:eles
%       principal_stress_set(:,i)=tri_main_sigma(stress_zero_point(:,i));
%   end
%  principal_stress=principal_stress_set(1:2,:)
%  principal_direction=principal_stress_set(3,:)
%  principal_strain_set=zeros(3,eles);    %�����Ӧ����������
%   for i=1:1:eles
%       principal_strain_set(:,i)=rect_main_strain(strain_zero_point(:,i));
%   end
%  principal_strain=principal_strain_set(1:2,:)
%  principal_direction=principal_strain_set(3,:)


 
