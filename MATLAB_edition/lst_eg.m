%���ڵ�������Ԫ������

s=xlsread('lst_eg_input.xlsx');
s(all(isnan(s),2),:)=[];
eles=s(1,1);
nodes=s(1,2);
e=s(2,1);
nu=s(2,2);
t_thickness=s(2,3);
node_number=s(3:2+eles,1:6)';   %ÿ����Ԫ��Ӧ�ڵ�ı�š�

node_coordinate=zeros(2,3,eles);    %node-coordinate��ʱ��Ӧ���˵�ڵ�����ꡣ
for i=1:eles
    node_coordinate(:,:,i)=s(3*i+eles:3*i+2+eles,1:2)';
    end
node_coordinate;
node_coordinate_all=zeros(2,6,eles);    %node_coordinate_all��Ӧ����Ԫ�����������ڵ㡣
for i=1:eles
    node_coordinate_all(:,1:3,i)=node_coordinate(:,:,i);
    node_coordinate_all(:,4,i)=0.5*(node_coordinate(:,1,i)+node_coordinate(:,2,i));
    node_coordinate_all(:,5,i)=0.5*(node_coordinate(:,2,i)+node_coordinate(:,3,i));
    node_coordinate_all(:,6,i)=0.5*(node_coordinate(:,1,i)+node_coordinate(:,3,i));
end
node_coordinate=node_coordinate_all;    %node_coordinate��ʱ��Ӧ����Ԫ�����������ڵ㡣

unknown_u_index=s(3+eles*4,:);
unknown_u_index(:,all(isnan(unknown_u_index),1))=[];
unknown_u_index=unknown_u_index';       %unknown_u_index��U��������δ֪λ��λ�õ�������
known_f=s(4+eles*4,1:size(unknown_u_index,1))';     %f_knownΪ��֪�Ľڵ�����������Ч�ڵ��غɡ�
known_u=s(5+eles*4,1:2*nodes-size(unknown_u_index,1))';    %u_knownΪ��֪�Ľڵ�λ�ơ�

syms x y

l=sym(zeros(3,eles));       %lΪ����Ԫl1\l2\l3�ı��ʽ���ϡ�
n=sym(zeros(6,eles));       %nΪ����Ԫn1\...n6�ı��ʽ�ļ��ϡ�
for i=1:eles
    l(:,i)=lst_l(node_coordinate(:,1:3,i));
end

for i=1:eles
    n(:,i)=lst_n(l(:,i));
end

b=sym(zeros(3,12,eles));      %bΪ����ԪB����ļ��ϡ�
for i=1:1:eles
    for j=1:6
    b(:,2*j-1:2*j,i)=lst_ele_b(n(j,i));
    end
end     %���ɸ���ԪB����ļ��ϡ�

d=e/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
kk=zeros(12,12,eles);     %kkΪ����Ԫ�նȾ���ļ��ϡ�
for i=1:1:eles
    kk(:,:,i)=lst_ele_stiff_mat(node_coordinate(:,1:3,i),t_thickness,b(:,:,i),d);
end     %��ɵ�Ԫ�նȾ��󼯺ϡ�

k=zeros(2*nodes,2*nodes);   %kΪ����նȾ���
for i=1:1:eles
    k=lst_assembly_mat(k,kk(:,:,i),node_number(:,i));
end     %�������նȾ���

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

u=zeros(12,eles);    %uΪ����Ԫ��Ӧ�ڵ��λ�����󼯺ϡ�
for i=1:1:eles
     u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2);U(node_number(3,i)*2-1:node_number(3,i)*2);
        U(node_number(4,i)*2-1:node_number(4,i)*2);U(node_number(5,i)*2-1:node_number(5,i)*2);
        U(node_number(6,i)*2-1:node_number(6,i)*2)];
end

stress=sym(zeros(3,eles));    %stressΪÿ��Ԫ��Ӧ�����С�
 for i=1:1:eles
    stress(:,i)=d*b(:,:,i)*u(:,i);
 end
 stress=vpa(stress,5);
 
 stress_center=zeros(3,eles);   %stress_center����ÿ��Ԫ���ĵ�Ӧ��ֵ��
 for i=1:eles
     stress_center(:,i)=lst_stress_center(node_coordinate(:,1:3,i),stress(:,i));
 end
 
 main_stress=zeros(3,eles);    %�����Ӧ����������
 for i=1:1:eles
     main_stress(:,i)=tri_main_sigma(stress_center(:,i));
 end
principal_stress=main_stress(1:2,:);
principal_direction=main_stress(3,:);


%����д�ļ���
f=fopen('lst_eg_output.xls','w');

fprintf(f, 'element stiffness matrices:\r\n\r\n');
for i=1:eles
    for j=1:12
        for m=1:12
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
for i=1:size(stress,1)
    for j=1:size(stress,2)
        fprintf(f,'%s\t\t',char(stress(i,j)));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nstress_center= \r\n\r\n');
for i=1:size(stress_center,1)
    for j=1:size(stress_center,2)
        fprintf(f,'%8.4f\t\t',stress_center(i,j));
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
