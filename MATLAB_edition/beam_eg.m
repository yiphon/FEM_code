%  ��Ԫ�㷨ʵ����
s=xlsread('beam_eg_input.xlsx');
s(all(isnan(s),2),:)=[];
eles=s(1,1);
nodes=s(2,1);
e=s(3,:);
i=s(4,:);
l=s(5,:);
node_number=s(6:5+eles,1:2)';
ele_info=[e;i;l];       %���ɵ�Ԫ��Ϣ����
dis_f_val=s(6+eles,1);    %����ֲ������õ���Ϣ�� 
dis_f_ele=s(6+eles,2);
dis_f_range=s(6+eles,3);
cen_f_val= s(7+eles,1);    %����ǽڵ㼯����������Ϣ��
cen_f_ele= s(7+eles,2);
cen_f_pos = s(7+eles,3);

 kk=zeros(4,4,eles); %kk�Ǹ���Ԫ�նȾ���ļ��ϡ�
 for i=1:1:eles
     kk(:,:,i)=beam_ele_stiff_mat(ele_info(:,i));
 end     %���ɵ�Ԫ�նȾ���
 kk

 k=zeros(nodes*2,nodes*2);   %kΪ����նȾ���
 for i=1:1:eles
     k=beam_assembly_mat(k,kk(:,:,i),node_number(:,i));
 end     %��������նȾ���
 k
 
equi_dis=zeros(4,1);
equi_cen=zeros(4,1);    %equi_Ϊ��Ч�ڵ��غɾ���
equi_dis=beam_work_equi_distributed(dis_f_range, dis_f_val);
equi_cen=beam_work_equi_central(cen_f_pos,l(cen_f_ele),cen_f_val);

unknown_u_index=[4;6;8];    %unknown_u_index��U��������δ֪λ��λ�õ�������
f_known=[equi_dis(4);equi_cen(2);equi_cen(4)];  %f_knownΪ��֪�Ľڵ�����������Ч�ڵ��غɡ�
u_known=zeros(5,1);     %u_knownΪ��֪�Ľڵ�λ�ơ�
unknown_var=tri_solve_unknown(k,unknown_u_index,f_known,u_known,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
ua = unknown_var{1,1};
ua=double(ua)
fc = unknown_var{1,2};
fc=double(fc)
U = unknown_var{1,3}
F = unknown_var{1,4}
kcc = unknown_var{1,5}
kca = unknown_var{1,6}
kac = unknown_var{1,7}
kaa = unknown_var{1,8}

 u=zeros(4,eles);    %uΪ����Ԫ���ڵ��λ�ƾ���
 for i=1:1:eles
     u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2)];
 end
 
 for i=1:1:eles
    beam_deflection_plot(ele_info(3,i),u(:,i),eles,i);
 end     %�������Ԫ�Ӷ���ת��ͼ��

 f_nodal=zeros(4,eles);    %fΪ����Ԫ���ڵ������������Ԫ�ڲ�����ʱ����Ҫ��ȥ��Ч�ڵ��غɡ�
 for i=1:1:eles
     f_nodal(:,i)=beam_node_force(kk(:,:,i),u(:,i));
 end
 f_nodal(:,dis_f_ele)=f_nodal(:,dis_f_ele)-equi_dis;
 f_nodal(:,cen_f_ele)=f_nodal(:,cen_f_ele)-equi_cen;
 f_nodal

 
 %����д�ļ���
 f=fopen('beam_eg_output.xls','w');
    
 fprintf(f, 'element stiffness matrices:\r\n\r\n');   
 for i=1:eles
     for j=1:4
        for m=1:4
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

fprintf(f,'\r\n\r\nf_nodal= \r\n\r\n');
for i=1:size(f_nodal,1)
    for j=1:size(f_nodal,2)
        fprintf(f,'%8.4f\t',f_nodal(i,j));
    end
    fprintf(f,'\r\n');
end

fclose(f);






% e=210e9;
% I=50e-6;
% ele_info=[e,e,e;I,I,I;3,3,4];   %���ɵ�Ԫ��Ϣ����
% eles=3;
% nodes=4;
% node_number=[1,2,3;2,3,4];
% kk=zeros(4*eles,4); %kk�Ǹ���Ԫ�նȾ���ļ��ϡ�
% for i=1:1:eles
%     kk(4*i-3:4*i,:)=beam_ele_stiff_mat(ele_info(:,i));
% end     %���ɵ�Ԫ�նȾ���
% k=zeros(nodes*2,nodes*2);   %kΪ����նȾ���
% for i=1:1:eles
%     k=beam_assembly_mat(k,kk(4*i-3:4*i,:),node_number(:,i));
% end     %��������նȾ���
% equi=zeros(4,2);    %equiΪ��Ч�ڵ��غɾ���
% equi(:,1)=beam_work_equi_distributed(3,-10000);
% equi(:,2)=beam_work_equi_central(2,4,-30e3);
% f=[equi(4,1);equi(2,2);equi(4,2)];
% m=k([4,6,8],[4,6,8]);
% u=m\f;  %���λ�õ�����λ�ơ�
% U=[0;0;0;u(1);0;u(2);0;u(3)]   %UΪ��λ�ƾ���
% F=k*U   %FΪ�ܽڵ�������,��δ��ȥ��Ч�ڵ��غɣ������Ժ󣬾��Ǹ��ڵ��ܵĺ�����
% u=zeros(4,eles);    %uΪ����Ԫ���ڵ��λ�ƾ���
% for i=1:1:eles
%     u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2)];
% end
% for i=1:1:eles
%    beam_deflection_plot(ele_info(3,i),u(:,i),eles,i);
% end     %�������Ԫ�Ӷ���ת��ͼ��
% f=zeros(4,eles);    %fΪ����Ԫ���ڵ������������Ԫ�ڲ�����ʱ����Ҫ��ȥ��Ч�ڵ��غɡ�
% for i=1:1:eles
%     f(:,i)=beam_node_force(kk(4*i-3:4*i,:),u(:,i));
% end
% f(:,1)=f(:,1)-equi(:,1);
% f(:,3)=f(:,3)-equi(:,2);
% f
% stress_mat=beam_strain(ele_info(1,2),-0.02,2,u(:,2),ele_info(3,2));
% strain=stress_mat(1)
% stress=stress_mat(2)    %��ĳλ�õ�Ӧ��Ӧ��ֵ��



    