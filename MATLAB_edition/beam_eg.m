%  梁元算法实例。
s=xlsread('beam_eg_input.xlsx');
s(all(isnan(s),2),:)=[];
eles=s(1,1);
nodes=s(2,1);
e=s(3,:);
i=s(4,:);
l=s(5,:);
node_number=s(6:5+eles,1:2)';
ele_info=[e;i;l];       %生成单元信息矩阵。
dis_f_val=s(6+eles,1);    %读入分布力作用的信息。 
dis_f_ele=s(6+eles,2);
dis_f_range=s(6+eles,3);
cen_f_val= s(7+eles,1);    %读入非节点集中力作用信息。
cen_f_ele= s(7+eles,2);
cen_f_pos = s(7+eles,3);

 kk=zeros(4,4,eles); %kk是各单元刚度矩阵的集合。
 for i=1:1:eles
     kk(:,:,i)=beam_ele_stiff_mat(ele_info(:,i));
 end     %生成单元刚度矩阵。
 kk

 k=zeros(nodes*2,nodes*2);   %k为整体刚度矩阵。
 for i=1:1:eles
     k=beam_assembly_mat(k,kk(:,:,i),node_number(:,i));
 end     %生成整体刚度矩阵。
 k
 
equi_dis=zeros(4,1);
equi_cen=zeros(4,1);    %equi_为等效节点载荷矩阵。
equi_dis=beam_work_equi_distributed(dis_f_range, dis_f_val);
equi_cen=beam_work_equi_central(cen_f_pos,l(cen_f_ele),cen_f_val);

unknown_u_index=[4;6;8];    %unknown_u_index是U列向量中未知位移位置的索引。
f_known=[equi_dis(4);equi_cen(2);equi_cen(4)];  %f_known为已知的节点力，包含等效节点载荷。
u_known=zeros(5,1);     %u_known为已知的节点位移。
unknown_var=tri_solve_unknown(k,unknown_u_index,f_known,u_known,nodes);     %利用分块矩阵计算未知的变量，力和位移。
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

 u=zeros(4,eles);    %u为各单元两节点的位移矩阵。
 for i=1:1:eles
     u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2)];
 end
 
 for i=1:1:eles
    beam_deflection_plot(ele_info(3,i),u(:,i),eles,i);
 end     %绘出各单元挠度与转角图像。

 f_nodal=zeros(4,eles);    %f为各单元两节点的力矩阵，做单元内部受力时，需要减去等效节点载荷。
 for i=1:1:eles
     f_nodal(:,i)=beam_node_force(kk(:,:,i),u(:,i));
 end
 f_nodal(:,dis_f_ele)=f_nodal(:,dis_f_ele)-equi_dis;
 f_nodal(:,cen_f_ele)=f_nodal(:,cen_f_ele)-equi_cen;
 f_nodal

 
 %以下写文件。
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
% ele_info=[e,e,e;I,I,I;3,3,4];   %生成单元信息矩阵。
% eles=3;
% nodes=4;
% node_number=[1,2,3;2,3,4];
% kk=zeros(4*eles,4); %kk是各单元刚度矩阵的集合。
% for i=1:1:eles
%     kk(4*i-3:4*i,:)=beam_ele_stiff_mat(ele_info(:,i));
% end     %生成单元刚度矩阵。
% k=zeros(nodes*2,nodes*2);   %k为整体刚度矩阵。
% for i=1:1:eles
%     k=beam_assembly_mat(k,kk(4*i-3:4*i,:),node_number(:,i));
% end     %生成整体刚度矩阵。
% equi=zeros(4,2);    %equi为等效节点载荷矩阵。
% equi(:,1)=beam_work_equi_distributed(3,-10000);
% equi(:,2)=beam_work_equi_central(2,4,-30e3);
% f=[equi(4,1);equi(2,2);equi(4,2)];
% m=k([4,6,8],[4,6,8]);
% u=m\f;  %求出位置的三个位移。
% U=[0;0;0;u(1);0;u(2);0;u(3)]   %U为总位移矩阵。
% F=k*U   %F为总节点力矩阵,但未减去等效节点载荷，整合以后，就是各节点受的合力。
% u=zeros(4,eles);    %u为各单元两节点的位移矩阵。
% for i=1:1:eles
%     u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2)];
% end
% for i=1:1:eles
%    beam_deflection_plot(ele_info(3,i),u(:,i),eles,i);
% end     %绘出各单元挠度与转角图像。
% f=zeros(4,eles);    %f为各单元两节点的力矩阵，做单元内部受力时，需要减去等效节点载荷。
% for i=1:1:eles
%     f(:,i)=beam_node_force(kk(4*i-3:4*i,:),u(:,i));
% end
% f(:,1)=f(:,1)-equi(:,1);
% f(:,3)=f(:,3)-equi(:,2);
% f
% stress_mat=beam_strain(ele_info(1,2),-0.02,2,u(:,2),ele_info(3,2));
% strain=stress_mat(1)
% stress=stress_mat(2)    %求某位置的应力应变值。



    