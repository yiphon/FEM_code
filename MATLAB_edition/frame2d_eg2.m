%ƽ��ռ�Ԫ������
e=70e9;
a=1e-2;
iz=1e-5;
eles_beam=1;
nodes_beam=2;
node_number=[1;2];  %ÿ����Ԫ��Ӧ�ڵ�ı�š�
node_coordinate=[0;0;4;0];  %ÿ����Ԫ��Ӧ�ڵ�����ꡣ
ele_info=zeros(6,eles_beam);    %��Ԫ��Ϣ����ÿ��������e\a\iz\l\cos\sin��
ele_info(1,:)=e*ones(1,eles_beam);
ele_info(2,:)=a*ones(1,eles_beam);
ele_info(3,:)=iz*ones(1,eles_beam);
ele_info(4:6,:)=frame2d_ele_info(node_coordinate);  %��ȫ��Ԫ��Ϣ�������Ϣ��
r=zeros(6*eles_beam,6);   %rΪ��Ԫת�ƾ���ļ��ϣ���eles�����С�
for i=1:1:eles_beam
    r(6*i-5:6*i,:)=frame2d_r_mat(ele_info(5:6,i));
end     %��ɵ�Ԫת�ƾ���ļ��ϡ�
k_ele=zeros(6*eles_beam,6);  %k_ele�ǵ�Ԫ�ھֲ�����ϵ�µĸնȾ�����eles�����С�
for i=1:1:eles_beam
    k_ele(6*i-5:6*i,:)=frame2d_k_ele_mat(ele_info(1:4,i));
end     %��ɾֲ�����ϵ�µĵ�Ԫ�նȾ��󼯺ϡ�
kk=zeros(6*eles_beam,6); %kk�����������µ�Ԫ�նȾ���ļ��ϣ���eles�����С�
for i=1:1:eles_beam
    kk(6*i-5:6*i,:)=frame2d_ele_stiff_mat(r(6*i-5:6*i,:),k_ele(6*i-5:6*i,:));
end     %�����������ϵ�µ�Ԫ�նȾ��󼯺ϡ�
k=zeros(nodes_beam*3+2,nodes_beam*3+2);   %k������նȾ���
for i=1:1:eles_beam
    k=frame2d_assembly_mat(k,kk(6*i-5:6*i,:),node_number(:,i));
end     %w�������նȾ���ռܲ��֡�
k_truss=ele_stiffness_mat_2dtruss(1,1,0,0,4,3);
k_truss=5*5e6*k_truss;
k([1,2,7,8],[1,2,7,8])=k([1,2,7,8],[1,2,7,8])+k_truss;  %�������նȾ�����ܲ��֡�

 m=k(1:3,1:3);
%  equi=beam_work_equi_distributed(5,-5000);   %equiΪ��Ч�ڵ��غ�����
 %f=[20000;equi(1);equi(2);0;equi(3);equi(4)];
 f=[0;-1e4;0];
 u=m\f;
 U=[u;0;0;0;0;0]   %UΪ�ڵ���λ������
 F=k*U   %FΪ�ڵ�����������
 
 u_frame=zeros(6,eles_beam);   %uΪÿ�ռܵ�Ԫ���ڵ��λ�����С�
 for i=1:1:eles_beam
    u_frame(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3)];
 end
f_frame=zeros(6,eles_beam);    %fΪÿ�ռܵ�Ԫ���ڵ���������У������¶�ÿ����Ԫ�����Ĳ���󣬻�Ӧ��ȥ��Ч�ڵ��غɡ�
for i=1:1:eles_beam
    f_frame(:,i)=k_ele(6*i-5:6*i,:)*r(6*i-5:6*i,:)*u_frame(:,i);    
end
f_frame

u_truss=U([1,2,7,8]);
f_truss=ele_inner_force_2dtruss(1,1,0,0,4,3,u_truss);
f_truss=5*5e6*f_truss;
f_truss

for i=1:1:eles_beam
    beam_deflection_plot(ele_info(4,i),u_frame([2,3,5,6],i),eles_beam,i);
end     %�������Ԫ�Ӷ�ת��ͼ��
