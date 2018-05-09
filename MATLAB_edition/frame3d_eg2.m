%����ռ�Ԫ������
e=207e9;
g=84e9;
a=18e-4;
iy=0.06*0.03^3/12;
iz=0.03*0.06^3/12;
j=5e-5;
eles=1;
nodes=2;
node_number=[1;2];  %ÿ����Ԫ��Ӧ�ڵ�ı�š�
node_coordinate=[0,0,0,1.5,0,0]';  %ÿ����Ԫ��Ӧ�ڵ�����ꡣ
ele_info=zeros(10,eles);    %��Ԫ��Ϣ����ÿ��������e\a\iy\iz\j\l\cXx\cYx\cZx\g��
ele_info(1,:)=e*ones(1,eles);
ele_info(2,:)=a*ones(1,eles);
ele_info(3,:)=iy*ones(1,eles);
ele_info(4,:)=iz*ones(1,eles);
ele_info(5,:)=j*ones(1,eles);
ele_info(10,:)=g*ones(1,eles);
ele_info(6:9,:)=frame3d_ele_info(node_coordinate);  %��ȫ��Ԫ��Ϣ�������Ϣ��
r=zeros(12*eles,12);   %rΪ��Ԫת�ƾ���ļ��ϣ�ʮ��eles��12�С�
for i=1:1:eles
    r(12*i-11:12*i,:)=frame3d_r_mat(ele_info(7:9,i));
end     %��ɵ�Ԫת�ƾ���ļ��ϡ�
k_ele=zeros(12*eles,12);  %k_ele�ǵ�Ԫ�ھֲ�����ϵ�µĸնȾ���ʮ��eles��12�С�
for i=1:1:eles
    k_ele(12*i-11:12*i,:)=frame3d_k_ele_mat(ele_info(:,i));
end     %��ɾֲ�����ϵ�µĵ�Ԫ�նȾ��󼯺ϡ�
kk=zeros(12*eles,12); %kk�����������µ�Ԫ�նȾ���ļ��ϣ�12eles��12�С�
for i=1:1:eles
    kk(12*i-11:12*i,:)=frame3d_ele_stiff_mat(r(12*i-11:12*i,:),k_ele(12*i-11:12*i,:));
end     %�����������ϵ�µ�Ԫ�նȾ��󼯺ϡ�
k=zeros(nodes*6,nodes*6);   %k������նȾ���
for i=1:1:eles
    k=frame3d_assembly_mat(k,kk(12*i-11:12*i,:),node_number(:,i));
end     %w�������նȾ���
 m=k(7:12,7:12);
 f=[0;500;300;0;0;0];
 u=m\f;
 U=[zeros(6,1);u]   %UΪ�ڵ���λ������
 F=k*U   %FΪ�ڵ�����������
 u=zeros(12,eles);   %uΪÿ��Ԫ���ڵ��λ�����С�
 for i=1:1:eles
    u(:,i)=[U(node_number(1,i)*6-5:node_number(1,i)*6);U(node_number(2,i)*6-5:node_number(2,i)*6)];
 end
f=zeros(12,eles);    %fΪÿ��Ԫ���ڵ���������У������¶�ÿ����Ԫ�����Ĳ���󣬻�Ӧ��ȥ��Ч�ڵ��غɡ�
for i=1:1:eles
    f(:,i)=k_ele(12*i-11:12*i,:)*r(12*i-11:12*i,:)*u(:,i);
end
f
max_bending_sigma=sym(zeros(2,eles));   %���������Ԫ���������Ӧ��ֵ��
for i=1:1:eles
    max_bending_sigma(:,i)=frame3d_max_bending_sigma(ele_info(1,i),ele_info(6,i),u(:,i));
end
max_bending_sigma=vpa(max_bending_sigma,5)
