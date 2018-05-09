%  ��Ԫ�㷨ʵ����
e1=207e9;
I=40e-3^4/12;
ele_info=[e1,e1,;I,I;0.3,0.3];   %���ɵ�Ԫ��Ϣ����
eles_beam=2;
nodes_beam=3;
node_number=[1,2;2,3];
kk=zeros(4*eles_beam,4); %kk�Ǹ���Ԫ�նȾ���ļ��ϡ�
for i=1:1:eles_beam
    kk(4*i-3:4*i,:)=beam_ele_stiff_mat(ele_info(:,i));
end     %���ɵ�Ԫ�նȾ���
kk_bar=69e9*78.54e-6/0.2*[1 -1;-1 1]; %kk_bar�ǵ��ɵĸնȾ���
k=zeros(nodes_beam*2+1,nodes_beam*2+1);   %kΪ����նȾ���
for i=1:1:eles_beam
    k=beam_assembly_mat(k,kk(4*i-3:4*i,:),node_number(:,i));
end
k(3,3)=k(3,3)+kk_bar(1,1);
k(3,7)=k(3,7)+kk_bar(1,2);
k(7,3)=k(7,3)+kk_bar(2,1);
k(7,7)=k(7,7)+kk_bar(2,2);  %��������նȾ���
% equi=zeros(4,2);    %equiΪ��Ч�ڵ��غɾ���
% equi(:,1)=beam_work_equi_distributed(3);
% equi(:,2)=beam_work_equi_central(2,4,-30e3);
% f=[equi(4,1);equi(2,2);equi(4,2)];
m=k(2:6,2:6);
f=[0;0;0;-1e4;0];
u=m\f;  %���λ�õ�����λ�ơ�
U=[0;u;0]   %UΪ��λ�ƾ���
F=k*U   %FΪ�ܽڵ�������
u=zeros(4,eles_beam);    %uΪ����Ԫ���ڵ��λ�ƾ���
for i=1:1:eles_beam
    u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2)];
end
for i=1:1:eles_beam
   beam_deflection_plot(ele_info(3,i),u(:,i),3,i);
end     %�������Ԫ�Ӷ���ת��ͼ��
f_beam=zeros(4,eles_beam);    %fΪ����Ԫ���ڵ��������
for i=1:1:eles_beam
    f_beam(:,i)=beam_node_force(kk(4*i-3:4*i,:),u(:,i));
end
f_bar=kk_bar*[U(3);U(7)];
f_beam
f_bar


    