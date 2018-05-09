%ƽ��ռ�Ԫ������
e=1e7;
a=8;
iz=2*4^3/12;
eles=2;
nodes=3;
node_number=[3,2;1,2]';  %ÿ����Ԫ��Ӧ�ڵ�ı�š�
node_coordinate=[0,30,30,30;0,0,30,30]';  %ÿ����Ԫ��Ӧ�ڵ�����ꡣ
ele_info=zeros(6,eles);    %��Ԫ��Ϣ����ÿ��������e\a\iz\l\cos\sin��
ele_info(1,:)=e*ones(1,eles);
ele_info(2,:)=a*ones(1,eles);
ele_info(3,:)=iz*ones(1,eles);
ele_info(4:6,:)=frame2d_ele_info(node_coordinate);  %��ȫ��Ԫ��Ϣ�������Ϣ��
r=zeros(6*eles,6);   %rΪ��Ԫת�ƾ���ļ��ϣ���eles�����С�
for i=1:1:eles
    r(6*i-5:6*i,:)=frame2d_r_mat(ele_info(5:6,i));
end     %��ɵ�Ԫת�ƾ���ļ��ϡ�
k_ele=zeros(6*eles,6);  %k_ele�ǵ�Ԫ�ھֲ�����ϵ�µĸնȾ�����eles�����С�
for i=1:1:eles
    k_ele(6*i-5:6*i,:)=frame2d_k_ele_mat(ele_info(1:4,i));
end     %��ɾֲ�����ϵ�µĵ�Ԫ�նȾ��󼯺ϡ�
kk=zeros(6*eles,6); %kk�����������µ�Ԫ�նȾ���ļ��ϣ���eles�����С�
for i=1:1:eles
    kk(6*i-5:6*i,:)=frame2d_ele_stiff_mat(r(6*i-5:6*i,:),k_ele(6*i-5:6*i,:));
end     %�����������ϵ�µ�Ԫ�նȾ��󼯺ϡ�
kk

k=zeros(nodes*3,nodes*3);   %k������նȾ���
for i=eles:-1:1
    k=frame2d_assembly_mat(k,kk(6*i-5:6*i,:),node_number(:,i));
end     %w�������նȾ���

k
equi=beam_work_equi_distributed(30,-20); 
equi
unknown_u_index=[4,5,6]';   %unknown_u_index��δ֪λ����ŵ�������
 known_f=[0;-1200+equi(3);-1500+equi(4)];
 known_u=zeros(6,1);
 unknown_var=frame2d_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
 unknown_u=unknown_var{1,1}
 unknown_f=unknown_var{1,2}
U=unknown_var{1,3}  %UΪ��λ������
F=unknown_var{1,4}   %FΪ����������
F([8,9,5,6])=F([8,9,5,6])-equi;
f_node_reaction=F([5,6,8,9])


 u=zeros(6,eles);   %uΪÿ��Ԫ���ڵ��λ�����С�
 for i=1:1:eles
    u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3)];
 end
f=zeros(6,eles);    %fΪÿ��Ԫ���ڵ���������У������¶�ÿ����Ԫ�����Ĳ���󣬻�Ӧ��ȥ��Ч�ڵ��غɡ�
for i=1:1:eles
    f(:,i)=k_ele(6*i-5:6*i,:)*r(6*i-5:6*i,:)*u(:,i);
    if i==1
        equi=[0;equi(1);equi(2);0;equi(3);equi(4)];
        f(:,i)=f(:,i)-equi;
    end
end
f


max_sigma=sym(zeros(2,eles));   %���������Ԫ�����Ӧ��ֵ��
for i=1:1:eles
    max_sigma(:,i)=frame2d_max_sigma(ele_info(1,i),ele_info(4,i),u(:,i),f(:,i),ele_info(2,i));
end
max_sigma=vpa(max_sigma,5)


 
 



