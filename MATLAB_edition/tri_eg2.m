%������Ԫ������
e=200e9;
nu=0.3;
t=0.01;
eles=4;
nodes=5;
node_number=[1,5,2;
    1,3,5;
    3,4,5;
    2,5,4]';  %ÿ����Ԫ��Ӧ�ڵ�ı�ţ�����ʱ�����С�
a1=[0,0.4];
a2=[0.7,0.4];
a3=[0,0];
a4=[0.7,0];
a5=[0.35,0.2];
node_coordinate=[a1,a5,a2;
    a1,a3,a5;
    a3,a4,a5;
    a2,a5,a4]';  %ÿ����Ԫ��Ӧ�ڵ�����ꡣ
a=zeros(eles,1);
for i=1:1:eles
    a(i)=tri_ele_square(node_coordinate(:,i));
end     %���ɸ���Ԫ�������
b=zeros(3*eles,6);  %bΪ����ԪB����ļ��ϡ�
for i=1:1:eles
    b(3*i-2:3*i,:)=tri_ele_b(a(i),node_coordinate(:,i));
end     %���ɸ���ԪB����ļ��ϡ�
d=e/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
kk=zeros(6*eles,6);     %kkΪ����Ԫ�նȾ���ļ��ϡ�
for i=1:1:eles
    kk(6*i-5:6*i,:)=tri_ele_stiff_mat(a(i),t,b(3*i-2:3*i,:),d);
end     %��ɵ�Ԫ�նȾ��󼯺ϡ�
k=zeros(2*nodes+2,2*nodes+2);   %kΪ����նȾ���
for i=1:1:eles
    k=tri_assembly_mat(k,kk(6*i-5:6*i,:),node_number(:,i));
end     %�������նȾ���������Ԫ���֡�
k1=4e6*[1 -1;-1,1];
k2=k1;
dof_k1=[6,11];
dof_k2=[8,12];
for p=1:1:2
    for q=1:1:2
        k(dof_k1(p),dof_k1(q))=k(dof_k1(p),dof_k1(q))+k1(p,q);
        k(dof_k2(p),dof_k2(q))=k(dof_k2(p),dof_k2(q))+k2(p,q);
    end
end     %�������նȾ��󵯻�Ԫ���֡�
 m=k(1:10,1:10);
 f=[0;17.5e3;0;17.5e3;zeros(6,1)];
 u=m\f;
 U=[u;0;0]     %UΪ��λ������
 F=k*U   %FΪ����������
 u=zeros(6,eles);    %uΪ����Ԫ��Ӧ�ڵ��λ�����󼯺ϡ�
 for i=1:1:eles
      u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2);U(node_number(3,i)*2-1:node_number(3,i)*2)];
   end
 sigma=zeros(3,eles);    %sigmaΪÿ��Ԫ��Ӧ�����С�
  for i=1:1:eles
     sigma(:,i)=d*b(3*i-2:3*i,:)*u(:,i);
  end
  sigma
  principal_stress=zeros(3,eles);    %�����Ӧ����������
  for i=1:1:eles
      principal_stress(:,i)=tri_main_sigma(sigma(:,i));
  end
  principal_stress_ele=principal_stress(1:2,:)
  principal_direction=principal_stress(3,:)
 f_bar1=k1*[U(6);U(11)]
 f_bar2=k2*[U(8);U(12)]









