%������ܽṹ����Ԫ�㷨ʵ����
eles_str=['Input the number of the truss elements:' 10];
eles=input(eles_str)     %eles��ʾ�ṹ�ĵ�Ԫ������
nodes_str=['Input the number of nodes in the structure:' 10];
nodes=input(nodes_str)      %nodes��ʾ�ṹ�еĽڵ����

e_str=['Input the elastic modulus for each elements, split with comma:' 10];
E=input(e_str,'s');
E=str2num(E)        %EΪ����Ԫ����ģ���ļ��ϡ�
a_str=['Input the area of transverse cross section for each elements, split with comma:' 10];
A=input(a_str,'s');
A=str2num(A)        %AΪ����Ԫ�������ļ��ϡ�

n_n_str=['Input the node sequences for each elements, split with comma, ' 10 'and ' ...
            'separate the elements with semicolon. ' 10 'eg:  2,3;4,5 then press Enter.' 10];
node_number=input(n_n_str,'s');
node_number=str2num(node_number)'   %node_numberΪ����Ԫ��Ӧ�ڵ�ļ��ϣ�ÿ�ж�Ӧһ����Ԫ��

n_c_str=['Input the coordinate of each nodes in sequence, split with comma, ' 10 'and ' ...
            'separate the nodes with semicolon. ' 10 'eg:  0,0,0;0.5,2.5,3.2 then press Enter.' 10];
node_coordinate=input(n_c_str,'s');
node_coordinate=str2num(node_coordinate)'       %node_coordinateΪ���ڵ�����꼯�ϣ�ÿ�ж�Ӧһ��һ���ڵ㡣
node_coordinate=truss3d_ele_coordinate(node_coordinate,node_number);  %node_coordinate����ÿ����Ԫ��Ӧ�ڵ���������󼯺ϡ�

ele_info=zeros(4,eles);
 for i=1:1:eles
     ele_info(:,i)=truss3d_ele_info(node_coordinate(:,i));
 end     %��ɵ�Ԫ������Ϣ���󣬰������������볤�ȡ�

 kk=zeros(6,6,eles);   %��Ԫ�նȾ���ļ��ϡ�
 for i=1:1:eles
     kk(:,:,i)=truss3d_ele_stiffness_mat(E(i),A(i),ele_info(:,i));
 end     %��ɵ�Ԫ�նȾ���
 kk
 
 k=zeros(3*nodes,3*nodes);
 for i=1:1:eles
     k=truss3d_ele_mat_assemble(k,kk(:,:,i),node_number(:,i));
 end     %װ������նȾ���
 k
 
u_u_i_str=['Input the indexs of the unknown nodal displacement in the U vector ' 10 'of the ' ...
    'global stiffness equations K*U=F, split with semicolon.' 10]; 
unknown_u_index=input(u_u_i_str,'s');
unknown_u_index=str2num(unknown_u_index)       %unknown_u_index������նȾ�����δ֪�ڵ�λ����ŵ�������

k_u_str=['Input the known nodal displacement sequence in the U vector ' 10 'of the ' ...
    'global stiffness equations K*U=F, split with semicolon.' 10]; 
known_u=input(k_u_str,'s');
known_u=str2num(known_u)        %known_u������նȷ�������֪�ڵ�λ�Ƶļ��ϡ�


k_f_str=['Input the known nodal force sequence in the F vector ' 10 'of the ' ...
    'global stiffness equations K*U=F, split with semicolon.' 10]; 
known_f=input(k_f_str,'s');
known_f=str2num(known_f)       %known_f������նȷ�������֪�ڵ����ļ��ϡ�

unknown_var=tetra_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
unknown_u=unknown_var{1,1}
unknown_f=unknown_var{1,2}
U=unknown_var{1,3}  %UΪ��λ������
F=unknown_var{1,4}   %FΪ����������

u=zeros(6,eles);    %uΪ����Ԫ��Ӧ�ڵ��λ�����󼯺ϡ�
 for i=1:1:eles
       u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3)];
 end

 force=zeros(eles,1);
 stress=zeros(eles,1);
 for i=1:1:eles
     force(i)=truss3d_ele_inner_force(E(i),A(i),ele_info(:,i),u(:,i));
     stress(i)=force(i)/A(i);
 end
 force
 stress       %�����Ԫ��������Ӧ����


 

%   E= 2e11;
%   A=[0.003;0.003;0.003;0.003];
%   eles=4;
%   nodes=5
%   node_number=[1,2,3,4;5,5,5,5];  %ÿ����Ԫ��Ӧ������յ�Ľڵ�������
%   node_coordinate=[0,0,-3,0,5,0;
%                    -3,0,0,0,5,0;
%                   0,0,3,0,5,0;
%                    4,0,0,0,5,0]';       %ÿ�б�ʾÿ����Ԫ������յ�Ľڵ��������С�
%   ele_info=zeros(4,eles);
%   for i=1:1:eles
%      ele_info(:,i)=truss3d_ele_info(node_coordinate(:,i));
%   end     %��ɵ�Ԫ������Ϣ���󣬰������������볤�ȡ�
%  kk=zeros(6*eles,6);   %��Ԫ�նȾ���ļ��ϡ�
%   for i=1:1:eles
%       kk(6*i-5:6*i,:)=truss3d_ele_stiffness_mat(E,A(i),ele_info(:,i));
%   end     %��ɵ�Ԫ�նȾ���
%   k=zeros(3*nodes,3*nodes);
%   for i=1:1:eles
%       k=truss3d_ele_mat_assemble(k,kk(6*i-5:6*i,:),node_number(:,i));
%   end     %װ������նȾ���
%   m=k(13:15,13:15);
%   f=[15e3;0;-2e4];
%   u=m\f;          %������֪����ɾ������նȾ����Ա�����������е�δ֪����
%   U=[0;0;0;0;0;0;0;0;0;0;0;0;u]        %UΪ����λ�ƾ���
%   F=k*U      %FΪ����������
%   u=zeros(6,eles);
%   for i=1:1:eles
%      u(:,i)=[U((3*node_number(1,i)-2):(3*node_number(1,i)));U((3*node_number(2,i)-2):(3*node_number(2,i)))];
%   end         %ÿ����Ԫ��Ӧ�������ڵ��λ�Ʒ�����
%   force=zeros(eles,1);
%   stress=zeros(eles,1);
%   for i=1:1:eles
%       force(i)=truss3d_ele_inner_force(E,A(i),ele_info(:,i),u(:,i));
%       stress(i)=force(i)/A(i);
%   end
%   force
%   stress       %�����Ԫ��������Ӧ����
% 
% 
% 
% 
% 
% 
% 
