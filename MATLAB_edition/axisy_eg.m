%��Գƾ���Ԫ������
source=xlsread('axisy_eg_input.xlsx');
source(all(isnan(source),2),:)=[];
eles=source(1,1);
nodes=source(1,2);
e=source(2,1);
nu=source(2,2);
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
%   ȫ������r��[Ni]*[r1;r2;r3;r4]��ʾ��xi����֪��������ֵ��r\z����s\t��ʾ�ˡ�
%   ����ֻ��һ����Ԫ����������������ȡ�ľ���Ԫ������Ԫ�ߴ���ȫһ��������Ԫ��Ӧ��[Ni]��һ���ġ�
%   ���������������Ԫ�Ľڵ�ֵֻҪ���Ǻ;���Ԫһ����Ȼ���갴�ĸ���1ȡ��[Ni]Ҳ��һ���ġ�


syms s t u1 u2 u3 u4
r_expression=sym(zeros(eles,1));     %x�������������洢����Ԫ�ĳ�����x��s\t��ʾ�ı��ʽ��y���ž����Ӧ��
z_expression=sym(zeros(eles,1));


st_coordinate=sym([1,1;-1,1;-1,-1;1,-1]');       %st_coordinate���������ʾ����Ȼ�����£�ÿ����Ԫ��Ӧ�Ľڵ�����ꡣ
u=sym([u1;u2;u3;u4]);     %u�Է��ű�ʾ��Ԫ�ĸ��ڵ��uλ�ƣ�ֻ��Ϊ�����ڼ������������ֵ������û��ʵ�ʶ�Ӧ����ֵ֪��

st_ele_mat=sym(zeros(4,4));         %st_ele_mat����A*[a1;a2;a3;a4]=[u1;u2;u3;u4]�������е�A��ÿ����[1,si,ti,si*ti]��
for i=1:4
    si=st_coordinate(1,i);
    ti=st_coordinate(2,i);
    st_ele_mat(i,:)=[1,si,ti,si*ti];
end

 
ai=st_ele_mat\u;        %aiΪ[a1;a2;a3;a4]����
u_expression=ai'*sym([1;s;t;s*t]);       %u_expression����s\t��ʾ�ĳ�����u������
Ni=rect_coeffs_ui(u_expression,u);  %����Ni��������ui��ϵ������ȡ��

r_ele=sym(zeros(4,eles));
z_ele=sym(zeros(4,eles));   %r_ele�Ǹ���Ԫ�ڵ����������r�ļ��ϣ�����֪����ֵ��
for i=1:eles
    r_ele(:,i)=node_coordinate([1,3,5,7],i);
    z_ele(:,i)=node_coordinate([2,4,6,8],i);
end


for i=1:eles
    r_expression(i)=Ni'*r_ele(:,i);
    z_expression(i)=Ni'*z_ele(:,i);
end     %r��z��s\t���ʽ�Ѿ��õ���


b=sym(zeros(4,8,eles));     %bΪ����ԪB����ļ��ϡ�
for i=1:eles
    for j=1:4
        b(:,2*j-1:2*j,i)=axisy_ele_b_1(Ni(j),r_expression(i),z_expression(i));
    end
end

 d=e*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),nu/(1-nu),0;
                               nu/(1-nu),1,nu/(1-nu),0;
                               nu/(1-nu),nu/(1-nu),1,0;
                               0,0,0,(1-2*nu)/(2-2*nu)];
 d=sym(d);
 
 kk=sym(zeros(8,8,eles));    %kk�Ǹ���Ԫ�նȾ���ļ��ϡ�
 for i=1:1:eles
      kk(:,:,i)=axisy_ele_stiff_mat_1(b(:,:,i),d,r_expression(i),z_expression(i));
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
 
 stress=sym(zeros(4,eles));   %sigmaΪ��ԪӦ�����󣬴�ʱ������ǹ�������ı����������ǳ�����strainΪӦ�䡣
 strain=sym(zeros(4,eles));
 for i=1:1:eles
     stress(:,i)=d*b(:,:,i)*u(:,i);
     strain(:,i)=b(:,:,i)*u(:,i);
 end
 num_stress=vpa(stress,6);
 num_strain=vpa(strain,6);
 
stress_zero_point=double(subs(stress,{s,t},{0,0}));
strain_zero_point=double(subs(strain,{s,t},{0,0}));

%����д�ļ���
f=fopen('axisy_eg_output.xls','w');

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


fclose(f);

  








% e=200e9;
% nu=0.3;
% eles=4;
% nodes=9;
% node_number=[6,1,5,7;
%     7,5,2,8;
%     9,7,8,3;
%     4,6,7,9]';  %ÿ����Ԫ��Ӧ�ڵ�ı�ţ�����ʱ�롢һ�����������С�
% node_coordinate=0.001*[40,10;40,0;60,0;60,10;
%     40,5;50,10;50,5;50,0;60,5]';    %ÿ���ڵ�����ꡣ
% 
% syms u1 v1 u2 v2 u3 v3 u4 v4 u5 v5 u6 v6 u7 v7 u8 v8 u9 v9 r z
% u_sym=sym([u1;u2;u3;u4;u5;u6;u7;u8;u9]);
% v_sym=sym([v1;v2;v3;v4;v5;v6;v7;v8;v9]);
% 
% node_info_for_n=zeros(nodes,4);   %node_info_for_n�Ǹ��ڵ�[1;ri;zi;ri*zi]�ļ��ϣ���������״��������Nʱ�á�
% for i=1:1:nodes
%     node_info_for_n(i,:)=axisy_node_info_for_n(node_coordinate(:,i));
% end
% 
% ele_info_for_n=zeros(4*eles,4);     %ele_info_for_n�Ǹ���Ԫ��Ӧ�ڵ�[1;ri;zi;ri*zi]�ļ��ϣ���������״��������Nʱ�á�
% for i=1:1:eles
%     ele_info_for_n(4*i-3:4*i,:)=axisy_ele_info_for_n(node_info_for_n,node_number(:,i));
% end
% 
% ele_u_sym=sym(zeros(4,eles));  %ele_u_symΪ����ԪΪ���N�к���ϵ�����õĸ��ڵ����꣬sym�࣬������֪����
% ele_v_sym=sym(zeros(4,eles));
% for i=1:1:eles
%     ele_u_sym(:,i)=axisy_ele_uv_sym(u_sym,node_number(:,i));
%     ele_v_sym(:,i)=axisy_ele_uv_sym(v_sym,node_number(:,i));
% end
% 
% n_coeffs_for_rz=sym(zeros(8,eles));     %n_coeffs_for_rz�Ǹ���Ԫ��������λ�Ʋ�ֵ��ʽ��r\zΪ����õ���ϵ����ǰ���ж�Ӧuλ�ƣ������ж�Ӧrλ�ơ�
% for i=1:1:eles
%     n_coeffs_for_rz(1:4,i)=axisy_n_coeffs_for_rz(ele_u_sym(:,i),ele_info_for_n(4*i-3:4*i,:));
%     n_coeffs_for_rz(5:8,i)=axisy_n_coeffs_for_rz(ele_v_sym(:,i),ele_info_for_n(4*i-3:4*i,:));
% end
% 
% uv_functions=sym(zeros(eles,1));  %uv_functions�����յõ��ĵ�Ԫ������һ��λ�Ƶĵĺ�����
% for i=1:1:eles
%     uv_functions(i)=sym([1,r,z,r*z])*n_coeffs_for_rz(1:4,i);
% end
% uv_functions
% 
% n=sym(zeros(eles,4));   %nΪ����Ԫ��Ӧ��λ�ƺ�������ļ��ϣ���uv_functions��ui\vi����õ���ϵ��������ֵ������
% for i=1:1:eles
%     p=ele_u_sym(:,i);
%     uv_functions(i)=collect(uv_functions(i),p);
%     p=p(end:-1:1);
%     n(i,:)=coeffs(uv_functions(i),p);
% end
% n
% 
% jacob_eles=sym(zeros(4,2*eles));    %���ɸ���Ԫ��Ӧ�Ĳ�ֵ��������ſɱȾ��󼯺ϡ�
%  for i=1:1:eles
%      jacob_eles(:,2*i-1:2*i)=jacobian(n(i,:));
%  end
%  jacob_eles
%  
%  b=sym(zeros(4*eles,8));     %bΪ����ԪB����ļ��ϡ�
%  for i=1:1:eles
%      b(4*i-3:4*i,:)=axisy_ele_b(jacob_eles(:,2*i-1:2*i),n(i,:));
%  end     %����B���󼯺ϡ�
%  b
%  
%  d=e*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),nu/(1-nu),0;
%                                nu/(1-nu),1,nu/(1-nu),0;
%                                nu/(1-nu),nu/(1-nu),1,0;
%                                0,0,0,(1-2*nu)/(2-2*nu)];
%  d=sym(d);
%  
%  kk=sym(zeros(8*eles,8));    %kk�Ǹ���Ԫ�նȾ���ļ��ϡ�
%  for i=1:1:eles
%       kk(8*i-7:8*i,:)=axisy_ele_stiff_mat(node_number(:,i),node_coordinate,b(4*i-3:4*i,:),d);
%  end
%  kk
%  kk=double(kk)
%  
%  k=zeros(2*nodes,2*nodes);   %kΪ����նȾ���
%  for i=1:1:eles
%      k=rect_assembly_mat(k,kk(8*i-7:8*i,:),node_number(:,i));
%  end
% k
% 
%  unknown_u_index=[1;3;9;11;13;15];   %unknown_u_index��δ֪λ����ŵ�������
%  known_f=[2513/2;2513/2;2513;0;0;0];
%  known_u=zeros(12,1);
%  unknown_var=tri_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
%  unknown_u=unknown_var{1,1}
%  unknown_f=unknown_var{1,2}
%  U=unknown_var{1,3}  %UΪ��λ������
%  F=unknown_var{1,4}   %FΪ����������
%   
%  u=zeros(8,eles);    %uΪ����Ԫ��Ӧ�ڵ��λ�����󼯺ϡ�
%  for i=1:1:eles
%        u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2);U(node_number(3,i)*2-1:node_number(3,i)*2);
%            U(node_number(4,i)*2-1:node_number(4,i)*2)];
%  end
%  
%  stress=sym(zeros(4,eles));   %sigmaΪ��ԪӦ�����󣬴�ʱ������ǹ�������ı����������ǳ�����strainΪӦ�䡣
%  strain=sym(zeros(4,eles));
%  for i=1:1:eles
%      stress(:,i)=d*b(4*i-3:4*i,:)*u(:,i);
%      strain(:,i)=b(4*i-3:4*i,:)*u(:,i);
%  end
%  num_stress=vpa(stress,6)
%  num_strain=vpa(strain,6)
%  
%  eles_center=zeros(2,eles);  %ȷ��ÿ����Ԫ���ĵ����ꡣ
%   for i=1:1:eles
%      eles_center(:,i)=[(node_coordinate(1,node_number(1,i))+node_coordinate(1,node_number(3,i)))/2;
%          (node_coordinate(2,node_number(1,i))+node_coordinate(2,node_number(3,i)))/2];
%   end
%   eles_center
%   
%   stress_zero_point=sym(zeros(4,eles));
%   strain_zero_point=sym(zeros(4,eles)); %������Ԫ���ĵ��Ӧ��Ӧ�䡣
%   for i=1:1:eles
%       for j=1:1:4
%       if ~isempty(findstr('r',findsym(stress(j,i)))) & ~isempty(findstr('z',findsym(stress(j,i))))
%           stress_zero_point(j,i)=subs(stress(j,i),{r,z},{eles_center(1,i),eles_center(2,i)});
%       elseif isempty(findstr('r',findsym(stress(j,i)))) & isempty(findstr('z',findsym(stress(j,i))))
%           stress_zero_point(j,i)=stress(j,i);
%       elseif isempty(findstr('z',findsym(stress(j,i))))
%           stress_zero_point(j,i)=subs(stress(j,i),r,eles_center(1,i));
%       elseif isempty(findstr('r',findsym(stress(j,i))))
%           stress_zero_point(j,i)=subs(stress(j,i),z,eles_center(2,i));
%       end
%       end
%   end
%   for i=1:1:eles
%       for j=1:1:4
%       if ~isempty(findstr('r',findsym(strain(j,i)))) & ~isempty(findstr('z',findsym(strain(j,i))))
%           strain_zero_point(j,i)=subs(strain(j,i),{r,z},{eles_center(1,i),eles_center(2,i)});
%       elseif isempty(findstr('r',findsym(strain(j,i)))) & isempty(findstr('z',findsym(strain(j,i))))
%           strain_zero_point(j,i)=strain(j,i);
%       elseif isempty(findstr('z',findsym(strain(j,i))))
%           strain_zero_point(j,i)=subs(strain(j,i),r,eles_center(1,i));
%       elseif isempty(findstr('r',findsym(strain(j,i))))
%           strain_zero_point(j,i)=subs(strain(j,i),z,eles_center(2,i));
%       end
%       end
%   end
%   
%    stress_zero_point=double(stress_zero_point)       
%    strain_zero_point=double(strain_zero_point)
