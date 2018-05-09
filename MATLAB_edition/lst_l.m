function q = lst_l( node_coordinate )
%LST_L Summary of this function goes here
%   Detailed explanation goes here
%   ���ɸ���Ԫ��l���꣬n_cΪ��Ԫ�˵����ꡣ
syms x y
n_c=sym(node_coordinate);
q=sym(zeros(3,1));
aa=sym(zeros(3,1));
aa(1)=det([1 x y;1 n_c(1,2) n_c(2,2);1 n_c(1,3) n_c(2,3)]);
aa(2)=det([1 x y;1 n_c(1,3) n_c(2,3);1 n_c(1,1) n_c(2,1)]);
aa(3)=det([1 x y;1 n_c(1,1) n_c(2,1);1 n_c(1,2) n_c(2,2)]);
a=det([1 n_c(1,1) n_c(2,1);1 n_c(1,2) n_c(2,2);1 n_c(1,3) n_c(2,3)]);
for i=1:3
    q(i)=aa(i)/a;
end
end

