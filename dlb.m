%系统参数
M = 2; m = 0.1; l = 0.5;
I = 0;
b1 = 0.2;b2 = 0.1;
g = 9.8;

%状态方程
A = [0 1 0 0
     0 (-(I+m*l^2)*b1)/(I*(M+m)+M*m*l^2) -(m^2*g*l^2)/(I*(M+m)+M*m*l^2) (m*l*b2)/(I*(M+m)+M*m*l^2)
     0 0 0 1
     0 (m*l*b1)/(I*(M+m)+M*m*l^2) (m*g*l*(m+M))/(I*(M+m)+M*m*l^2) (-(M+m)*b2)/(I*(M+m)+M*m*l^2)];
B = [0; (I+(m*l)^2)/(I*(M+m)+M*m*l^2); 0; (-m*l)/(I*(M+m)+M*m*l^2)];

%状态反馈
alpha = poly(A);
alpha_s = [1 24 196 720 1600];
K_ = fliplr(alpha_s - alpha);
K_ = K_(1:4);
Q = [B A*B A*A*B A*A*A*B];
rQ = rank(Q);%4,完全可控
P = Q*[alpha(4) alpha(3) alpha(2) 1;alpha(3) alpha(2) 1 0;alpha(2) 1 0 0;1 0 0 0];
K = K_/P;

%状态观测器
Ac = [0 0 1 0
      0 0 0 1
      A(2,1) A(2,3) A(2,2) A(2,4)
      A(4,1) A(4,3) A(4,2) A(4,4)];
Bc = [0; 0; B(2); B(4)];
Cc = [1 0 0 0; 0 1 0 0];
Qc = [Cc;Cc*A;Cc*A*A;Cc*A*A*A];
rQc = rank(Qc);%4,完全可观
C1 = Cc(:,1:2); C2 = Cc(:,3:4);
P_ = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
A_ = inv(P_)*Ac*P_;
B_ = inv(P_)*Bc;
C_ = Cc*P_;