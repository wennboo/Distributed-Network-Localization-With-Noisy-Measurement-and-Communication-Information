%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%验证增益a(t)的gama对仿真影响（仅仅考虑通信噪声，注意到初值很影响迭代收敛，本文假定初值等于真实值
%%%%%%%%%%%%
clc;
clear;
close all;
flag=1;
if(flag)
    format long;
else
    format short;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%
flag_ini=0;%%%%%%%%%先1后0
flag1=0;
flag_ite=1;
flag_fig=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%图形及各点间角度的关系
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A0=[0 0 0 0 0;
   1 0 0 0 0;
   0 1 0 0 0;
   0 0 1 0 0;
   0 0 0 1 0];
flag=1;
if(flag)
A0=zeros(16,16);
A0(2,1)=1;
A0(3,2)=1;
A0(4,3)=1;
A0(5,4)=1;
A0(6,5)=1;

A0(7,2)=1;
A0(8,3)=1;
A0(9,4)=1;
A0(10,5)=1;
A0(11,6)=1;

A0(12,7)=1;
A0(13,8)=1;
A0(14,9)=1;
A0(15,10)=1;
A0(16,11)=1;

end
N0=length(A0);
A=A0+A0';
A=A>0;%%%%%%%%%%%%%%%%%变成对称阵
% A0=A;
% A0(1,2)=0;
% A0(1,6)=0;
L=(diag(sum(A0,2))-A0);
M=16000;%迭代次数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%增益系数
im=sqrt(-1);
gama0=0.85;
gama1=0.7;
gama2=0.55;
c=1/max(sum(A0,2))-0.1;
c=0.8;
%c=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigma0=pi/6;%方位角测量噪声界限
%sigma0=0.017;

sigma0=0;
sigma1=0;%测量噪声方差
sigma2=0.5;%通信噪声方差

flag=0;
if(flag)
sigma0=0;
sigma1=0;%测量噪声方差
sigma2=50;%通信噪声方差 
end



%%%迭代初始值
%%%初始值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%目标配置q0
posi00=12;
theta00=[pi/4];
%yibs=0.01

qs=unifrnd (-8,8,1,N0-1);
qs=qs+im*unifrnd (-8,8,1,N0-1);
ori=unifrnd (-3*pi/4,3*pi/4,1,N0-1);
p_ini=-12*ones(1,N0-1);

%%%%%%%%%%%%%%%%%%%精确的位置及初始值
q0=[posi00 qs].';

% yibs1=[0;yibs*ones(N0-1,1)];
% p0=q0+yibs1;%%%%%%%%假定初值等于真实值
p0=[posi00,p_ini]';
rel_theta0=[theta00,ori]';
star_z=exp(im*rel_theta0);

% hat_theta0=exp(im*rel_theta0)+yibs1%%%%%%%%%%%%%%%%%%%%%%
hat_theta0=zeros(N0,1);
hat_theta0(1)=exp(im*rel_theta0(1));%%%%%%%%%%%%%%%%%复数形式的初值
%%%%%%%%%邻居节点间互相测量
Ang0=zeros(N0);
P_r0=zeros(N0);
Dis0=zeros(N0);
Or_0=zeros(N0);
Q_0=zeros(N0);
%Err=zeros(N0);
for i0=1:N0
    for j0=1:N0
        P_r0(i0,j0)=q0(j0)-q0(i0);
        Ang0(i0,j0)=phase(P_r0(i0,j0));
        Or_0(i0,j0)=rel_theta0(i0)-rel_theta0(j0);
        Dis0(i0,j0)=norm(P_r0(i0,j0));
        Q_0(i0,j0)=q0(j0)-q0(i0);
    end
end
ro0=exp(-im*rel_theta0);%%%%%%
Ro0=-rel_theta0*ones(1,N0);
L_ang0=Ro0+Ang0;%%%%%%%本地坐标下的角度
v10=unifrnd (-sigma0,sigma0,N0,N0);%%%%%%%%平均分布的测量噪声
v20=unifrnd (-sigma1,sigma1,N0,N0);
B_dis0=Dis0+v20;
B_bea0=L_ang0+v10;
B_b0=B_dis0.*exp(im*B_bea0);
B0=Dis0.*exp(im*L_ang0);
%Noi0=zeros(N0);
Ang0=A.*(L_ang0+v10);%%%%%%%%%注意此处为A（A对称）
%%%%%%%%%%mu_ji的求解
mu0=(A.*L_ang0)'-(A.*L_ang0);
B_mu0=Ang0'-Ang0;%%%%%%%%%%%第一次测量误差
T_mu0=B_mu0;%%%波浪mu
H_mu0=B_mu0;%%帽mu
H_v10=zeros(N0,1);
H_mu00=B_mu0;
H00=A0.*exp(-im*H_mu00');%%%%%%%%%%%Q矩阵的计算
%tem0=A0.*exp(-im*Or_0)
H0=diag(sum(A0,2))+H00;

B_b00=B_b0;
B_dis00=B_dis0;
B_bea00=B_bea0;
%H_v100=H_v10;
%B_mu00=B_mu0;
E=eye(N0);
D0=[];%%%%%%%%%%%%方位估计中D
D=[];%%%%%%%%%%定位中D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%邻接矩阵的块对角阵
for i=1:N0;
D0=blkdiag(D0,H00(i,:));
end
size(D0)
for i=1:N0;
D=blkdiag(D,A0(i,:));
end
%%%%%%

%%%%%%%%%%%%%%%%%%%%
%%%%%方位角初始值
t_theta00=hat_theta0;
t_theta10=hat_theta0;
t_theta20=hat_theta0;

%%%%%%%%%%%%方位误差存储
err_o0=norm(t_theta00-star_z)/norm(hat_theta0-star_z);
err_o1=err_o0;
err_o2=err_o0;

%%%%%%%%%%%方位值存储
t_theta_t0=hat_theta0;
t_theta_t1=hat_theta0;
t_theta_t2=hat_theta0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%坐标初始值
p00=p0-q0;
p10=p00;
p20=p00;
%%%%%%%%%%%迭代误差存储
M_p0=norm(p00)/norm(p0-q0);
M_p1=M_p0;
M_p2=M_p0;
%p00(1)=5+im*5;
%%%%%%%%%%%坐标存储
p_t0=p00;
p_t1=p00;
p_t2=p00;

h_v20=unifrnd (-sigma2,sigma2,N0^2,1)+im*unifrnd (-sigma2,sigma2,N0^2,1);;
h_v200=h_v20;

%%%%%%%%%%%%
w20=unifrnd (-sigma2,sigma2,N0^2,1)+im*unifrnd (-sigma2,sigma2,N0^2,1);

if(flag_ite)
for t=1:M
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    at0=c/(1+t)^gama0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %at1=c/(1+t)^gama1;
    %at2=c/(1+t)^gama2;
    t_theta01=(E-at0*H0)*t_theta00-at0*D0*h_v200;
    t_theta_t0=[t_theta_t0 t_theta01];%%%%%%存储估计值
    t_theta00=t_theta01;
    tem_err_o=norm(t_theta01-star_z)/norm(hat_theta0-star_z);
    err_o0=[err_o0 tem_err_o];
    %p1=zeros(N0,1);
    H_w10=[];
    for m1=1:1:N0
        tem1=0;
        for m2=1:1:N0
            tem0=t_theta01(m1)*B_b00(m1,m2);
            tem0=A0(m1,m2)*(Q_0(m1,m2)-tem0);
            tem1=tem0+tem1;
        end
        H_w10=[H_w10 tem1];
    end
    H_w10=H_w10.';
  
    p01=(E-at0*L)*p00+at0*D*w20+at0*H_w10;
    M_p0=[M_p0 norm(p01)/norm(p0-q0)];
    p_t0=[p_t0 p01];%%%%%%%%%%%存储真实位置
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p00=p01; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    at1=c/(1+t)^gama1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %at1=c/(1+t)^gama1;
    %at2=c/(1+t)^gama2;
    t_theta11=(E-at1*H0)*t_theta10-at1*D0*h_v200;
    t_theta_t1=[t_theta_t1 t_theta11];%%%%%%存储估计值
    t_theta10=t_theta11;
    tem_err_o=norm(t_theta11-star_z)/norm(hat_theta0-star_z);
    err_o1=[err_o1 tem_err_o];
    %p1=zeros(N0,1);
    H_w10=[];
    for m1=1:1:N0
        tem1=0;
        for m2=1:1:N0
            tem0=t_theta11(m1)*B_b00(m1,m2);
            tem0=A0(m1,m2)*(Q_0(m1,m2)-tem0);
            tem1=tem0+tem1;
        end
        H_w10=[H_w10 tem1];
    end
    H_w10=H_w10.';
  
    p11=(E-at1*L)*p10+at1*D*w20+at1*H_w10;
    M_p1=[M_p1 norm(p11)/norm(p0-q0)];
    p_t1=[p_t1 p11];%%%%%%%%%%%存储真实位置
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p10=p11; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    at2=c/(1+t)^gama2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %at1=c/(1+t)^gama1;
    %at2=c/(1+t)^gama2;
    t_theta21=(E-at2*H0)*t_theta20-at2*D0*h_v200;
    t_theta_t2=[t_theta_t2 t_theta21];%%%%%%存储估计值
    t_theta20=t_theta21;
    tem_err_o=norm(t_theta21-star_z)/norm(hat_theta0-star_z);
    err_o2=[err_o2 tem_err_o];
    %p1=zeros(N0,1);
    H_w10=[];
    for m1=1:1:N0
        tem1=0;
        for m2=1:1:N0
            tem0=t_theta21(m1)*B_b00(m1,m2);
            tem0=A0(m1,m2)*(Q_0(m1,m2)-tem0);
            tem1=tem0+tem1;
        end
        H_w10=[H_w10 tem1];
    end
    H_w10=H_w10.';
  
    p21=(E-at2*L)*p20+at2*D*w20+at2*H_w10;
    M_p2=[M_p2 norm(p21)/norm(p0-q0)];
    p_t2=[p_t2 p21];%%%%%%%%%%%存储真实位置
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p20=p21; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%不同增益用公共噪声%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v10=unifrnd (-sigma0,sigma0,N0,N0);%%%%%%%%平均分布的测量噪
    v20=unifrnd (-sigma1,sigma1,N0,N0);
    B_dis1=Dis0+v20;
    B_bea1=L_ang0+v10;
    B_dis00=(t-1)/t*B_dis00+1/t*B_dis1;
    B_bea00=(t-1)/t*B_bea00+1/t*B_bea1;
    B_b00=B_dis00.*exp(im*B_bea00);
    Ang1=A.*B_bea1;
%%%%%%%%%%tilde mu_ji的求解
    B_mu1=Ang1'-Ang1;%%%%%%%%%%%方位角bar_mu
    H_mu00=(t-1)/t*H_mu00+1/t*B_mu1;
    
    H00=A0.*exp(-im*H_mu00');%%%%%%%%%%%Q矩阵的计算
    %tem0=A0.*exp(-im*Or_0)
    H0=diag(sum(A0,2))+H00;
    D0=[];%%%%%%%%%%%%方位估计中D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%邻接矩阵的块对角阵
    for i=1:N0;
        D0=blkdiag(D0,H00(i,:));
    end
    %%%%%%hat v1的求解   
    h_v200=unifrnd (-sigma2,sigma2,N0^2,1)+im*unifrnd (-sigma2,sigma2,N0^2,1);%%%%%方位角估计通信噪声
    w20=unifrnd (-sigma2,sigma2,N0^2,1)+im*unifrnd (-sigma2,sigma2,N0^2,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
end



% p_end=p(:,end)+q0;
% t_theta_end=t_theta(:,end)
% phase(t_theta_end)
% phase(star_z)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(flag_fig)
subplot(2,1,1)    
plot(err_o0,'b');
hold on;
plot(err_o1,'g');
hold on;
plot(err_o2,'r');
%legend('big gama', 'middle gama','small gama');


subplot(2,1,2) 
plot(M_p0,'b');  
hold on;
plot(M_p1,'g'); 
hold on;
plot(M_p2,'r'); 
%legend('big gama', 'middle gama','small gama');


    
% if(flag1)
% hold on;
% load err_lee;
% plot(M_p,'r');
% hold on;
% load err_lin;
% plot(M_p,'b');
% end
end    
%axis equal;





