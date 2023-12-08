%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%û�г�ʼ��λԼ���Ķ�λ�������Լ���λ�㷨����֪�ڵ��λ�úͳ�ʼ��λ��sensor_initʵ��
%%%%%%%%%%%%��֤���еĶ��й����ĳ�ʼֵ
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
flag_ini=0;%%%%%%%%%��1��0
flag1=1;
flag_ite=1;
flag_fig=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ͼ�μ������ǶȵĹ�ϵ
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
A=A>0;%%%%%%%%%%%%%%%%%��ɶԳ���
% A0=A;
% A0(1,2)=0;
% A0(1,6)=0;
L=(diag(sum(A0,2))-A0);
M=16000;%��������

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%����ϵ��
im=sqrt(-1);
delta=0.7;
c=1/max(sum(A0,2))-0.1;
c=0.4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigma0=pi/6;%��λ�ǲ�����������
%sigma0=0.017;

sigma0=pi/6;
sigma1=0.2;%������������
sigma2=0.2;%ͨ����������

flag=0;
if(flag)
sigma0=0;
sigma1=0;%������������
sigma2=0;%ͨ���������� 
end



%%%������ʼֵ
%%%��ʼֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Ŀ������q0
posi00=12;
theta00=[pi/4];


load ini_q0
load ini_pin
load ini_ori
%%%%%%%%%%%%%%%%%%%��ȷ��λ�ü���ʼֵ
q0=[posi00 qs].';
p0=[posi00,p_ini]';
rel_theta0=[theta00,ori]';
star_z=exp(im*rel_theta0);
hat_theta0=zeros(N0,1);
hat_theta0(1)=exp(im*rel_theta0(1));%%%%%%%%%%%%%%%%%������ʽ�ĳ�ֵ
%%%%%%%%%�ھӽڵ�以�����
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
L_ang0=Ro0+Ang0;%%%%%%%���������µĽǶ�
v10=unifrnd (-sigma0,sigma0,N0,N0);%%%%%%%%ƽ���ֲ��Ĳ�������
v20=unifrnd (-sigma1,sigma1,N0,N0);
B_dis0=Dis0+v20;
B_bea0=L_ang0+v10;
B_b0=B_dis0.*exp(im*B_bea0);
B0=Dis0.*exp(im*L_ang0);
%Noi0=zeros(N0);
Ang0=A.*(L_ang0+v10);%%%%%%%%%ע��˴�ΪA��A�Գƣ�
%%%%%%%%%%mu_ji�����
mu0=(A.*L_ang0)'-(A.*L_ang0);
B_mu0=Ang0'-Ang0;%%%%%%%%%%%��һ�β������
T_mu0=B_mu0;%%%����mu
H_mu0=B_mu0;%%ñmu
H_v10=zeros(N0,1);
H_mu00=B_mu0;
H00=A0.*exp(-im*H_mu00');%%%%%%%%%%%Q����ļ���
%tem0=A0.*exp(-im*Or_0)
H0=diag(sum(A0,2))+H00;
% for i0=1:N0
%     for j0=1:N0
%         tem0=mod(H_mu00(i0,j0)+2*pi,2*pi)-pi;
%         tem1=A0(i0,j0)*(Or_0(i0,j0)-tem0);
%         H_v10(i0)=H_v10(i0)+tem1;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�������
%t_theta0=hat_theta0-rel_theta0;
B_b00=B_b0;
B_dis00=B_dis0;
B_bea00=B_bea0;
%H_v100=H_v10;
%B_mu00=B_mu0;
E=eye(N0);
D0=[];%%%%%%%%%%%%��λ������D
D=[];%%%%%%%%%%��λ��D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%�ڽӾ���Ŀ�Խ���
for i=1:N0;
D0=blkdiag(D0,H00(i,:));
end
size(D0)
for i=1:N0;
D=blkdiag(D,A0(i,:));
end
%%%%%%

%%%%%%%%%%%%%%%%%%%%
%%%%%��ʼֵ
t_theta00=hat_theta0;
t_theta=hat_theta0;
p00=p0-q0;
M_p=norm(p00)/norm(p0-q0);
%p00(1)=5+im*5;
p=p00;
h_v20=unifrnd (-sigma2,sigma2,N0^2,1)+im*unifrnd (-sigma2,sigma2,N0^2,1);;
h_v200=h_v20;
err_o=norm(t_theta00-star_z)/norm(hat_theta0-star_z);
%%%%%%%%%%%%


if(flag_ite)
for t=1:M
    %%%%%%%%%%%%%%%%%�Ƕ�%%%%%%%%%%%%
    at=c/(1+t)^delta;
    t_theta1=(E-at*H0)*t_theta00-at*D0*h_v200;
    t_theta=[t_theta t_theta1];%%%%%%�洢����ֵ
    tem_err_o=norm(t_theta1-star_z)/norm(hat_theta0-star_z);
    err_o=[err_o tem_err_o];
    p1=zeros(N0,1);
    H_w10=[];
    for m1=1:1:N0
        tem1=0;
        for m2=1:1:N0
            tem0=t_theta00(m1)*B_b00(m1,m2);
            tem0=A0(m1,m2)*(Q_0(m1,m2)-tem0);
            tem1=tem0+tem1;
        end
        H_w10=[H_w10 tem1];
    end
    H_w10=H_w10.';
    w20=unifrnd (-sigma2,sigma2,N0^2,1)+im*unifrnd (-sigma2,sigma2,N0^2,1);
    p1=(E-at*L)*p00+at*D*w20+at*H_w10;
    M_p=[M_p norm(p1)/norm(p0-q0)];
    p=[p p1];%%%%%%%%%%%�洢��ʵλ��
    v10=unifrnd (-sigma0,sigma0,N0,N0);%%%%%%%%ƽ���ֲ��Ĳ�����
    v20=unifrnd (-sigma1,sigma1,N0,N0);
    B_dis1=Dis0+v20;
    B_bea1=L_ang0+v10;
    B_dis00=(t-1)/t*B_dis00+1/t*B_dis1;
    B_bea00=(t-1)/t*B_bea00+1/t*B_bea1;
    B_b00=B_dis00.*exp(im*B_bea00);
    Ang1=A.*B_bea1;
%%%%%%%%%%tilde mu_ji�����
    B_mu1=Ang1'-Ang1;%%%%%%%%%%%��λ��bar_mu
    H_mu00=(t-1)/t*H_mu00+1/t*B_mu1;
    
    H00=A0.*exp(-im*H_mu00');%%%%%%%%%%%Q����ļ���
    %tem0=A0.*exp(-im*Or_0)
    H0=diag(sum(A0,2))+H00;
    D0=[];%%%%%%%%%%%%��λ������D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%�ڽӾ���Ŀ�Խ���
    for i=1:N0;
        D0=blkdiag(D0,H00(i,:));
    end
    %%%%%%hat v1�����   
    h_v200=unifrnd (-sigma2,sigma2,N0^2,1)+im*unifrnd (-sigma2,sigma2,N0^2,1);%%%%%��λ�ǹ���ͨ������
    t_theta00=t_theta1;
    p00=p1; 
end



p_end=p(:,end)+q0;
t_theta_end=t_theta(:,end)
phase(t_theta_end)
phase(star_z)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(flag_fig)
z0=real(q0);
z1=imag(q0);
%%%%%%%%%��ȷ��λ��
scatter(z0(1),z1(1),'s','r');
hold on;
z0=z0(2:N0);
z1=z1(2:N0);
scatter(z0,z1,'o','r');
    hold on;
    z0=real(p0);
    z0=z0(2:N0);
    z1=imag(p0);
    z1=z1(2:N0);
    scatter(z0,z1,'^','r');
    hold on;
    z0=real(p_end);
    z0=z0(2:N0);
    z1=imag(p_end);
    z1=z1(2:N0);
%%%%%%%%%���������λ��
    scatter(z0,z1,'*','b');
    legend('The Anchor node', 'The accurate position of the sensor node','The initial position estimate of the sensor node','The final position estimate of the sensor node');
   %axis equal; 

hold on;
N=length(p);
p1=real(p)+real(q0)*ones(1,N);
p2=imag(p)+imag(q0)*ones(1,N);
for i=2:N0
    plot(p1(i,:),p2(i,:),'k');
    hold on;
end
    axis([-25,15,-15,15]);
    %hold on;
    
if(flag1)
subplot(2,1,1);
plot(err_o,'k');
save un_err_o4 err_o;
subplot(2,1,2);
plot(M_p,'k');      
save un_M_p4 M_p;

end

% subplot(2,1,1);
% plot(err_o,'k');
% hold on;
% load gauss_err_o1;
% plot(err_o,'r');
% hold on;
% load gauss_err_o4;
% plot(err_o,'b');
% subplot(2,1,2);
% plot(M_p,'k');
% hold on;
% load gauss_M_p1;
% plot(M_p,'r');
% hold on;
% load gauss_M_p4;
% plot(M_p,'b');
  end    
%axis equal;





