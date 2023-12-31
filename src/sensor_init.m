
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
delta=0.7;
c=1/max(sum(A0,2))-0.1;
c=0.4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigma0=pi/6;%方位角测量噪声界限
%sigma0=0.017;

sigma0=pi/6;
sigma1=0.1;%测量噪声方差
sigma2=0.1;%通信噪声方差

flag=0;
if(flag)
sigma0=0;
sigma1=0;%测量噪声方差
sigma2=0;%通信噪声方差 
end



%%%迭代初始值
%%%初始值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%目标配置q0
posi00=12;
theta00=[pi/4];



if(flag_ini)
qs=unifrnd (-8,8,1,N0-1);
qs=qs+im*unifrnd (-8,8,1,N0-1);
ori=unifrnd (-3*pi/4,3*pi/4,1,N0-1);
p_ini=-12*ones(1,N0-1);
save ini_q0 qs;
save ini_pin p_ini;
save ini_ori ori;
end