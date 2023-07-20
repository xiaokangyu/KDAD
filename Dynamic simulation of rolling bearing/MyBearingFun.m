
function dyx = MyBearingFun(t,dx)
% function ZHAI(tm,dt)
%tm为仿真时间，dt为仿真时间步长
tm=1.0;
dt=0.00001;
m_rp=32.1;
m_bR =20;
m_bL =20;
m_rR=4;
m_rL=4;
m_wR = 2.0;
m_wL =2.0;
mc = 50;

M = diag([m_rp,m_rp,m_bR,m_bR,m_bL,m_bL,m_rR,m_rR,m_rL,m_rL,m_wR,m_wR,m_wL,m_wL,mc,mc]);

k = 2.5e7;
kfRH = 7.5e8;
ktRH = 2.5e8;
kfLH = 7.5e8;
ktLH = 2.5e8;
kcH = 2.5e9;
K =zeros(16,16);
K(1,1)=2*k;
K(1,7)=-k;
K(1,9)=-k;

K(2,2)=2*k;
K(2,8)=-k;
K(2,10)=-k;

K(3,3)=kfRH+ktRH;
K(3,15)=-kfRH;
K(3,11)=-ktRH;

K(4,4)=kfRH+ktRH;
K(4,16)=-kfRH;
K(4,12)=-ktRH;

K(5,5)=kfLH+ktLH;
K(5,15)=-kfLH;
K(5,13)=-ktLH;
K(6,6)=kfLH+ktLH;
K(6,16)=-kfLH;
K(6,14)=-ktLH;

K(7,7)=k;
K(7,1)=-k;
K(8,8)=k;
K(8,2)=-k;
K(9,9)=k;
K(9,1)=-k;
K(10,10)=k;
K(10,2)=-k;
K(11,11)=ktRH;
K(11,3)=-ktRH;
K(12,12)=ktRH;
K(12,4)=-ktRH;
K(13,13)=ktLH;
K(13,5)=-ktLH;
K(14,14)=ktLH;
K(14,6)=-ktLH;

K(15,15)=kcH+kfRH+kfLH;
K(15,3)=-kfRH;
K(15,5)=-kfLH;
K(16,16)=kcH+kfRH+kfLH;
K(16,4)=-kfRH;
K(16,6)=-kfLH;

c=2100;
cfRH = 2100;
ctRH = 1050;
cfLH = 2100;
ctLH = 1050;
crb=1050;
ccH = 2100;
C = zeros(16,16);
C(1,1)=c;
C(2,2)=c;
C(3,3)=cfRH+ctRH;
C(3,11)=-ctRH;
C(3,15)=-cfRH;

C(4,4)=cfRH+ctRH;
C(4,12)=-ctRH;
C(4,16)=-cfRH;

C(5,5)=cfLH+ctLH;
C(5,13)=-ctLH;
C(5,15)=-cfLH;
C(6,6)=cfLH+ctLH;
C(6,14)=-ctLH;
C(6,16)=-cfLH;

C(7,7)=crb;
C(8,8)=crb;
C(9,9)=crb;
C(10,10)=crb;

C(11,11)=ctRH;
C(11,3)=-ctRH;
C(12,12)=ctRH;
C(12,4)=-ctRH;

C(13,13)=ctLH;
C(13,5)=-ctLH;
C(14,14)=ctLH;
C(14,6)=-ctLH;

C(15,15)=ccH+cfRH+cfLH;
C(15,3)=-cfRH;
C(15,5)=-cfLH;
C(16,16)=ccH+cfRH+cfLH;
C(16,4)=-cfRH;
C(16,6)=-cfLH;


e=0.05e-3;%偏心距
casing_width=0.01e-3;%机匣间隙
w=200;%角速度
r_0=4e-5;%轴承间隙
g=9.8;%重力加速度
f=0.3;%摩擦系数

C_b=13.34e9;%轴承接触刚度
N_b=8; %轴承滚珠个数
r=40.1e-3;%轴承内圈
R=63.9e-3;%轴承外圈




%-----------------------------------------------------------------------
x( : ,1)=dx(1:16);    %初始位移
v( : ,1)=dx(17:32);    %初始速度

%计算初始力
t0 = 0;
xrR0= 0;
xwR0= 0;
r_d = 0;
[fxbR0,fybR0,xxx] = ComBearForce01(w,r,R,N_b,t0,C_b,0,0,r_0,1);
[fxbL0,fybL0,xxx] = ComBearForce01(w,r,R,N_b,t0,C_b,0,0,r_0,r_d);


f0(1)=m_rp*e*w^2*cos(w*t0);
f0(2)=m_rp*e*w^2*sin(w*t0)-m_rp*g;
f0(3)=0;
f0(4)=-m_bR*g;
f0(5)=0;
f0(6)=-m_bL*g;
f0(7)=fxbR0;
f0(8)=-m_rR*g+fybR0;
f0(9)=fxbL0;
f0(10)=-m_rL*g+fybL0;
f0(11)=-fxbR0;
f0(12)=-m_wR*g-fybR0;
f0(13)=-fxbL0;
f0(14)=-m_wL*g-fybL0;
f0(15)=0;
f0(16)=-mc*g;

    t0 = t;
    x0= x(7,end)-x(11,end);
    y0= x(8,end)-x(12,end);
    x10= x(9,end)-x(13,end);
    y10= x(10,end)-x(14,end);
    [fxbR0,fybR0,r_d] = ComBearForce01(w,r,R,N_b,t0,C_b,x0,y0,r_0,1);
    %res(js) = r_d;
    [fxbL0,fybL0,r_d] = ComBearForce01(w,r,R,N_b,t0,C_b,x10,y10,r_0,0);

    f0(1)=m_rp*e*w^2*cos(w*t0);
    f0(2)=m_rp*e*w^2*sin(w*t0)-m_rp*g;
    f0(3)=0;
    f0(4)=-m_bR*g;
    f0(5)=0;
    f0(6)=-m_bL*g;
    f0(7)=fxbR0;
    f0(8)=-m_rR*g+fybR0;
    f0(9)=fxbL0;
    f0(10)=-m_rL*g+fybL0;
    f0(11)=-fxbR0;
    f0(12)=-m_wR*g-fybR0;
    f0(13)=-fxbL0;
    f0(14)=-m_wL*g-fybL0;
    f0(15)=0;
    f0(16)=-mc*g;
    
    a( : ,1) = inv(M)*(f0'-K*x(:,1)-C*v( : ,1));
    dyx = [v(:,1);a(:,1)];
    
 end   
    

