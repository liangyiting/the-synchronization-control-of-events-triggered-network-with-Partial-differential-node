function []=erreerer();
close all
lw=2;%w的维度
global D C;%w''+cw=0  w'=w' w''=-cw； D参数较大时容易实现同步
D=[0.1,0;0,0.2];C=[0.1,1;-0.5,0.2];%模型线性常数%Rk=max(abs(xl),abs(xr));%
global ffun gfun;
L=2;sigma=2;
ffun=@(w)L/2*(w+sin(1*w));    %定义f函数
gfun=@(w)sigma/2*(w+sin(5*w)+0*w./sqrt(w.^2+1));%定义g函数
%时间
tf=0.2;%终点时间

N=6;%节点个数
global beta;%事件触发器参数， beta越小越稳定
beta=0.1*ones(N,1);
%定义邻接矩阵 随机生成的
%邻接越紧密越容易同步.所以你这里面应该对邻接有要求，比如满足某种条件的小世界网络之类的。
if 1
    %随机生成邻接矩阵
    adj=rand(N,N);adj=adj+adj';adj=adj/2;adj=[adj>0.5];
else %全连接邻接矩阵
    adj=ones(N,N);
end
for i=1:N;adj(i,i)=0;end
%生成连接权值,%连接权值越大，则越容易同步
if 0
    %随机生成连接权值 随机权值有可能出现无法同步的现象
    B=10*adj.*(1+rand(N,N));%相互影响的系数
else
    %固定连接权值
    B=5*adj;
end
for i=1:N;
    B(i,i)=-sum(B(i,:));
end;

%控制器参数
delta=15;%控制器参数
mu=5;%控制器参数
%x网格
xl=0;%左边界
xr=1;%右边界
nx=20;%nx是x分段数-- nx不能太大，否则会造成数值不稳定（偏微分方程求解中，dt/dx必须小于某一个值才会稳定)
dx=(xr-xl)/nx;%步长
maxdt=dx*0.1;%偏微分方程求解中，dt/dx必须小于某一个值才会稳定)

%初始状态w 平衡状态wzero  这里我为了简化，让wzero全部为0了。这个在设置f函数和g函数的时候必须要保证零点是平衡点，因此f被设置成了L*(w+sin(w)),g被设置成了sigma*(w+sin(w))。
%w=ones(N,nx);%w(i*+j) 随机初始化
w=zeros(lw,nx,N);
wzero=0*w;
for j=1:lw;
    %w(j,:,:)=0.5*rand(N,1)*[sin(pi*((1:nx)/nx))];%保证两端为0
    w(j,:,:)=0.15*[sin(pi*((1:nx)/nx))]'*[(1:N)-N/2];%保证两端为0
    wzero(j,:,:)=0*[sin(pi*((1:nx)/nx))]'*rand(1,N);%zeros(N,nx);%零点，这里由于f是线性无偏的，且两端为0，因此平衡点直接为0
end

global wold;%wold是最近一次触发时间点的状态
wold=w;
%%
%求解偏微分方程组
%定义仿真模型
odefun=@(t,z)fun(t,z,adj,B,delta,mu,dx,wzero,lw);
if 0
    opt=odeset('maxstep',maxdt);
    [tt,www]=ode15s(odefun,[0,tf],w(:),opt);
else
    dt=maxdt/10;
    Nt=ceil(tf/dt);
    tt=0:dt:Nt*dt;
    www=zeros(Nt+1,lw*nx*N);
    www(1,:)=w(:);
    z=w(:);
    ckxt=www*0;
    cu=www*0;
    for i=1:Nt;
        t=dt*i;
        t,
        tt(i+1)=t;
        [dz,kxtall,u]=odefun(t,z);
        z=z+dt*dz;
        ckxt(i+1,:)=kxtall(:)';
        www(i+1,:)=z(:)';
        cu(i+1,:)=u(:)';
    end
end
%%
if 1
    ttnew=linspace(0,tf,150);
    wwwnew=interp1(tt,www,ttnew);
    ckxtnew=interp1(tt,ckxt,ttnew);
    cunew=interp1(tt,cu,ttnew);
    tt=ttnew;www=wwwnew;
    ckxt=ckxtnew;cu=cunew.*(ckxt>0);
end
%%
%作图部分
ind=reshape(1:lw*nx*N,lw,nx,N);
for kss=1:lw;%只取第一个维度的
    %ind1=ind(kss,:,:);ind1=ind1(:)';ww=www(:,ind1);%
    cww=cellmat(N,1);
    m=ceil(sqrt(N));
    [xmesh,tmesh]=meshgrid(xl+dx:dx:xr,tt);
    for i=1:N;
        i1=ind(kss,:,i);
        %cww{i}=ww(:,i:N:end);
        cww{i}=www(:,i1);
        %figure(1);subplot(li1,li2,i*2-1);pcolor(cww{i});xlabel('x轴');ylabel('t轴');title(strcat('节点',num2str(i)));shading interp;
        figure(kss+lw);
        subplot(m,m,i),
        %figure(i);
        surf(xmesh,tmesh,cww{i});xlabel('x轴');ylabel('t轴');title(strcat('节点',num2str(i),',w_',num2str(kss)));shading interp;
        figure(kss);hold on;
        %         facecolor={'interp','flat','texturemap','none','flat','interp'};
        %         edgecolor=linspace(0.5,1,N);
        surf(xmesh,tmesh,cww{i},cww{i}*10,'facecolor','interp');%,'facealpha',edgecolor(i));
        xlabel('x轴');ylabel('t轴');title(strcat('所有节点：w_',num2str(kss)));%shading interp;
        hold off;
        grid on;
    end
    
    %状态作图-分位置作图
    m=ceil(sqrt(nx));
    figure(lw*2+kss);
    names=cellmat(N,1);for i=1:N;names{i}=strcat('节点',num2str(i));end
    for i=1:nx
        i1=ind(kss,i,:);
        li=www(:,i1(:)');
        %li=ww(:,i*N-N+1:i*N);
        subplot(m,m,i);plot(tt,li);
        if i==1;legend(names);end
        if i>=1;title(strcat('x=',num2str(i/nx*(xr-xl)+xl),':w_',num2str(kss)));end
        if i==nx;xlabel('时间t(s)');end
    end
    if 1
        %触发事件作图-分位置分节点作图
        m=ceil(sqrt(nx));
        for j=1:min(3,N);%最多画3个节点的
            figure(lw*8+kss+j);
            names=cellmat(N,1);for i=1:N;names{i}=strcat('节点',num2str(i));end
            for i=1:nx
                i1=ind(kss,i,j);
                li=ckxt(:,i1(:)');
                subplot(m,m,i);plot(tt,li>0);axis([-Inf,Inf,-0.2,1.2]);
                if i==1;legend(names{j});end
                if i>=1;title(strcat('节点',num2str(j),',x=',num2str(i/nx*(xr-xl)+xl),':触发时间(K(x,t)>0)'));end
                if i==nx;xlabel('时间t(s)');end
            end
        end
        %K(x,t)作图-分位置作图
        m=ceil(sqrt(nx));
        figure(lw*6+kss+1);
        names=cellmat(N,1);for i=1:N;names{i}=strcat('节点',num2str(i));end
        for i=1:nx
            i1=ind(kss,i,:);
            li=ckxt(:,i1(:)');
            %li=ww(:,i*N-N+1:i*N);
            ep=1e-3;
            subplot(m,m,i);plot(tt,li>0);%axis([-Inf,Inf,-0.2,1.2]);
            if i==1;legend(names);end
            if i>=1;title(strcat('x=',num2str(i/nx*(xr-xl)+xl),':(K(x,t)>0)'));end
            if i==nx;xlabel('时间t(s)');end
        end
        
        %控制作图
        m=ceil(sqrt(nx));
        figure(lw*3+kss+1);
        names=cellmat(N,1);for i=1:N;names{i}=strcat('节点',num2str(i));end
        for i=1:nx
            i1=ind(kss,i,:);
            li=cu(:,i1(:)');
            %li=ww(:,i*N-N+1:i*N);
            subplot(m,m,i);plot(tt,li);
            if i==1;legend(names);end
            if i>=1;title(strcat('x=',num2str(i/nx*(xr-xl)+xl),':u(x,t)_',num2str(kss)));end
            if i==nx;xlabel('时间t(s)');end
        end
        
        %控制-分节点作图
        m=ceil(sqrt(N));[xmesh,tmesh]=meshgrid(xl+dx:dx:xr,tt);
        for i=1:N;
            i1=ind(kss,:,i);
            if 1;
                figure(lw*4+kss+1);
                subplot(m,m,i),
                %figure(i);
                surf(xmesh,tmesh,cu(:,i1));xlabel('x轴');ylabel('t轴');title(strcat('节点',num2str(i),',u(x,t)_',num2str(kss)));shading interp;
            end
            figure(lw*4+kss+lw+1);hold on;
            surf(xmesh,tmesh,cu(:,i1),cu(:,i1)*10,'facecolor','interp');%,'facealpha',edgecolor(i));
            xlabel('x轴');ylabel('t轴');title(strcat('所有节点：u(x,t)_',num2str(kss)));%shading interp;
            hold off;
            grid on;
        end
    end
end
if 1
    wf=www(end,:);
    odefun(0,wf);
end

%%
%函数
function [dw,Kxtall,uall]=fun(t,win,adj,B,delta,mu,dx,wzero,lw)
global D C;
global ffun gfun;
global wold;%wold是最近一次触发时间点的状态
global beta;%触发器参数
t,
N=size(adj,1);%节点个数
nx=numel(win)/N/lw;%x分段数目
w=reshape(win,lw,nx,N);
%x边界条件
wr=zeros(lw,N);wl=zeros(lw,N);%各个节点在x轴两端的值
%状态方程
dw=w;%导数
Kxtall=0*w;
uall=0*w;
for i=1:N;
    %控制
    ep=1e-2;
    phi=zeros(1,nx);
    for i1=1:N;
        erw=w(:,:,i)-w(:,:,i1);
        phi=phi+sqrt(sum((erw.*erw)*adj(i,i1),1)/lw);
    end
    %phi=sqrt(phi);
    phi_=ones(lw,1)*phi;
    %phi_1=sum(phi,1);
    zi=w(:,:,i)-wzero(:,:,i);
    avg_z2=0;
    for j=1:lw;
        avg_z2=avg_z2+(zi(j,:)*zi(j,:)')/nx/lw;
    end
    q=0.5;
    ui=-delta*zi./max(ep,abs(zi)).*phi_-mu*(avg_z2).^q*phi_.*zi./max(ep^1,abs(zi).^3);
    
    %事件触发
    phi=zeros(1,nx);
    for i1=1:N;
        erw=w(:,:,i)-wold(:,:,i1);
        phi=phi+sqrt(sum((erw.*erw)*adj(i,i1),1)/lw);%
    end
    er=wold(:,:,i)-w(:,:,i);%第i个节点的状态改变量
    Ni=sum(adj(i,:));
    Kxt=sqrt(sum(er.*er,1)/lw)-beta(i)*phi/Ni;
    ep=1e-3;Kxt=Kxt-ep;%即Kxt小于一定程度情况下可以不触发
    Kxt_=ones(lw,1)*Kxt;
    if 1
        %考虑事件触发
        uikxt=ui.*(Kxt_>=0);
    else%不考虑事件触发
        uikxt=ui;
    end
    %更新wold
    wold(:,:,i)=(Kxt_>0).*w(:,:,i)+(Kxt_<=0).*wold(:,:,i);%Kxt<=0的点不更新节点状态，>0的点更新
    %状态方程
    dw(:,:,i)=(C+D/dx^2)*[w(:,2:nx,i),wr(:,i)]+...%-(C+L)*[wzero(:,2:nx,i),wr(:,i)]+...
        +ffun([w(:,2:nx,i),wr(:,i)])+...
        -2*D/dx^2*[w(:,1:nx,i)]+...
        D/dx^2*[wl(:,i),w(:,1:nx-1,i)]+...
        uikxt;
    for j=1:N;
        dw(:,:,i)=dw(:,:,i)+B(i,j)*gfun([w(:,2:nx,j),wr(:,j)]);
    end
    Kxtall(:,:,i)=Kxt_;
    uall(:,:,i)=ui;
end
dw=dw(:);