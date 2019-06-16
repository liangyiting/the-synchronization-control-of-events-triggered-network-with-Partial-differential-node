function []=erreerer();
close all
lw=2;%w��ά��
global D C;%w''+cw=0  w'=w' w''=-cw�� D�����ϴ�ʱ����ʵ��ͬ��
D=[0.1,0;0,0.2];C=[0.1,1;-0.5,0.2];%ģ�����Գ���%Rk=max(abs(xl),abs(xr));%
global ffun gfun;
L=2;sigma=2;
ffun=@(w)L/2*(w+sin(1*w));    %����f����
gfun=@(w)sigma/2*(w+sin(5*w)+0*w./sqrt(w.^2+1));%����g����
%ʱ��
tf=0.2;%�յ�ʱ��

N=6;%�ڵ����
global beta;%�¼������������� betaԽСԽ�ȶ�
beta=0.1*ones(N,1);
%�����ڽӾ��� ������ɵ�
%�ڽ�Խ����Խ����ͬ��.������������Ӧ�ö��ڽ���Ҫ�󣬱�������ĳ��������С��������֮��ġ�
if 1
    %��������ڽӾ���
    adj=rand(N,N);adj=adj+adj';adj=adj/2;adj=[adj>0.5];
else %ȫ�����ڽӾ���
    adj=ones(N,N);
end
for i=1:N;adj(i,i)=0;end
%��������Ȩֵ,%����ȨֵԽ����Խ����ͬ��
if 0
    %�����������Ȩֵ ���Ȩֵ�п��ܳ����޷�ͬ��������
    B=10*adj.*(1+rand(N,N));%�໥Ӱ���ϵ��
else
    %�̶�����Ȩֵ
    B=5*adj;
end
for i=1:N;
    B(i,i)=-sum(B(i,:));
end;

%����������
delta=15;%����������
mu=5;%����������
%x����
xl=0;%��߽�
xr=1;%�ұ߽�
nx=20;%nx��x�ֶ���-- nx����̫�󣬷���������ֵ���ȶ���ƫ΢�ַ�������У�dt/dx����С��ĳһ��ֵ�Ż��ȶ�)
dx=(xr-xl)/nx;%����
maxdt=dx*0.1;%ƫ΢�ַ�������У�dt/dx����С��ĳһ��ֵ�Ż��ȶ�)

%��ʼ״̬w ƽ��״̬wzero  ������Ϊ�˼򻯣���wzeroȫ��Ϊ0�ˡ����������f������g������ʱ�����Ҫ��֤�����ƽ��㣬���f�����ó���L*(w+sin(w)),g�����ó���sigma*(w+sin(w))��
%w=ones(N,nx);%w(i*+j) �����ʼ��
w=zeros(lw,nx,N);
wzero=0*w;
for j=1:lw;
    %w(j,:,:)=0.5*rand(N,1)*[sin(pi*((1:nx)/nx))];%��֤����Ϊ0
    w(j,:,:)=0.15*[sin(pi*((1:nx)/nx))]'*[(1:N)-N/2];%��֤����Ϊ0
    wzero(j,:,:)=0*[sin(pi*((1:nx)/nx))]'*rand(1,N);%zeros(N,nx);%��㣬��������f��������ƫ�ģ�������Ϊ0�����ƽ���ֱ��Ϊ0
end

global wold;%wold�����һ�δ���ʱ����״̬
wold=w;
%%
%���ƫ΢�ַ�����
%�������ģ��
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
%��ͼ����
ind=reshape(1:lw*nx*N,lw,nx,N);
for kss=1:lw;%ֻȡ��һ��ά�ȵ�
    %ind1=ind(kss,:,:);ind1=ind1(:)';ww=www(:,ind1);%
    cww=cellmat(N,1);
    m=ceil(sqrt(N));
    [xmesh,tmesh]=meshgrid(xl+dx:dx:xr,tt);
    for i=1:N;
        i1=ind(kss,:,i);
        %cww{i}=ww(:,i:N:end);
        cww{i}=www(:,i1);
        %figure(1);subplot(li1,li2,i*2-1);pcolor(cww{i});xlabel('x��');ylabel('t��');title(strcat('�ڵ�',num2str(i)));shading interp;
        figure(kss+lw);
        subplot(m,m,i),
        %figure(i);
        surf(xmesh,tmesh,cww{i});xlabel('x��');ylabel('t��');title(strcat('�ڵ�',num2str(i),',w_',num2str(kss)));shading interp;
        figure(kss);hold on;
        %         facecolor={'interp','flat','texturemap','none','flat','interp'};
        %         edgecolor=linspace(0.5,1,N);
        surf(xmesh,tmesh,cww{i},cww{i}*10,'facecolor','interp');%,'facealpha',edgecolor(i));
        xlabel('x��');ylabel('t��');title(strcat('���нڵ㣺w_',num2str(kss)));%shading interp;
        hold off;
        grid on;
    end
    
    %״̬��ͼ-��λ����ͼ
    m=ceil(sqrt(nx));
    figure(lw*2+kss);
    names=cellmat(N,1);for i=1:N;names{i}=strcat('�ڵ�',num2str(i));end
    for i=1:nx
        i1=ind(kss,i,:);
        li=www(:,i1(:)');
        %li=ww(:,i*N-N+1:i*N);
        subplot(m,m,i);plot(tt,li);
        if i==1;legend(names);end
        if i>=1;title(strcat('x=',num2str(i/nx*(xr-xl)+xl),':w_',num2str(kss)));end
        if i==nx;xlabel('ʱ��t(s)');end
    end
    if 1
        %�����¼���ͼ-��λ�÷ֽڵ���ͼ
        m=ceil(sqrt(nx));
        for j=1:min(3,N);%��໭3���ڵ��
            figure(lw*8+kss+j);
            names=cellmat(N,1);for i=1:N;names{i}=strcat('�ڵ�',num2str(i));end
            for i=1:nx
                i1=ind(kss,i,j);
                li=ckxt(:,i1(:)');
                subplot(m,m,i);plot(tt,li>0);axis([-Inf,Inf,-0.2,1.2]);
                if i==1;legend(names{j});end
                if i>=1;title(strcat('�ڵ�',num2str(j),',x=',num2str(i/nx*(xr-xl)+xl),':����ʱ��(K(x,t)>0)'));end
                if i==nx;xlabel('ʱ��t(s)');end
            end
        end
        %K(x,t)��ͼ-��λ����ͼ
        m=ceil(sqrt(nx));
        figure(lw*6+kss+1);
        names=cellmat(N,1);for i=1:N;names{i}=strcat('�ڵ�',num2str(i));end
        for i=1:nx
            i1=ind(kss,i,:);
            li=ckxt(:,i1(:)');
            %li=ww(:,i*N-N+1:i*N);
            ep=1e-3;
            subplot(m,m,i);plot(tt,li>0);%axis([-Inf,Inf,-0.2,1.2]);
            if i==1;legend(names);end
            if i>=1;title(strcat('x=',num2str(i/nx*(xr-xl)+xl),':(K(x,t)>0)'));end
            if i==nx;xlabel('ʱ��t(s)');end
        end
        
        %������ͼ
        m=ceil(sqrt(nx));
        figure(lw*3+kss+1);
        names=cellmat(N,1);for i=1:N;names{i}=strcat('�ڵ�',num2str(i));end
        for i=1:nx
            i1=ind(kss,i,:);
            li=cu(:,i1(:)');
            %li=ww(:,i*N-N+1:i*N);
            subplot(m,m,i);plot(tt,li);
            if i==1;legend(names);end
            if i>=1;title(strcat('x=',num2str(i/nx*(xr-xl)+xl),':u(x,t)_',num2str(kss)));end
            if i==nx;xlabel('ʱ��t(s)');end
        end
        
        %����-�ֽڵ���ͼ
        m=ceil(sqrt(N));[xmesh,tmesh]=meshgrid(xl+dx:dx:xr,tt);
        for i=1:N;
            i1=ind(kss,:,i);
            if 1;
                figure(lw*4+kss+1);
                subplot(m,m,i),
                %figure(i);
                surf(xmesh,tmesh,cu(:,i1));xlabel('x��');ylabel('t��');title(strcat('�ڵ�',num2str(i),',u(x,t)_',num2str(kss)));shading interp;
            end
            figure(lw*4+kss+lw+1);hold on;
            surf(xmesh,tmesh,cu(:,i1),cu(:,i1)*10,'facecolor','interp');%,'facealpha',edgecolor(i));
            xlabel('x��');ylabel('t��');title(strcat('���нڵ㣺u(x,t)_',num2str(kss)));%shading interp;
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
%����
function [dw,Kxtall,uall]=fun(t,win,adj,B,delta,mu,dx,wzero,lw)
global D C;
global ffun gfun;
global wold;%wold�����һ�δ���ʱ����״̬
global beta;%����������
t,
N=size(adj,1);%�ڵ����
nx=numel(win)/N/lw;%x�ֶ���Ŀ
w=reshape(win,lw,nx,N);
%x�߽�����
wr=zeros(lw,N);wl=zeros(lw,N);%�����ڵ���x�����˵�ֵ
%״̬����
dw=w;%����
Kxtall=0*w;
uall=0*w;
for i=1:N;
    %����
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
    
    %�¼�����
    phi=zeros(1,nx);
    for i1=1:N;
        erw=w(:,:,i)-wold(:,:,i1);
        phi=phi+sqrt(sum((erw.*erw)*adj(i,i1),1)/lw);%
    end
    er=wold(:,:,i)-w(:,:,i);%��i���ڵ��״̬�ı���
    Ni=sum(adj(i,:));
    Kxt=sqrt(sum(er.*er,1)/lw)-beta(i)*phi/Ni;
    ep=1e-3;Kxt=Kxt-ep;%��KxtС��һ���̶�����¿��Բ�����
    Kxt_=ones(lw,1)*Kxt;
    if 1
        %�����¼�����
        uikxt=ui.*(Kxt_>=0);
    else%�������¼�����
        uikxt=ui;
    end
    %����wold
    wold(:,:,i)=(Kxt_>0).*w(:,:,i)+(Kxt_<=0).*wold(:,:,i);%Kxt<=0�ĵ㲻���½ڵ�״̬��>0�ĵ����
    %״̬����
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