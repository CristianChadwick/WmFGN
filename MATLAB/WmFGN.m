

close all
clear all
clc





Data=importdata('historicos_biobio_CEN.csv');


Q=Data.data(1:end-1,4:end);

Q_3D=permute(reshape(Q,12,[],9),[2,1,3]);

mu_x=mean(Q_3D,1);
std_x=std(Q_3D,[],1);
std_y=sqrt(log(1+(std_x./mu_x).^2));
mu_y=log(mu_x)-std_y.^2./2;

std_y_3D=reshape(ones(size(Q_3D,1),1)*reshape(std_y,1,[]),size(Q_3D));
mu_y_3D=reshape(ones(size(Q_3D,1),1)*reshape(mu_y,1,[]),size(Q_3D));


Q_3D_N=(log(Q_3D)-mu_y_3D)./std_y_3D;

% Number of rundom numbers generated
n=10000;
Zc_3D=zeros(n,size(Q_3D,2),size(Q_3D,3));
XX=zeros(n+1,size(Q_3D,2),size(Q_3D,3));

for i=1:size(Q_3D,3)
    Q_aux=Q_3D_N(:,:,i);
    Q_aux_p=[Q_aux(1:end-1,end/2+1:end),Q_aux(1+1:end,1:end/2)];
    Q=chol(corr(Q_aux));
    Q_p=chol(corr(Q_aux_p));
    X=normrnd(0,1,n+1,size(Q_3D,2));
    XX(:,:,i)=X;
    X_p=[X(1:end-1,end/2+1:end),X(1+1:end,1:end/2)];
    Z=X*Q;
    Z_p=X_p*Q_p;
    Zc=[Z_p(:,end/2+1:end),Z(1+1:end,end/2+1:end)];
    Zc_3D(:,:,i)=Zc;
end


Zc_3D_M=zeros(size(Zc_3D));

for j=1:size(Q_3D,2)
    Q_aux=squeeze(Q_3D_N(:,j,:));
    Q=chol(corr(Q_aux));
    Q_aux_s=squeeze(Zc_3D(:,j,:));
    Zc_3D_c=Q_aux_s*Q;

    Zc_3D_M(:,j,:)=permute(Zc_3D_c,[1,3,2]);
end



XX_C=zeros(size(XX));

for j=1:size(Q_3D,2)
    Q_aux=squeeze(Q_3D_N(:,j,:));
    Q=chol(corr(Q_aux));
    Q_aux_s=squeeze(XX(:,j,:));
    XX_aux=Q_aux_s*Q;

    XX_C(:,j,:)=permute(XX_aux,[1,3,2]);
end

Zc_3DI=zeros(n,size(Q_3D,2),size(Q_3D,3));

for i=1:size(Q_3D,3)
    Q_aux=Q_3D_N(:,:,i);
    Q_aux_p=[Q_aux(1:end-1,end/2+1:end),Q_aux(1+1:end,1:end/2)];
    Q=chol(corr(Q_aux));
    Q_p=chol(corr(Q_aux_p)); 
    X=XX_C(:,:,i);
    X_p=[X(1:end-1,end/2+1:end),X(1+1:end,1:end/2)];
    Z=X*Q;
    Z_p=X_p*Q_p;
    Zc=[Z_p(:,end/2+1:end),Z(1+1:end,end/2+1:end)];
    Zc_3DI(:,:,i)=Zc;
end

Zc_3D_Ave=(Zc_3D_M+Zc_3DI)./2;

Zc_3D_SH=zeros(n,size(Q_3D,2),size(Q_3D,3));
ID_X=randi(size(Q_3D_N,1),n+1,size(Q_3D,2));

for i=1:size(Q_3D,3)
    Q_aux=Q_3D_N(:,:,i);
    Q_aux_p=[Q_aux(1:end-1,end/2+1:end),Q_aux(1+1:end,1:end/2)];
    Q=chol(corr(Q_aux));
    Q_p=chol(corr(Q_aux_p));
    X=zeros(n+1,size(Q_3D,2));
    for j=1:size(X,2)
        for k=1:size(X,1)
            X(k,j)=Q_3D_N(ID_X(k,j),j,i);
        end
    end
    X_p=[X(1:end-1,end/2+1:end),X(1+1:end,1:end/2)];
    Z=X*Q;
    Z_p=X_p*Q_p;
    Zc=[Z_p(:,end/2+1:end),Z(1+1:end,end/2+1:end)];
    Zc_3D_SH(:,:,i)=Zc;
end



% WmFGN:


alpha=0:0.01:1;
Err_Temp_Med=zeros(size(alpha));
Err_Spa_Med=zeros(size(alpha));
Err_Temp_Max=zeros(size(alpha));
Err_Spa_Max=zeros(size(alpha));
for k=1:length(alpha)
    
    Zc_3D_Ave=((1-alpha(k)).*Zc_3D_M+alpha(k).*Zc_3DI);
    
    E_corr_Temp=zeros(size(Zc_3DI,2),size(Zc_3DI,2),size(Zc_3DI,3));
    for i=1:size(Q_3D,3)
        Q_aux=Q_3D_N(:,:,i);
        Q_aux_S5=Zc_3D_Ave(:,:,i);
        E_corr_Temp(:,:,i)=abs(corr(Q_aux)-corr(Q_aux_S5));
    end
    Err_Temp_Med(k)=mean(mean(mean(E_corr_Temp)));
    Err_Temp_Max(k)=max(max(max(E_corr_Temp)));
    
    E_corr_Spa=zeros(size(Zc_3DI,3),size(Zc_3DI,3),size(Zc_3DI,2));
    for j=1:size(Q_3D,2)
        Q_aux=squeeze(Q_3D_N(:,j,:));
        Q_aux_S5=squeeze(Zc_3D_Ave(:,j,:));
        E_corr_Spa(:,:,j)=abs(corr(Q_aux)-corr(Q_aux_S5));
    end
    Err_Spa_Med(k)=mean(mean(mean(E_corr_Spa)));
    Err_Spa_Max(k)=max(max(max(E_corr_Spa)));
end


E_corr_Temp_mFGN=zeros(size(Zc_3DI,2),size(Zc_3DI,2),size(Zc_3DI,3));
for i=1:size(Q_3D,3)
    Q_aux=Q_3D_N(:,:,i);
    Q_aux_S3=Zc_3D_SH(:,:,i);
    E_corr_Temp_mFGN(:,:,i)=abs(corr(Q_aux)-corr(Q_aux_S3));
end
Err_Temp_Med_mFGN=mean(mean(mean(E_corr_Temp_mFGN)));
Err_Temp_Max_mFGN=max(max(max(E_corr_Temp_mFGN)));

E_corr_Spa_mFGN=zeros(size(Zc_3DI,3),size(Zc_3DI,3),size(Zc_3DI,2));
for j=1:size(Q_3D,2)
    Q_aux=squeeze(Q_3D_N(:,j,:));
    Q_aux_S3=squeeze(Zc_3D_SH(:,j,:));
    E_corr_Spa_mFGN(:,:,j)=abs(corr(Q_aux)-corr(Q_aux_S3));
end
Err_Spa_Med_mFGN=mean(mean(mean(E_corr_Spa_mFGN)));
Err_Spa_Max_mFGN=max(max(max(E_corr_Spa_mFGN)));



lambda=alpha;
A1=zeros(size(lambda));
B1=zeros(size(lambda));
A2=zeros(size(lambda));
B2=zeros(size(lambda));
Obj1_mFGN=zeros(size(lambda));
Obj2_mFGN=zeros(size(lambda));
for k=1:length(lambda)
    Obj1=lambda(k).*Err_Temp_Med+(1-lambda(k)).*Err_Spa_Med;
    [A1(k),B1(k)]=min(Obj1);
    Obj2=lambda(k).*Err_Temp_Max+(1-lambda(k)).*Err_Spa_Max;
    [A2(k),B2(k)]=min(Obj2);
    Obj1_mFGN(k)=lambda(k).*Err_Temp_Med_mFGN+(1-lambda(k)).*Err_Spa_Med_mFGN;
    Obj2_mFGN(k)=lambda(k).*Err_Temp_Max_mFGN+(1-lambda(k)).*Err_Spa_Max_mFGN;
end





















% Q_3D es el histórico
% Adoptando un alpha de 0.562 (que iguala los errores):
ALPHA=0.562;
Zc_3D_Ave=((1-ALPHA).*Zc_3D_M+ALPHA.*Zc_3DI);


std_y_3D_2=reshape(ones(size(Zc_3D_Ave,1),1)*reshape(std_y,1,[]),size(Zc_3D_Ave));
mu_y_3D_2=reshape(ones(size(Zc_3D_Ave,1),1)*reshape(mu_y,1,[]),size(Zc_3D_Ave));

Q_3D_gen=exp(Zc_3D_Ave.*std_y_3D_2+mu_y_3D_2);



Hsize= [[0.075   0.79    0.25    0.18];...
        [0.39    0.79    0.25    0.18];...
        [0.705   0.79    0.25    0.18];...
        [0.075   0.555   0.25    0.18];...
        [0.39    0.555   0.25    0.18];...
        [0.705   0.555   0.25    0.18];...
        [0.075   0.32    0.25    0.18];...
        [0.39    0.32    0.25    0.18];...
        [0.705   0.32    0.25    0.18];...
        [0.075   0.095   0.25    0.18];...
        [0.39    0.095   0.25    0.18];...
        [0.705   0.095   0.25    0.18];...
        [0.755   0.02    0.16    0.02]];

Prob_obs=(1:size(Q_3D,1))/(size(Q_3D,1)+1);
Prob_gen=(1:size(Q_3D_gen,1))/(size(Q_3D_gen,1)+1);

Let=[{'a)'},{'b)'},{'c)'},{'d)'},{'e)'},{'f)'},{'g)'},{'h)'},{'i)'},{'j)'},{'k)'},{'l)'}];
for i=1:size(Q_3D_gen,3)
    AUX_gen= Q_3D_gen(:,:,i);
    AUX_obs= Q_3D(:,:,i);
   
    figure('position',[0 0 1040 1250])
    h1=zeros(size(AUX_obs,2));
    for j=1:size(AUX_obs,2)
        h1(j)=subplot(4,3,j);
        plot(Prob_obs,sort(AUX_obs(:,j),'descend'),'LineWidth',2)
        hold on
        plot(Prob_gen,sort(AUX_gen(:,j),'descend'),'--','LineWidth',2)
        set(gca,'FontSize',12)
        AUX_y=max(AUX_gen(:,j))*1.2;
        axis([0 1 0 AUX_y])
        text(0.1*(1-(0))+0,1.08*(AUX_y-(0))+0,Let{j},'FontSize',12)
        if j==7
            text(-0.23*(1-(0))+0,0.85*(AUX_y-(0))+0,'Steamflow [m^3/s]','rotation',90,'FontSize',12)
        end
        if j==11
            xlabel('Probability','FontSize',12)
        end
        if j==12
            legend1=legend('Obs.','WmFGN','Orientation','horizontal');
        end
    end
    
    for j=1:size(AUX_obs,2)
        set(h1(j), 'Position', Hsize(j,:))
    end
    set(legend1, 'Position', Hsize(13,:))
    
end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %






SPsz=[[0.090    0.58    0.33    0.35];...
    [0.565     0.58     0.33    0.35];...
    [0.090     0.095    0.33    0.35];...
    [0.565     0.095    0.33    0.35];];

figure('position',[0,0,800 600])
H1=subplot(2,2,1);
plot(lambda,alpha(B1),'LineWidth',2,'col',[0.25,0.78,.37])
ylabel('\alpha')
xlabel('\lambda')
axis([0 1 0 1])
text(0.1*(1-(0))+0,1.08*(1-(0))+0,'a)','FontSize',12)
set(gca,'FontSize',12)
legend(strcat('Mean Error'),'location','best')


H2=subplot(2,2,2);
plot(lambda,alpha(B2),'LineWidth',2,'col',[0.75,0.78,.37])
ylabel('\alpha')
xlabel('\lambda')
axis([0 1 0 1])
text(0.1*(1-(0))+0,1.08*(1-(0))+0,'b)','FontSize',12)
set(gca,'FontSize',12)
legend(strcat('Max Error'),'location','NorthWest')


H3=subplot(2,2,3);
plot(lambda,A1,'LineWidth',2,'col',[0.02,0.46,.75])
hold on
plot(lambda,Obj1_mFGN,'--','col',[0.86,0.35,.13],'LineWidth',2)
ylabel('Objective Function 1')
xlabel('\lambda')
axis([0 1 0 0.08])
text(0.1*(1-(0))+0,1.08*(0.08-(0))+0,'c)','FontSize',12)
set(gca,'FontSize',12)
legend('WmFGN Mean Error','mFGN Mean Error','location','best')


H4=subplot(2,2,4);
plot(lambda,A2,'LineWidth',2,'col',[0.02,0.46,.75])
hold on
plot(lambda,Obj2_mFGN,'--','col',[0.86,0.35,.13],'LineWidth',2)
ylabel('Objective Function 2')
xlabel('\lambda')
axis([0 1 0 0.35])
text(0.1*(1-(0))+0,1.08*(0.35-(0))+0,'d)','FontSize',12)
set(gca,'FontSize',12)
legend('WmFGN Max Error','mFGN Max Error','location','best')

set(H1,'position',SPsz(1,:))
set(H2,'position',SPsz(2,:))
set(H3,'position',SPsz(3,:))
set(H4,'position',SPsz(4,:))


SPsz=[[0.090    0.58    0.28    0.35];...
    [0.465     0.58     0.28    0.35];...
    [0.090     0.095    0.28    0.35];...
    [0.465     0.095    0.28    0.35];...
    [0.837      0.74     0.08    0.08];...
    [0.837      0.2      0.08    0.08]];

figure('position',[0,0,900 600])
H1=subplot(2,2,1);
plot(alpha,Err_Temp_Med,'--','LineWidth',2,'color',[0.02,0.46,.75])
hold on
plot(alpha,Err_Spa_Med,'LineWidth',2,'color',[0.02,0.46,.75])
plot(alpha,Err_Temp_Med_mFGN*ones(1,length(alpha)),'--','LineWidth',2,'color',[0.86,0.35,.13])
plot(alpha,Err_Spa_Med_mFGN*ones(1,length(alpha)),'-','LineWidth',2,'color',[0.86,0.35,.13])
xlabel('\alpha')
ylabel('Mean Error')
axis([0 1 0 0.08])
text(0.1*(1-(0))+0,1.08*(0.08-(0))+0,'a)','FontSize',12)
set(gca,'FontSize',12)
% legend('WmFGN Temporal Error','WmFGN Spatial Error','mFGN Temporal Error','mFGN Error Spatial','location','best')

H2=subplot(2,2,2);
plot(alpha,Err_Temp_Max,'--','LineWidth',2,'color',[0.02,0.46,.75])
hold on
plot(alpha,Err_Spa_Max,'LineWidth',2,'color',[0.02,0.46,.75])
plot(alpha,Err_Temp_Max_mFGN*ones(1,length(alpha)),'--','LineWidth',2,'color',[0.86,0.35,.13])
plot(alpha,Err_Spa_Max_mFGN*ones(1,length(alpha)),'-','LineWidth',2,'color',[0.86,0.35,.13])
xlabel('\alpha')
ylabel('Maximum Error')
axis([0 1 0 0.35])
text(0.1*(1-(0))+0,1.08*(0.35-(0))+0,'b)','FontSize',12)
legend1=legend('WmFGN Temporal Error','WmFGN Spatial Error','mFGN Temporal Error','mFGN Spatial Error');
set(gca,'FontSize',12)

H3=subplot(2,2,3);
plot(Err_Spa_Med,Err_Temp_Med,'-.','LineWidth',2,'color',[0.02,0.46,.75])
hold on
plot(Err_Spa_Med_mFGN,Err_Temp_Med_mFGN,'o','LineWidth',2,'color',[0.86,0.35,.13])
axis([0 0.08 0 0.08])
text(0.1*(0.08-(0))+0,1.08*(0.08-(0))+0,'c)','FontSize',12)
xlabel('Mean Spatial Error')
ylabel('Mean Temporal Error')
set(gca,'FontSize',12)
% legend('WmFGN','mFGN','location','best')

H4=subplot(2,2,4);
plot(Err_Spa_Max,Err_Temp_Max,'-.','LineWidth',2,'color',[0.02,0.46,.75])
hold on
plot(Err_Spa_Max_mFGN,Err_Temp_Max_mFGN,'o','LineWidth',2,'color',[0.86,0.35,.13])
axis([0 0.35 0 0.35])
text(0.1*(0.35-(0))+0,1.08*(0.35-(0))+0,'d)','FontSize',12)
xlabel('Max Spatial Error')
ylabel('Max Temporal Error')
legend2=legend('WmFGN','mFGN');
set(gca,'FontSize',12)

set(H1,'position',SPsz(1,:))
set(H2,'position',SPsz(2,:))
set(H3,'position',SPsz(3,:))
set(H4,'position',SPsz(4,:))
set(legend1,'Position',SPsz(5,:))
set(legend2,'Position',SPsz(6,:))





Hsize= [[0.1200    0.7184    0.1550    0.2457];...
        [0.4100    0.7184    0.1550    0.2457];...
        [0.7000    0.7184    0.1550    0.2457];...
        [0.03      0.3862    0.1550    0.2457];...
        [0.2200    0.3862    0.1550    0.2457];...
        [0.4100    0.3862    0.1550    0.2457];...
        [0.6000    0.3862    0.1550    0.2457];...
        [0.7900    0.3862    0.1550    0.2457];...
        [0.03      0.054     0.1550    0.2457];...
        [0.2200    0.054     0.1550    0.2457];...
        [0.4100    0.054     0.1550    0.2457];...
        [0.6000    0.054     0.1550    0.2457];...
        [0.7900    0.054     0.1550    0.2457]];


Sel_lam = [0,0.25,0.5,0.75,1];
Aux=zeros(length(Sel_lam),1);
for i=1:length(Sel_lam)
    Aux(i)=find(lambda==Sel_lam(i));
end

for kkk=1:2
Zc_3D_Ave_Mean=zeros([size(Zc_3D_M),length(Sel_lam)]);
for i=1:length(Sel_lam)
    if kkk==1
        Zc_3D_Ave_Mean(:,:,:,i)=((1-alpha(B1(Aux(i)))).*Zc_3D_M+alpha(B1(Aux(i))).*Zc_3DI);
    elseif kkk==2
        Zc_3D_Ave_Mean(:,:,:,i)=((1-alpha(B2(Aux(i)))).*Zc_3D_M+alpha(B2(Aux(i))).*Zc_3DI);
    end
end

for i=1:size(Q_3D,3)

    Q_aux=Q_3D_N(:,:,i);
    Q_aux_S=Zc_3D_SH(:,:,i);
    Q_aux_S1=Zc_3D_Ave_Mean(:,:,i,1);
    Q_aux_S2=Zc_3D_Ave_Mean(:,:,i,2);
    Q_aux_S3=Zc_3D_Ave_Mean(:,:,i,3);
    Q_aux_S4=Zc_3D_Ave_Mean(:,:,i,4);
    Q_aux_S5=Zc_3D_Ave_Mean(:,:,i,5);
    
    NN= size(Q_aux,2);
    x=1:NN;
    y=1:NN;
    [X,Y]=meshgrid(x,y);
    
    
    figure('position',[0 0 1570 850])
    h1=subplot(3,5,1.5);
%     contourf(corr(Q_aux))
    Z=corr(Q_aux);
    surf(X,Y,Z)
    view(2)
    title('obs')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
    colorbar
    caxis([-1 1])
    text(1,0,'a)','FontSize',12)
    
    h2=subplot(3,5,3);
%     contourf(corr(Q_aux_S))
    Z=corr(Q_aux_S);
    surf(X,Y,Z)
    view(2)
    title('mFGN')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
    colorbar
    caxis([-1 1])
    text(1,0,'b)','FontSize',12)
    
    h3=subplot(3,5,4.5);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S));
    surf(X,Y,Z)
    view(2)
    title('|obs- mFGN|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
    colorbar
    caxis([0 0.4])
    text(1,0,'c)','FontSize',12)
    
    
    h4=subplot(3,5,6);
%     contourf(corr(Q_aux_S1))
    Z=corr(Q_aux_S1);
    surf(X,Y,Z)
    view(2)
    title('WmFGN \lambda = 0.00')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([-1 1])
    text(1,0,'d)','FontSize',12)
    
    h5=subplot(3,5,7);
%     contourf(corr(Q_aux_S2))
    Z=corr(Q_aux_S2);
    surf(X,Y,Z)
    view(2)
    title('WmFGN \lambda = 0.25')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([-1 1])
    text(1,0,'e)','FontSize',12)
    
    h6=subplot(3,5,8);
%     contourf(corr(Q_aux_S3))
    Z=corr(Q_aux_S3);
    surf(X,Y,Z)
    view(2)
    title('WmFGN \lambda = 0.50')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([-1 1])
    text(1,0,'f)','FontSize',12)
    
    h7=subplot(3,5,9);
%     contourf(corr(Q_aux_S4))
    Z=corr(Q_aux_S4);
    surf(X,Y,Z)
    view(2)
    title('WmFGN \lambda = 0.75')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([-1 1])
    text(1,0,'g)','FontSize',12)
    
    h8=subplot(3,5,10);
%     contourf(corr(Q_aux_S5))
    Z=corr(Q_aux_S5);
    surf(X,Y,Z)
    view(2)
    title('WmFGN \lambda = 1.00')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
    colorbar
    caxis([-1 1])
    text(1,0,'h)','FontSize',12)
    
    
    h9=subplot(3,5,11);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S1)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S1));
    surf(X,Y,Z)
    view(2)
    title('|obs- WmFGN \lambda = 0.00|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([0 0.4])
    text(1,0,'i)','FontSize',12)
    
    h10=subplot(3,5,12);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S2)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S2));
    surf(X,Y,Z)
    view(2)
    title('|obs- WmFGN \lambda = 0.25|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([0 0.4])
    text(1,0,'j)','FontSize',12)
    
    h11=subplot(3,5,13);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S3)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S3));
    surf(X,Y,Z)
    view(2)
    title('|obs- WmFGN \lambda = 0.50|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([0 0.4])
    text(1,0,'k)','FontSize',12)
    
    h12=subplot(3,5,14);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S4)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S4));
    surf(X,Y,Z)
    view(2)
    title('|obs- WmFGN \lambda = 0.75|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([0 0.4])
    text(1,0,'l)','FontSize',12)
    
    h13=subplot(3,5,15);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S5)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S5));
    surf(X,Y,Z)
    view(2)
    title('|obs- WmFGN \lambda = 1.00|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
    colorbar
    caxis([0 0.4])
    text(1,0,'m)','FontSize',12)
    
    
    set(h1, 'Position', Hsize(1,:))
    set(h2, 'Position', Hsize(2,:))
    set(h3, 'Position', Hsize(3,:))
    set(h4, 'Position', Hsize(4,:))
    set(h5, 'Position', Hsize(5,:))
    set(h6, 'Position', Hsize(6,:))
    set(h7, 'Position', Hsize(7,:))
    set(h8, 'Position', Hsize(8,:))
    set(h9, 'Position', Hsize(9,:))
    set(h10, 'Position', Hsize(10,:))
    set(h11, 'Position', Hsize(11,:))
    set(h12, 'Position', Hsize(12,:))
    set(h13, 'Position', Hsize(13,:))
    
end





for j=1:size(Q_3D,2)
    
    
    Q_aux=squeeze(Q_3D_N(:,j,:));
    Q_aux_S=squeeze(Zc_3D_SH(:,j,:));
    Q_aux_S1=squeeze(Zc_3D_Ave_Mean(:,j,:,1));
    Q_aux_S2=squeeze(Zc_3D_Ave_Mean(:,j,:,2));
    Q_aux_S3=squeeze(Zc_3D_Ave_Mean(:,j,:,3));
    Q_aux_S4=squeeze(Zc_3D_Ave_Mean(:,j,:,4));
    Q_aux_S5=squeeze(Zc_3D_Ave_Mean(:,j,:,5));
    
    NN= size(Q_aux,2);
    x=1:NN;
    y=1:NN;
    [X,Y]=meshgrid(x,y);

    figure('position',[0 0 1570 850])
    h1=subplot(3,5,1.5);
%     contourf(corr(Q_aux))
    Z=corr(Q_aux);
    surf(X,Y,Z)
    view(2)
    title('obs')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
    colorbar
    caxis([-1 1])
    text(1,0,'a)','FontSize',12)
    
    h2=subplot(3,5,3);
%     contourf(corr(Q_aux_S))
    Z=corr(Q_aux_S);
    surf(X,Y,Z)
    view(2)
    title('mFGN')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
    colorbar
    caxis([-1 1])
    text(1,0,'b)','FontSize',12)
    
    h3=subplot(3,5,4.5);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S));
    surf(X,Y,Z)
    view(2)
    title('|obs- mFGN|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
    colorbar
    caxis([0 0.4])
    text(1,0,'c)','FontSize',12)
    
    
    h4=subplot(3,5,6);
%     contourf(corr(Q_aux_S1))
    Z=corr(Q_aux_S1);
    surf(X,Y,Z)
    view(2)
    title('WmFGN \lambda = 0.00')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([-1 1])
    text(1,0,'d)','FontSize',12)
    
    h5=subplot(3,5,7);
%     contourf(corr(Q_aux_S2))
    Z=corr(Q_aux_S2);
    surf(X,Y,Z)
    view(2)
    title('WmFGN \lambda = 0.25')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([-1 1])
    text(1,0,'e)','FontSize',12)
    
    h6=subplot(3,5,8);
%     contourf(corr(Q_aux_S3))
    Z=corr(Q_aux_S3);
    surf(X,Y,Z)
    view(2)
    title('WmFGN \lambda = 0.50')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([-1 1])
    text(1,0,'f)','FontSize',12)
    
    h7=subplot(3,5,9);
%     contourf(corr(Q_aux_S4))
    Z=corr(Q_aux_S4);
    surf(X,Y,Z)
    view(2)
    title('WmFGN \lambda = 0.75')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([-1 1])
    text(1,0,'g)','FontSize',12)
    
    h8=subplot(3,5,10);
%     contourf(corr(Q_aux_S5))
    Z=corr(Q_aux_S5);
    surf(X,Y,Z)
    view(2)
    title('WmFGN \lambda = 1.00')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
    colorbar
    caxis([-1 1])
    text(1,0,'h)','FontSize',12)
    
    
    h9=subplot(3,5,11);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S1)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S1));
    surf(X,Y,Z)
    view(2)
    title('|obs- WmFGN \lambda = 0.00|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([0 0.4])
    text(1,0,'i)','FontSize',12)
    
    h10=subplot(3,5,12);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S2)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S2));
    surf(X,Y,Z)
    view(2)
    title('|obs- WmFGN \lambda = 0.25|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([0 0.4])
    text(1,0,'j)','FontSize',12)
    
    h11=subplot(3,5,13);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S3)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S3));
    surf(X,Y,Z)
    view(2)
    title('|obs- WmFGN \lambda = 0.50|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([0 0.4])
    text(1,0,'k)','FontSize',12)
    
    h12=subplot(3,5,14);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S4)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S4));
    surf(X,Y,Z)
    view(2)
    title('|obs- WmFGN \lambda = 0.75|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
%     colorbar
    caxis([0 0.4])
    text(1,0,'l)','FontSize',12)
    
    h13=subplot(3,5,15);
%     contourf(abs(corr(Q_aux)-corr(Q_aux_S5)))
    Z=abs(corr(Q_aux)-corr(Q_aux_S5));
    surf(X,Y,Z)
    view(2)
    title('|obs- WmFGN \lambda = 1.00|')
    set(gca, 'YDir','reverse','FontSize',12)
    axis([1 NN 1 NN])
    colorbar
    caxis([0 0.4])
    text(1,0,'m)','FontSize',12)
    
    
    set(h1, 'Position', Hsize(1,:))
    set(h2, 'Position', Hsize(2,:))
    set(h3, 'Position', Hsize(3,:))
    set(h4, 'Position', Hsize(4,:))
    set(h5, 'Position', Hsize(5,:))
    set(h6, 'Position', Hsize(6,:))
    set(h7, 'Position', Hsize(7,:))
    set(h8, 'Position', Hsize(8,:))
    set(h9, 'Position', Hsize(9,:))
    set(h10, 'Position', Hsize(10,:))
    set(h11, 'Position', Hsize(11,:))
    set(h12, 'Position', Hsize(12,:))
    set(h13, 'Position', Hsize(13,:))
end


end





