clear;clc;close all;

%% Paramètres
Fe=2e7;
Ts = 1e-6;
Fse=Ts*Fe;
bk = [1,0,0,1,0];

%% Emetteur

p_t = 0.5*([-1*ones(Fse/2,1);ones(Fse/2,1)]); %Filtre

A_k=real(pskmod(bk, 2)); %Symboles

sl_t=A_k(1).*p_t;
for i=A_k(:,2:5)
    sl_t = [sl_t ; i.*p_t];
end
sl_t = sl_t+0.5;

%% Recepteur

rl_t= conv(sl_t,(-p_t));
Debut=20;
Fin=length(bk)*Fse;
rm= rl_t(Debut:Fse:Fin);

%rm=downsample(rl_t,Fse);

bk_chap=zeros(5,1);

plot(rm)
for j=1:5
    if rm(j)>1
        bk_chap(j)=0;
    else
        bk_chap(j)=1;
    end
end

        
