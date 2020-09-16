clear;clc;close all;

%% Paramètres
Fe=2e7;
Ts = 1e-6;
Nb=1000;
Fse=Ts*Fe;


eb_n0_db=0:1:10;
eb_n0=10.^(eb_n0_db/10);
TEB=zeros(size(eb_n0));
Pb=qfunc(sqrt(2*eb_n0));
%% Emetteur
for k=1:length(eb_n0_db)
    error_cnt=0;
    bit_cnt=0;
    while error_cnt<100
        
p_t = 0.5*([-1*ones(Fse/2,1);ones(Fse/2,1)]); %Filtre
%bk = [1,0,0,1,0];
bk = randi([0,1],[1,Nb]);
A_k=real(pskmod(bk, 2)); %Symboles

sl_t=A_k(1).*p_t;
for i=A_k(:,2:Nb)
    sl_t = [sl_t ; i.*p_t];
end
sl_t = sl_t+0.5;
nl_t=sqrt(1/2)*randn(20000,1);
yl_t=sl_t+nl_t;
%% Recepteur

rl_t= conv(yl_t,(-p_t));
Debut=20;
Fin=length(bk)*Fse;
rm= rl_t(Debut:Fse:Fin);

%rm=downsample(rl_t,Fse);

bk_chap=zeros(Nb,1);
erreur=0;
for j=1:Nb
    if rm(j)>1
        bk_chap(j)=0;
        if bk(j)==1
            erreur=erreur+1;
        end
    else
        bk_chap(j)=1;
        if bk(j)==0
            erreur=erreur+1;
        end
    end
end
error_cnt=error_cnt+erreur
bit_cnt=bit_cnt+Nb
    end
    TEB(k)=error_cnt/bit_cnt;
end
plot(TEB)
%plot(Pb)

