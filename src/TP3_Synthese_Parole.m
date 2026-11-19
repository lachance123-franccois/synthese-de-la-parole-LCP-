%% ANALYSE SPECTRALE DES SIGNAUX ET DES SYSTEMES
%% Analyse, codage et synthèse de la parole
%% I. Analyse
% tracon la reponse temporelle du signal du signal audio

%lecture du signal
son = audioread('sig.wav');
Fe = 11000;
N = length(son);
t = (0:(N-1))/Fe;

%affichage
figure;
plot(t, son);
title('audiogramme du signal');


%ecoute du signal
%sound(son);
%pause(2);


%suppression des zones de silence
n1 = round(0.3*Fe);
n2 = round(1.03*Fe);

x = son(n1:n2);
Nn = length(x);
tn = (0:(Nn-1))/Fe;

%affichage
figure;
plot(tn, x);
title('audiogramme signal silence supprimé');

%ecoute du signal
%sound(x);
pause(2);

% zomme sur chaque phoneme





%phoneme 1
%suppression des zones de silence
n1 = round(0.3*Fe);
n2 = round(0.5*Fe);

x1 = son(n1:n2);
Nn1 = length(x1);
tn = (0:(Nn1-1))/Fe;

%affichage
figure;
plot(tn, x1);
title('phoneme a');

%ecoute du signal
%sound(x1);
pause(2);

% zomme sur chaque phoneme

%phoneme 2
%suppression des zones de silence
n1 = round(0.5*Fe);
n2 = round(0.6*Fe);

x2 = son(n1:n2);
Nn2 = length(x2);
tn = (0:(Nn2-1))/Fe;

%affichage
figure;
plot(tn, x2);
title('phoneme ch');

%ecoute du signal
%sound(x2);
pause(2);

% zomme sur le phoneme é
close all
%phoneme 3
%suppression des zones de silence
n1 = round(0.61*Fe);
n2 = round(0.72*Fe);

x3 = son(n1:n2);
Nn3 = length(x3);
tn = (0:(Nn3-1))/Fe;

%affichage
figure;
plot(tn, x3);
title('phoneme é');

%ecoute du signal
%sound(x3);
pause(2);

% zomme sur le phoneme ch

%phoneme 4
%suppression des zones de silence
n1 = round(0.72*Fe);
n2 = round(0.86*Fe);

x4 = son(n1:n2);
Nn4 = length(x4);
tn = (0:(Nn4-1))/Fe;

%affichage
figure;
plot(tn, x4);
title('phoneme ch');

%ecoute du signal
%sound(x4);
pause(2);

% zomme sur le phoneme é

%phoneme 5
%suppression des zones de silence
n1 = round(0.86*Fe);
n2 = round(1.03*Fe);

x5 = son(n1:n2);
Nn5 = length(x5);
tn = (0:(Nn5-1))/Fe;

%affichage
figure;
plot(tn, x5);
title('phoneme é');

%ecoute du signal
%sound(x5);
pause(2);

%tracons le module de la trensformer de fourier 
N = length(x);
X = fftshift(fft(x, N));

f = ((-N/2):(N/2)-1)*(Fe/N);

%affichage
figure;
plot(f, abs(X));
title('module de la TF du signal X')

% en observant le spectre l'on se rend compte que l'on perd l'information
% temporelle car le signal est un signal non-stationnaire donc le module de
% sa transforme de fourier ne peut plus nous renseigner sur les frequences
% respective des different phoneme, d'ou l'interet de faire l'analyse temps
% frequence

% analyse temps frequence

%spectrogram a bande etroite
figure;
spectrogram(x,0.05*Fe,0.025*Fe,0.05*Fe,Fe,'yaxis');
title('spectrogram a bande etroite 0.05')


%spectrogram a large bande
figure;
spectrogram(x,0.01*Fe,0.005*Fe,0.01*Fe,Fe,'yaxis');
title('spectrogram a bande large 0.01');


% pour l'analyse des different spectrogramme nous nous rendons que
% lorsqu'on prend une large fenetre on a une bonne resolution fréquentielle
% et une mauvaise resolution temporelle, donc il y'a un comprime a faire.
% sur le spectrogramme a bande etroit nous reperons la frequence
% fondemental qui nous renseigne sur la periode pitch et le spectrograme a
% large bande nous observons la localisation temporelle des different
% phoneme que nous avions au depart. les phoneme non-voise n'ont pas de
% frequence fondamentale puisqu'elle sont considere comme un bruit.

%% 2 Codage LPC

% estimons des coeficien du modele PLC
           
M = round(0.03 * Fe);  
tranches = floor(length(x)/M);

matrice_tranche = zeros(tranches, M);

puissance = zeros(tranches, 1);
pitch = zeros(tranches, 1);

%pour le calcul des coeficient
P=12;
tmin = round((1/600)*Fe);
tmax = round((1/70)*Fe);


Coef = zeros(P, tranches);
ind_n1 = round(0.18*Fe);
ind_n2 = round(0.3*Fe);

ind_n3 = round(0.41*Fe);
ind_n4 = round(0.56*Fe);
for i = 1:tranches
    tranche = x(((i-1)*M+1) : i*M);
    matrice_tranche(i, :) = tranche;

    puissance(i) = mean(tranche.^2);

    % autocorrelation
    r = xcorr(tranche, M,'biased'); 

    intervale = ((1/600):(1/11000):(1/70));


    rp= r(M+1:end);

    [~, ind] = max(rp(tmin:tmax));


 

    %condition si la tranche se trouve dans un phoneme voise ou non voise
    if (((i-1)*M+1) > ind_n1)  && (i*M<ind_n2) ||  (((i-1)*M+1) > ind_n3)  && (i*M<ind_n4) 
        pitch(i) = 0;
    else
        pitch(i) = tmin + ind(1)-1;
    end
    
   

    % calcul des coeficien LPC d'un modele d'odre P = 12
 
    rp_lp = r(M+1:M+P);
    rh_lp= -r(M+2:M+P+1);
    R = toeplitz(rp_lp);
    Coef(:, i) = R\rh_lp;

end
%
figure; 
stem(pitch);
title('Pitch estimE par tranche');
ylabel('echantillons');
xlabel('Tranche');



% le debit pour transmetre le signal de depart est : 11000 * 8= 88000
% bits/s
% en faisant le codage avec le LPC
% on code les coefitien(12), la puissance, et le pitch: (12+1+1) * 8 = 112
% pour tout les tranches on a : 112 * 24 = 2688 bits/s
% nous constatons que avec le codage LPC on a besoin de moins de debit pour
% transmetre les donnees

%% 3 Synthese
% Synthese du signal 
%signal syntetiser
x_synth = zeros(tranches *M, 1);
zi = zeros(P, 1);
fenetre = hann(M);
for i =1:tranches
   
    if pitch(i) == 0
       source = randn(M, 1);
       source=source/std(source);
    else
        source=zeros(M,1);
        source(1:pitch(i):end) = 1;
    end
    
    % filtrage
    a = [1; Coef(:, i)];

    [x_tranche, zi] = filter(1, a, source, zi);

    % Ajustement de la puissance
    p_tranche = mean(x_tranche .^ 2);
    if p_tranche > 0 && puissance(i) > 0
        x_tranche = x_tranche * sqrt(puissance(i) / p_tranche);
    end
  
    x_synth((i-1)*M+1 : i*M) = x_tranche;
end
% %% reponses en frequence
N_fft = 512;
[H, f_h] = freqz(1, a, N_fft/2, Fe);
plot(f_h, 20*log10(abs(H)), 'r', 'LineWidth', 1.5);
%
%puissance du signal
figure;
plot(puissance);
xlabel('Tranche');
ylabel('Puissance');
title('Puissance par tranche');

% Ecoute du son synthetiser
sound(x_synth, Fe);
pause(3);

figure;
subplot(2,1,1);
plot(x);      
title('Signal original');
subplot(2,1,2); 
plot(x_synth); 
title('Signal synthétisé');

figure;
subplot(2,1,1);
spectrogram(x, 0.01*Fe, 0.005*Fe, 0.01*Fe, Fe, 'yaxis');
title('Spectrogramme original');
subplot(2,1,2);
spectrogram(x_synth, 0.01*Fe,0.005*Fe,0.01*Fe, Fe, 'yaxis');
title('Spectrogramme synthétisé');


M = round(0.03 * Fe);       
hop= M / 2;               
fenetre= hann(M);                

% Nombre de tranches avec recouvrement
tranches = floor((length(x) - M) / hop) + 1;

puissance = zeros(tranches, 1);
pitch= zeros(tranches, 1);
P = 12;
Coef= zeros(P, tranches);

tmin = round(Fe / 600);
tmax = round(Fe / 70);

ind_n1 = round(0.18 * Fe);
ind_n2 = round(0.30 * Fe);
ind_n3 = round(0.41 * Fe);
ind_n4 = round(0.56 * Fe);

for i = 1:tranches

    debut = (i-1)*hop + 1;
    fin   = debut + M - 1;
    tranche = x(debut:fin) .* fenetre;   
    puissance(i) = mean(tranche .^ 2);

    r = xcorr(tranche, M, 'biased');

    % Pitch
    rp_pitch = r(M+1:end);
    [~, ind] = max(rp_pitch(tmin:tmax));

    non_voise = (debut > ind_n1 && fin < ind_n2) || ...
                (debut > ind_n3 && fin < ind_n4);
    if non_voise
        pitch(i) = 0;
    else
        pitch(i) = tmin + ind(1) - 1;
    end

    % Coefficients LPC
    rp_lpc = r(M+1 : M+P);
    rh_lpc = -r(M+2 : M+P+1);
    R      = toeplitz(rp_lpc);
    Coef(:, i) = R \ rh_lpc;

end

x_synth2 = zeros(length(x) + M, 1);   
zi       = zeros(P, 1);               
for i = 1:tranches

    if pitch(i) == 0
        source = randn(M, 1);
        source = source / std(source);
    else
        source = zeros(M, 1);
        source(1:pitch(i):end) = 1;
    end

    a = [1; Coef(:, i)];
    [x_tranche, zi] = filter(1, a, source, zi);

    % Ajustement puissance
    p_tr = mean(x_tranche .^ 2);
    if p_tr > 0 && puissance(i) > 0
        x_tranche = x_tranche * sqrt(puissance(i) / p_tr);
    end

    % Application fenêtre de Hanning +  recouvrement
    debut = (i-1)*hop + 1;
    x_synth2(debut : debut+M-1) = x_synth2(debut : debut+M-1) + ...
                                  x_tranche .* fenetre;
end

% Tronquer au longueur originale
x_synth2= x_synth2(1:length(x));

% Normalisation 
x_synth2 = x_synth2 / max(abs(x_synth2)) * max(abs(x));

sound(x_synth2, Fe);

% comprarison entre les signal original et celui avec la methode de hanning
figure;
subplot(2,1,1);
plot(x);      
title('Signal original');
subplot(2,1,2); 
plot(x_synth2); 
title('Signal synthétisé');

figure;
subplot(2,1,1);
spectrogram(x, 0.01*Fe, 0.005*Fe, 0.01*Fe, Fe, 'yaxis');
title('Spectrogramme original');
subplot(2,1,2);
spectrogram(x_synth2, 0.01*Fe,0.005*Fe,0.01*Fe, Fe, 'yaxis');
title('Spectrogramme synthétisé');

