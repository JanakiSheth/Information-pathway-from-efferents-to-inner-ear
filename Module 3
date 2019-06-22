clear variables;
%close all;

%% Parameters
%constants
F = 96490;				% C.mol-1; Faraday constant
kb = 1.381*10^-23;		%J.K-1; Boltzmann constant
T = 295;				%K; temperature
Ta = 1.5;
z = 2;					%valence of Ca2+
qe = 1.6*10^-19;		%C; charge of an electron
dca = 800.*10^-12;      %m2.s-1; diffusion coefficient of free Ca2+

pca = 1.*10^-18;		%m^3.s; Ca2+ permeability of the channel
vm0 = -55.*10^-3;		%V; resting membrane potential of hair cell
ccaext = 0.25*10^-3;     %M; extracellular Ca2+ concentration
deltaE = 65*10^(-21);

gamma = 0.14;			%geometric factor
n = 45;					%number of channels
d = 7.*10^-9;           %m; swing of transduction channel's gate
xihb = 200*10^-9;       %N.s.m-1; drag coefficient of hair bundle
xisf = 00*10^-9;        %N.s.m-1; drag coefficient of stimulus fiber
xia = 5*10^-6;
ksp = 300.*10^-6;		%N.m-1; combined stiffness of stereociliary pivots
xsp = 200*10^-9;        %offset deflection of stereociliary pivots
kes = 140.*10^-6;   	%N.m-1; stiffness of adaptation extent spring
xes = 40.*10^-9;         %m; offset of adaptation extent spring

p0 = 0.5;				%channel opening probability at rest
xc = 14*10^-9;          %m; extension of gating spring with channel closed
kgs0 = 2650*10^-6;   %N.m-1; stiffness in the absence of calcium
kgs1 = 2100*10^-6;       %N.m-1; slope of stiffness with slow relaxation
kre = 0*10^-6;     %N.m-1; slope of stiffness with fast relaxation
kgson = 5*10^6;        %s-1.m-1; rate of calcium binding to variable gating spring
kgsoff = 3.3;            %s-1.m-1; rate of calcium unbinding from the variable gating spring 

rm = 20.*10^-9;			%m; distance to adaptation motor's Ca2+-binding site
kmon = 10*10^6;          %s-1.M-1; rate constant for Ca2+ binding to adaptation motor
kmoff = 5*10^3;		%s-1; rate constant for Ca2+ unbinding from adaptation motor
smin = 0;               %m.s-1.N-1; minimum rate constant for motor slippage
smax = 400*10^3;       %m.s-1.N-1; maximum rate constant for motor slippage  
cmin = 0;               %m.s-1; minimum rate constant for motor climbing
cmax = 0.02*10^-6;      %m.s-1; maximum rate constant for motor climbing

rr = 10*10^-9;			%m; distance to relaxation element's Ca2+-binding site
kron = 2000.*10^6;		%s-1.M-1; rate constant for Ca2+ binding to relaxation element
kroff = 75*10^3;		%s-1; rate constant for Ca2+ unbinding from relaxation element

duration = 7;           %s; duration of the simulation
h = 0.0001;
ksf = 100*10^-6;                %N.m-1; stiffness of stimulus fiber

%projections and derived terms; all variables x, xa, etc are along x-axis	

A = exp(deltaE/(kb*T));	
B = pca*z*qe*ccaext/(kb*T*2*pi*dca*rm); 				%terms in calcium influx
delV = kb*T/(z*qe);
ca1 = -B*vm0/(1-exp(vm0/delV)); 						%constant for calcium dynamics
ca1r = (rm/rr)*ca1;                                     %constant for calcium dynamics at relaxation element
    
%offset parameters

ccam0 = ca1*p0;									%calcium concentration at the motor at rest
pm0 = 1/(1 + (kmoff/kmon)/ccam0) ;				%binding probability of calcium to the motor at rest
pgs0 = 1/(1 + (kgsoff/kgson)/ccam0) ;	
A1 = exp((kgs0 - kgs1*pgs0)*d*(xc - d/2)/(kb*T));
c0 = (1 - pm0)*(cmax - cmin) + cmin;
ccar0 = ca1r*p0;								%calcium concentration at the relaxation element at rest

%% arrays of variables 

t = zeros(duration/h,1);
X = zeros(duration/h,1); 
Xa = zeros(duration/h,1); 
prob = zeros(duration/h,1);
ccam = zeros(duration/h,1);
pm = zeros(duration/h,1);
ccar = zeros(duration/h,1);
pr = zeros(duration/h,1);
s = zeros(duration/h,1);
pgs = zeros(duration/h,1);
kgs = zeros(duration/h,1);
delta = zeros(duration/h,1);
offset = zeros(duration/h,1);

X(1) = 10^-11; Xa(1) = 0;

noise_channelclatter_amp = 0*sqrt(2*kb*T/(xihb + xisf));  % noise for the displacement function
noise_activemotors_amp = 0*sqrt(2*kb*T/xia)*gamma;    %noise for active motors

%% code follows Yuttana Bozovic 2011 paper
for i = 0/h+1:duration/h-1

noise_channelclatter = noise_channelclatter_amp*(randn);
noise_activemotors = noise_activemotors_amp*(randn);

% triangular offset of stimulus fiber
%{
if i > duration/(2*h)
   offset(i) = 1000*10^(-9)*(i - duration/(2*h))*(2*h)/duration - 1000*10^(-9); 
else
   offset(i) = -1000*10^(-9)*i*(2*h)/duration + 00*10^-9; 
end
%}
% offset of stimulus fiber
offset(i) = -500*10^-9;

%long step leading to recovery or showing hydrodynamic forces or step
%offset
%{
if (0 < i) && (i < 10000) 
    offset(i) = 000*10^-9;
elseif (20000 < i) && (i < 20500)   
    offset(i) = 1000*10^-9;
elseif (50000 < i) && (i < 50020)   
    offset(i) = 000*10^-9;    
else
    offset(i) = 0*10^-9;
end
%}
%time variant parameters  
delta(i+1) = kb*T/(kgs(i)*d);                     %term in open probability
prob(i+1) = 1/(1 + A1*exp(-(gamma*X(i) - Xa(i) + xc - d/2)/delta(i)));
ccam(i+1) = ca1*prob(i);                          %calcium concentration at the motor
pm(i+1) = 1/(1 + (kmoff/kmon)/ccam(i));           %binding probability of calcium to the motor 
s(i+1) = pm(i)*(smax - smin) + smin;              %rate of myosin slippage
ccar(i+1) = ca1r*prob(i);							%calcium concentration the relaxation element
pr(i+1) = 1/(1 + (kroff/kron)/ccar(i)) ;			%binding probability of calcium to the relaxation element
kgs(i+1) = kgs0 - kgs1*pgs(i) - kre*pr(i);        %gating spring stiffness, calcium dependence

X(i+1) = X(i) + h*((1/(xihb+xisf))*(ksf*(offset(i) - X(i)) - n*gamma*kgs(i)*(gamma*X(i) - Xa(i) +xc - prob(i)*d) - ksp*(X(i)-xsp))) + sqrt(h)*noise_channelclatter;
Xa(i+1) = Xa(i)+ h*(-c0 + s(i)*(kgs(i)*(gamma*X(i) - Xa(i) + xc - prob(i)*d) - kes*(Xa(i) - xes))) + sqrt(h)*noise_activemotors;
pgs(i+1) = pgs(i) + h*(kgson*ccam(i)*(1 - pgs(i)) - kgsoff*pgs(i));
t(i+1) = t(i)+h;

end

%% figures
f = figure();
ax = axes(f);
plot (ax,t(1:end),X(1:end));
hold on;
plot (ax,t(1:end),offset(1:end)/10);
title('Displacement vector');
xlabel('Time (s)');
ylabel('Displacement (m)');
%savefig(f, [save_directory, filesep, 'Displacement vector of spontaneous oscillations.fig']);
%print(f, '-dpng', '-r300',[save_directory, filesep, 'Displacement vector of spontaneous oscillations.png']);

f = figure();
ax = axes(f);
plot (ax,t(1:end),Xa(1:end));
title('Motor activity');
xlabel('Time (s)');
ylabel('Displacement of motors (m)');
%savefig(f, [save_directory, filesep, 'Motor activity of spontaneous oscillations.fig']);
%print(f, '-dpng', '-r300',[save_directory, filesep, 'Motor activity of spontaneous oscillations.png']);
