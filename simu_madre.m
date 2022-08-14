%%  PRIMERA SIMULACION. 
tic
addpath('/home/rosario/');
addpath('/home/rosario/PlantPollinatorModel');
addpath('/home/rosario/PlantPollinatorModel/src/NetworkGeneration');
addpath('/home/rosario/form_redes_otros');

% Valor de semilla aleatoria
    replica = 0;
    SEED = 0; 

path_Files = '/home/rosario/form_redes_otros/REDES_FINAL/';
 
file_Network = strcat(path_Files,'net054.mat'); 

load(file_Network); 

 
for num_net = 1 %100

BIOMASA_FINAL_S1 = 0 ;
while sum(BIOMASA_FINAL_S1) == 0 

NET_ORIGINAL = Networks{num_net};

ph_L = 12; 
ph_C = 6; 
disp('Creando Fenologias ...')
   
algSteps = numel(NET_ORIGINAL)*7*ph_L;
    
[V_PLANT, V_POLLI] = creadorFenofases_v06(NET_ORIGINAL, algSteps, ph_L, ph_C);
    
disp('Fenologias OK')
 
Fenofases.I = V_POLLI;
Fenofases.P = V_PLANT;

[nPolinizadores, nPlantas] = size(NET_ORIGINAL);
ESPECIES.TOTAL = nPolinizadores+nPlantas;
ESPECIES.PLANTS = nPlantas;
ESPECIES.POLLINATORS = nPolinizadores;

RANGE_POP = 0.25;

BIOMASS = bioMassInit(nPolinizadores,nPlantas,RANGE_POP); 
 
[~,p] = size(NET_ORIGINAL);
ceros = zeros(1,p);
NET_FALSA = [NET_ORIGINAL; ceros];

RANGE_VAR = 0.25;
basalValues = meanParameterValues;
PARAMS_BORRADOR = distrMBeta(NET_FALSA, basalValues, RANGE_VAR, replica);


PARAMS_NAT = parNat(PARAMS_BORRADOR);

%% Definimos parametros de la simulacion
    STEPS = 5000; % Cien a~nos aprox.
    SIM = 1; % Forma en que cambian las fenofases (Solo si "shift" =/= 0)
    
 %% Definimos parametros para creacion de archivos
    
    % Indica si crear archivos de salida
    OUT = true;
    OUTPUT_DIR = '/home/rosario/Desktop/resu_parcial/';
    
    ID_SIM_1 = ['_eiaa_', num2str(num_net)];  
    %% Corremos el experimento
    
    % Tratamiento: Control (solo sps nativas)
    
    SHIFT = 0; % No queremos cambios en las fenofases
    datos_simulacion_final.id_sim = ID_SIM_1;
    datos_simulacion_final.biomasas = BIOMASS;
    datos_simulacion_final.steps = STEPS;
    datos_simulacion_final.red = NET_ORIGINAL;
    datos_simulacion_final.fenos = Fenofases;
    datos_simulacion_final.sim = SIM;
    datos_simulacion_final.shift = SHIFT;
    datos_simulacion_final.output_dir = OUTPUT_DIR;
    
     % SIMULACION (TRANSIENTE)
 [~,BIOMASA_FINAL_S1,  poliSuma_S1, plantasSuma_S1,interTotales_S1,...
    maxEnlacesPol_S1, maxEnlacesPlanta_S1, visitPoli_S1, visitPlant_S1] = coreSimPrueba(datos_simulacion_final, ...
 replica, OUT, PARAMS_NAT);


path_Files = '/home/rosario/Desktop/resu_parcial/';

ex_et = strcat([path_Files,'Ex_eiaa_', num2str(num_net), '.mat']); 
       
load(ex_et); 

[riqPoli_S1, riqPlan_S1] = riquezaok(nPolinizadores, nPlantas , lsExtinct);


% 8 dic, esto debería redefinir BIOMASA_FINAL_S1 y si es diferente de 0
% entonces el ciclo while se corta y sigue a las siguientes simus 
%VAMOS A CREAR UNA MATRIZ CON TODOS ESTOS DATOS Y DESPUES AL FINAL LA
%GUARDAMOS COMO .CSV
DataSim1(num_net,1:8) = [num_net, riqPoli_S1, riqPlan_S1, poliSuma_S1, plantasSuma_S1,interTotales_S1...
    maxEnlacesPol_S1, maxEnlacesPlanta_S1]  ;
%11 enero, las frecuencias de salida ahora estan en otra matrix y en otro
%csv
%visitPoli_S1, visitPlant_S1];
fVis_pol_Sim1(num_net, :) = visitPoli_S1;
fVis_Pla_Sim1(num_net, :) = visitPlant_S1;

end        
% check point: safe point (no bugs) %12 dic ok
%% SIMULACION 2 continuación sin la inclusión de una especie. 
display("down sim 1")
%% Para esta segunda simulacion tenemos que cargar la red final (resultado de la primera)
path_Files = '/home/rosario/Desktop/resu_parcial/';

file_Network = strcat([path_Files,'Redes_eiaa_', num2str(num_net), '.mat']); 
       
load(file_Network);  

%COPIAMOS LA RED FINAL PARA EVITAR ALTERACIONES. 

Network_final_control = Redes.Final(:,:); %copiamos toda la red porque tengo miedo de que 
[polinizadores2,plantas2] = size(Network_final_control);

RED_INV_control = Network_final_control(:,:);

BIOMASA_S2_control =  BIOMASA_FINAL_S1(:,:);
Fenofases_control = Fenofases(:,:);

PARAMS_COMPLETOS_control = PARAMS_NAT(:,:); 

STEPS2 = 5000; % Cien anhos aprox. 
SIM = 1;
SHIFT = 0;

OUT = true; 
OUTPUT_DIR = '/home/rosario/Desktop/resu_parcial/';

ID_SIM_2_control = ['_ecaa_', num2str(num_net)]; %etapa final C control

datos_simulacion_final.id_sim = ID_SIM_2_control;
datos_simulacion_final.biomasas = BIOMASA_S2_control;
datos_simulacion_final.steps = STEPS2;
datos_simulacion_final.red = RED_INV_control;
datos_simulacion_final.fenos = Fenofases_control;
datos_simulacion_final.sim = SIM;
datos_simulacion_final.shift = SHIFT;
datos_simulacion_final.output_dir = OUTPUT_DIR;
%extraemos los datos
[~,~,poliSuma_S2, plantasSuma_S2,interTotales_S2, maxEnlacesPol_S2, maxEnlacesPlanta_S2,...
     visitPoli_S2, visitPlant_S2] = coreSimPrueba(datos_simulacion_final, replica, OUT, PARAMS_COMPLETOS_control);

ex_et_2 = strcat([path_Files,'Ex_ecaa_', num2str(num_net), '.mat']); 
       
load(ex_et_2); 

[riqPoli_S2, riqPlan_S2] = riquezaok(polinizadores2, plantas2 , lsExtinct);

%guardamos la simu
dataSim2(num_net,:) = [num_net,riqPoli_S2, riqPlan_S2 ,poliSuma_S2, plantasSuma_S2,interTotales_S2, maxEnlacesPol_S2, maxEnlacesPlanta_S2];...
%11 enero, las frecuencias de salida ahora estan en otra matrix y en otro
%csv
% 11 ENERO, HAY PROBLEMAS CON LAS DIMENSIONES DE LOS DOS LADOS DE LA MATRIZ
% CON LAS FRECUENCIAS, ASÍ QUE VAMOS A VER SI EL PROBLEMA ES QUE NO ESTAN
% TRANSPUESTAS


fVis_pol_Sim2(num_net, :) = visitPoli_S2;

fVis_Pla_Sim2(num_net,:) = visitPlant_S2;

display("down sim 2")
%% SIMULACION 3 continuación con la inclusión de una especie 11 ENE NO AGRE
%% Para esta segunda simulacion tenemos que cargar la red final (resultado de la primera)
file_Network = strcat([path_Files,'Redes_eiaa_', num2str(num_net), '.mat']); 
       
load(file_Network); 

%hacemos un loop hasta que se establezca la especie (por lo menos hasta el
%año 50 
abuinv = 0 ;
while  abuinv == 0

Network_final = Redes.Final(:,:); %copiamos toda la red porque tengo miedo de que 

[m,p] = size(Network_final);

%% DETERMINAMOS EL GRADO

grados = zeros(m,1); %para mejorar la eficiencia del script se crea una matriz del tamanho del num de iteracciones
for grado = 1:m
    grados(grado) = sum(Network_final(grado,:));    %13 nov modif para que use net final y no net original
end 
grados = grados(grados ~= 0); 
if grados == 0  %en el supuesto caso de que no haya interacciones pero si biomasa de plantas o algo así
    grado_sp_inv =3 ;
else
grado_sp_inv = max(grados) ;% modif 28 novceil(mean(grados))+3; %porque puede dar un numero decimal, lo redondeamos al mas bajo, (porque si)  
%17NOV arreglado para que el poli invasor tenga un grado de por lo menos 5
if grado_sp_inv < 5
    grado_sp_inv = 5;
end 
end 

grados_plantas = zeros(p,1); 
for grado = 1:p
    grados_plantas(grado) = sum(Network_final(:, grado)); %13 nov modif
end

maximos = zeros(1, grado_sp_inv);

for maximo = 1: grado_sp_inv
  
  [~, pos] = max(grados_plantas); 
  maximos(maximo) = pos;
  grados_plantas(pos) = NaN;
end 

nuevaSp = zeros(1,p); 

for inter = 1:p 
    if ismember(inter, maximos)
        nuevaSp(inter) =1; 
    end 
end
for pl = 1:p 
    nicho = find(nuevaSp==1);
end
% lo agrego a la red 
RED_INV = vertcat(Network_final,nuevaSp);
%copiamos la red para la simulacion con la especie agresiva (por las dudas)
RED_INV_AGRE = RED_INV(:,:);
%copiamos para agresivo sin competencia 
RED_INV_ASC = RED_INV(:,:);
%% Ahora modificamos la biomasa  (solo le tenemos que agregar la biomasa de una especie
%de polinizadores 

%modif 19 nov
%cargamos el vector que contiene los pesos posibles en mm³ de los
%polinizadores de una red. 

pesos_posibles =csvread('/home/rosario/PlantPollinatorModel/peso_insectos_ok', 0, 1);

%19o de todos esos pesos posibles, elegimos uno,  
p_s = pesos_posibles(randi(numel(pesos_posibles)));

%primero, vamos a SUPONER que la densidad del cuerpo del insecto es la
%densidad del agua, por lo tanto 1 mm³ es equivalente a 0.001. 
peso_seleccionado = p_s * 0.001;

%ahora tenemos el peso en GRAMOS. Esto es el peso de UN SOLO INSECTO
%y la densidad que se calcula acá es en g/m². 
% supongo, como es una especie nueva en la comunidad, hay 10 insectos por m²

biomasa_adultos = 3 * peso_seleccionado ;

%19o ahora creamos el vector que vamos a concatenar a las biomasas de las
%nativas 

Biomasa_insecto_invasor = [biomasa_adultos;0; 0];
BIOMASA_S2 =  vertcat(BIOMASA_FINAL_S1, Biomasa_insecto_invasor); 
%COPIO LA BIOMASA PARA LA SIMULACION CON LA AGRESIVA 
BIOMASA_AGRE = BIOMASA_S2(:,:);
BIOMASA_ASC = BIOMASA_S2(:,:);
%Checkpoint: safe (no bugs)

%% MODIFICAMOS LAS FENOFASES 
%12 Nov, como los insectos invasores no estan generando larvas (y pueden
%generarlas hipoteticamente si tienen los recursos florales) una hipotesis
%es que las fenofases producidas no se correspondan con las de las plantas
%que visitan (entonces no tienen alimento y solo se mueren sin producir
%nada). Como primer medida vamos a modificar la fenología del invasor para
%que se corresponda sí o sí con las plantas que consume 

%checkeamos con cuántas plantas interactúa
[~, plantas] = size(NET_ORIGINAL);

%Hacemos el array que contiene las fenologías de todas las plantas que 
%interactúan con los insectos, la hacemos con ceros y después la llenamos
plantas_invasora_bor = zeros(plantas,2);

for i= 1: length(nuevaSp)
    if nuevaSp(i) == 1
        plantas_invasora_bor(i,:) = V_PLANT(i, :);
    end 
end 
plantas_invasora_bor( all(~plantas_invasora_bor,2), : ) = [];
% Hacemos que la fenología sea "first in last out" 
comienzo_fen = min(plantas_invasora_bor(:,1));
fin_fen = max(plantas_invasora_bor(:,2));

fenos_inv = [comienzo_fen, fin_fen];

feno_originales = Fenofases.I ; 
fenofases_completas = vertcat(feno_originales, fenos_inv); 

Fenofases_final.I = fenofases_completas;
%11 ENE, COPIO PARA AGRESIVA
Fenofases_final_AGRE.I= fenofases_completas(:,:);
Fenofases_final_ASC.I = fenofases_completas(:,:);
Fenofases_final.P = V_PLANT;
%copio 
Fenofases_final_AGRE.P = V_PLANT(:,:);
Fenofases_final_ASC.P = V_PLANT(:,:);

%% PARAMETROS 

PARAMS_INVASORA = distrMBeta_Intr_sp(nuevaSp, basalValues, RANGE_VAR, SEED, NET_ORIGINAL);

PARAMS_COMPLETOS = parCompletos(PARAMS_BORRADOR, PARAMS_INVASORA); 
PARAMS_COMPLETOS_AGRE = PARAMS_COMPLETOS(:,:);

%% DATOS PARA LA SIMULACION

STEPS3 = 5000; % Cien anhos aprox.
SIM = 1;
SHIFT = 0;

OUT = true; 
  
OUTPUT_DIR = '/home/rosario/Desktop/resu_parcial/';

ID_SIM_2_COMP = ['_etcaa_', num2str(num_net)]; 

datos_simulacion_final.id_sim = ID_SIM_2_COMP;
datos_simulacion_final.biomasas = BIOMASA_S2;
datos_simulacion_final.steps = STEPS3;
datos_simulacion_final.red = RED_INV;
datos_simulacion_final.fenos = Fenofases_final;
datos_simulacion_final.sim = SIM;
datos_simulacion_final.shift = SHIFT;
datos_simulacion_final.output_dir = OUTPUT_DIR;

[~,~,poliSuma_S3, plantasSuma_S3,interTotales_S3, maxEnlacesPol_S3, maxEnlacesPlanta_S3,...
     visitPoli_S3, visitPlant_S3, abuinv] =... 
 coreSimPruebaDI2(datos_simulacion_final, replica, OUT, PARAMS_COMPLETOS);
%%vamos a calcular acá el solapameinto de nicho y demas

ex_c1 = strcat([path_Files,'Ex_etcaa_', num2str(num_net), '.mat']); 
       
load(ex_c1); 
[polinizadores3, plantas3] = size(RED_INV);
[riqPoli_S3, riqPlan_S3] = riquezaok(polinizadores3, plantas3, lsExtinct);
%calculamos el nicho para los tratamientos
[compPol,ncompPol, compPla, ncompPla] = solapNicho(RED_INV, lsExtinct);


%8 dic cambie para controlar la biomasa final de la comunidad 
dataSim3(num_net,:) = [ num_net, riqPoli_S3, riqPlan_S3 ,poliSuma_S3, ...
    plantasSuma_S3,interTotales_S3, maxEnlacesPol_S3, maxEnlacesPlanta_S3...
    compPol,ncompPol, compPla, ncompPla];

%11 enero, las frecuencias de salida ahora estan en otra matrix y en otro
%csv
%visitPoli_S1, visitPlant_S1];
fVis_pol_Sim3(num_net, :) = [visitPoli_S3];
fVis_Pla_Sim3(num_net, :) = [visitPlant_S3];
display("down sim 3")
end
%% SIMULACION 4 INTRODUCCION DE UN POLINIZADOR INVASOR CON COMPORTAMIENTO DE AGRESIVIDAD
%% DATOS PARA LA SIMULACION
abuinv2 = 0 ;
while  abuinv2 == 0
STEPS3 = 5000; % Cien anhos aprox.
SIM = 1;
SHIFT = 0;
OUT = true; 
OUTPUT_DIR = '/home/rosario/Desktop/resu_parcial/';

ID_SIM_2_AGRE = ['_etaaa_', num2str(num_net)]; %etapa tratamiento de agresividad

datos_simulacion_final.id_sim = ID_SIM_2_AGRE;
datos_simulacion_final.biomasas = BIOMASA_AGRE;
datos_simulacion_final.steps = STEPS3;
datos_simulacion_final.red = RED_INV_AGRE;
datos_simulacion_final.fenos = Fenofases_final_AGRE;
datos_simulacion_final.sim = SIM;
datos_simulacion_final.shift = SHIFT;
datos_simulacion_final.output_dir = OUTPUT_DIR;

[~,~,poliSuma_S4, plantasSuma_S4,interTotales_S4, maxEnlacesPol_S4, maxEnlacesPlanta_S4,...
     visitPoli_S4, visitPlant_S4,abuinv2] =... 
 coreSimPruebaDI(datos_simulacion_final, replica, OUT, PARAMS_COMPLETOS);
%%vamos a calcular acá el solapameinto de nicho y demas

ex_c1 = strcat([path_Files,'Ex_etaaa_', num2str(num_net), '.mat']); 
       
load(ex_c1); 
[polinizadores4, plantas4] = size(RED_INV_AGRE);
[riqPoli_S4, riqPlan_S4] = riquezaok(polinizadores4, plantas4, lsExtinct);
%calculamos el nicho para los tratamientos
[compPol4,ncompPol4, compPla4, ncompPla4] = solapNicho(RED_INV_AGRE, lsExtinct);


%8 dic cambie para controlar la biomasa final de la comunidad 
dataSim4(num_net,:) = [ num_net, riqPoli_S4, riqPlan_S4,poliSuma_S4, ...
    plantasSuma_S4,interTotales_S4, maxEnlacesPol_S4, maxEnlacesPlanta_S4...
    compPol4,ncompPol4, compPla4, ncompPla4];

%11 enero, las frecuencias de salida ahora estan en otra matrix y en otro
%csv
%visitPoli_S1, visitPlant_S1];
fVis_pol_Sim4(num_net, :) = [visitPoli_S4];
fVis_Pla_Sim4(num_net, :) = [visitPlant_S4];
display("down sim 4")
end
%% SIMULACION 5 INTRODUCCION DE UN POLINIZADOR AGRESIVO SIN COMPETENCIA 
abuinv3 = 0;
while  abuinv3 == 0

PARAMS_COMPLETOS_ASC = parCompletosASC(PARAMS_BORRADOR, PARAMS_INVASORA);

STEPS4 = 5000; % Cien anhos aprox.
SIM = 1;
SHIFT = 0;
OUT = true; 
OUTPUT_DIR = '/home/rosario/Desktop/resu_parcial/';

ID_SIM_2_ASC = ['_etascaa_', num2str(num_net)]; %etapa tratamiento de agresividad

datos_simulacion_final.id_sim = ID_SIM_2_ASC;
datos_simulacion_final.biomasas = BIOMASA_ASC;
datos_simulacion_final.steps = STEPS4;
datos_simulacion_final.red = RED_INV_ASC;
datos_simulacion_final.fenos = Fenofases_final_ASC;
datos_simulacion_final.sim = SIM;
datos_simulacion_final.shift = SHIFT;
datos_simulacion_final.output_dir = OUTPUT_DIR;

[~,~,poliSuma_S5, plantasSuma_S5,interTotales_S5, maxEnlacesPol_S5, maxEnlacesPlanta_S5,...
     visitPoli_S5, visitPlant_S5,abuinv3] =... 
 coreSimPruebaDI(datos_simulacion_final, replica, OUT, PARAMS_COMPLETOS_ASC);
%%vamos a calcular acá el solapameinto de nicho y demas

ex_c1 = strcat([path_Files,'Ex_etascaa_', num2str(num_net), '.mat']); 
       
load(ex_c1); 
[polinizadores5, plantas5] = size(RED_INV_ASC);
[riqPoli_S5, riqPlan_S5] = riquezaok(polinizadores5, plantas5, lsExtinct);
%calculamos el nicho para los tratamientos
[compPol5,ncompPol5, compPla5, ncompPla5] = solapNicho(RED_INV_ASC, lsExtinct);

%8 dic cambie para controlar la biomasa final de la comunidad 
dataSim5(num_net,:) = [ num_net, riqPoli_S5, riqPlan_S5,poliSuma_S5, ...
    plantasSuma_S5,interTotales_S5, maxEnlacesPol_S5, maxEnlacesPlanta_S5...
    compPol5,ncompPol5, compPla5, ncompPla5];

%11 enero, las frecuencias de salida ahora estan en otra matrix y en otro
%csv
%visitPoli_S1, visitPlant_S1];
fVis_pol_Sim5(num_net, :) = [visitPoli_S5];
fVis_Pla_Sim5(num_net, :) = [visitPlant_S5];

display("c'est fini, <3, I'm proud of you kiddo.")
end
end 
dlmwrite('/home/rosario/Desktop/resu_parcial/ST_aa.csv',DataSim1,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/FV_PO_T_aa.csv',fVis_pol_Sim1,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/FV_PL_T_aa.csv',fVis_Pla_Sim1,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/SC_aa.csv',dataSim2,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/FV_PO_C_aa.csv',fVis_pol_Sim2,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/FV_PL_C_aa.csv',fVis_Pla_Sim2,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/SCO_aa.csv',dataSim3,'precision','%i')  ;
dlmwrite('/home/rosario/Desktop/resu_parcial/FV_PO_CO_aa.csv',fVis_pol_Sim3,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/FV_PL_CO_aa.csv',fVis_Pla_Sim3,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/SAG_aa.csv',dataSim4,'precision','%i')  ;
dlmwrite('/home/rosario/Desktop/resu_parcial/FV_PO_AG_aa.csv',fVis_pol_Sim4,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/FV_PL_AG_aa.csv',fVis_Pla_Sim4,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/SASC_aa.csv',dataSim5,'precision','%i')  ;
dlmwrite('/home/rosario/Desktop/resu_parcial/FV_PO_ASC_aa.csv',fVis_pol_Sim5,'precision','%i');
dlmwrite('/home/rosario/Desktop/resu_parcial/FV_PL_ASC_aa.csv',fVis_Pla_Sim5,'precision','%i');

toc
%end 
 
