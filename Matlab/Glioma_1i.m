%Paciente 1
%Glioma_1i
%Corrijo ubicacion del tumor y nivel de deteccion del tumor para
%diagnostico y letalidad
dt = 0.1; %dias
h = 1; %mm, malla de 181x217x181 mm (18x22x18 cm)
ii = 181; %En la imagen final se ven como columnas x) 
jj = 217; %En la imagen final se ven como filas y)
kk = 181; % (slices,z)
C0 = 0; %concentracion inicial de celulas tumorales (por mm3)
C_max = 1e8; %concentracion maxima de celulas (por mm3)
C_mig = 1e7; %concentracion de celulas a la que comienza la migracion (por mm3)
IM = 5; %indice mitotico
nn = 5000; %corrida de 500 dias
max_iter = 100;

%Hago mapa de coeficientes
disp(datestr(now,'HH:MM AM'));
cerebro=importdata('Cerebro.mat');
talairach=importdata('Talairach_conv.mat');
for k=1:kk
    for j=1:jj
        for i=1:ii
            if ((cerebro(i,j,k)>110) && (cerebro(i,j,k)<=225)) %sustancia blanca
                p(i,j,k) = 0.107; %proliferacion neta, 1/dia
                switch talairach(i,j,k)
                    case{508,509,510,514,515,651,652,657,694,695,711,712,713,714,715} %cuerpo calloso
                        D(i,j,k)=0.306; %un 20% mas
                    case{273,359,361,624,633,638,639,641} %tracto optico
                        D(i,j,k)=0.306; %un 20% mas
                    case{5,6,71,72,215,216,440,459,502,503,574,576} %tallo cerebral, medula, protuberancia, mesensefalo
                        D(i,j,k)=0.204; %un 20% menos
                    otherwise
                        D(i,j,k)=0.255; %migracion neta, mm2/dia
                end;
            else if ((cerebro(i,j,k)>=75) && (cerebro(i,j,k)<=110)) %sustancia gris
                    p(i,j,k)=0.107; 
                    switch talairach(i,j,k)
                        case{1,2,3,4,61,62,63,64,65,66,67,68,80,81,82,83,84,85,86,87,88,89,114,115,116,117,118,119,120,157,158,159,160,161,162,163,164,165,211,212,321,322,564} %cerebelo
                            D(i,j,k)=0;
                        case{349,350,351,435,436,464,465,492,493,494,495,496,571,598,601,692,693,379,382,596,597} %nucleo estriado (caudado + putamen)
                            D(i,j,k)=0.0408; %un 20% menos
                        case{452,455,456,457,593,594,595,341,342,453,454,379,382,452,455,456,457,501,640,181,182,183,219,356,358,449,450,513} %globo palido, sustancia nigra, nucleo subtalamico, nucleo lentiforme, amigdala, claustrum
                            D(i,j,k)=0.0408; %un 20% menos
                        case{5,6,71,72,215,216,341,342,343,344,354,355,437,438,440,453,454,498,499,574} %tallo cerebral, medula, protuberancia, mesensefalo
                            D(i,j,k)=0.0408; %un 20% menos
                        otherwise
                            D(i,j,k)=0.051;
                    end;
                 else
                    p(i,j,k) = 0;
                    D(i,j,k) = 0;
                 end;
            end;
        end;
    end;
end;
clear cerebro;

%Carga vectores
x(1)=0;
for i=2:ii
    x(i)=x(i-1)+h;
end;
y(1)=0;
for j=2:jj
    y(j)=y(j-1)+h;
end;
z(1)=0;
for k=2:kk
    z(k)=z(k-1)+h;
end;

%Condiciones iniciales
%Origen del tumor en:
io = 32;
jo = 98;
ko = 85;
for k=1:kk
   for j=1:jj
       for i=1:ii
       C(i,j,k) = C0;
       end;
   end;
end;
C = ones(ii,jj,kk)*C0;
C(io,jo,ko) = 1;
%C(i,j,k) = C(i,j,k);
dia = 1;
migracion = 0; % 0 para tumor benigno y 1 para tumor maligno
diagnostico = 3188; %Para diagnostico a un diametro de 18.26 mm

%Optimizaciones:
for k=2:kk-1
    for j=2:jj-1
        for i=2:ii-1
            P_optimizado(i,j,k) = dt * p(i,j,k);
            M_optimizado(i,j,k) = dt * D(i,j,k) / (h*h);
        end;
    end;
end;

info = fopen('Info_1i.txt','a+');
fprintf(info,'Simulacion i del paciente 1 \n');
status = fclose(info);

%%%%%%%%%%%%%%%%%%%%%%%
% Iteracion temporal
%%%%%%%%%%%%%%%%%%%%%%%
for n=1:nn
    
  %Calculo del volumen tumoral y chequeo de areas de Brodmann
   cantidad1 = 0; %para deteccion de celulas tumorales
   cantidad2 = 0; %para diagnostico
   cantidad3 = 0; %para letalidad
   B_T=[1200,4371,7001,8085,4237,39416,24752,15300,19996,24206,15580,0,13017,0,0,0,4793,24266,25120,11476,14008,12039,2720,8317,2592,0,309,2125,836,4030,11256,9586,196,1809,1274,3024,8910,10291,9508,22172,3032,2492,1676,3180,4156,7903,11145]; %vector con areas de Brodman totales del Talairach
   B=zeros(1,47); %vector que guarda areas de Brodmann absolutas
   B_R=zeros(1,47); %vector que guarda areas de Brodmann relativizadas
   for k=1:kk
        for j=1:jj
            for i=1:ii
                if (C(i,j,k) >= 1)
                    cantidad1=cantidad1+1; 
                end;
                if (C(i,j,k) >= 1e7)
                    cantidad2=cantidad2+1;
                    if ((i>80)&&(i<102)&&(j>90)&&(j<117)&&(k>61)&&(k<71)) %area del foramen del tentorio
                        cantidad3=cantidad3+2;
                    else
                        switch talairach(i,j,k)
                            case{1,2,3,4,61,62,63,64,65,66,67,68,80,81,82,83,84,85,86,87,88,89,114,115,116,117,118,119,120,157,158,159,160,161,162,163,164,165,211,212,321,322,564} %cerebelo
                                cantidad3=cantidad3+2;
                            case{5,6,71,72,215,216,440,459,502,503,574,576} %tallo cerebral, medula, protuberancia, mesensefalo
                                cantidad3=cantidad3+2;
                            otherwise
                                cantidad3=cantidad3+1;
                        end;
                    end;
                    if (D(i,j,k)==0.051) %si estoy en corteza (de cerebro o cerebelo)
                        switch talairach(i,j,k) %busco areas de Brodmann
                            case{893,894}
                                B(1)=B(1)+1;
                            case{891,892,978}
                                B(2)=B(2)+1;
                            case{895,896,897,898,938,939,940,1019,1088,1091,1092}
                                B(3)=B(3)+1;
                            case{804,805,806,808,835,836,1003,1089,1090,1093,1099,1100,1104,1105}
                                B(4)=B(4)+1;
                            case{1064,1067,1068,1070,1083,1086,1095,1096,1097,1098}
                                B(5)=B(5)+1;
                            case{682,683,809,810,1009,1010,1011,1012,1022,1023,1024,1075,1076,1077,1078,1094}
                                B(6)=B(6)+1;
                            case{952,953,968,969,971,972,1039,1040,1043,1044,1065,1066,1080,1081,1087,1101,1102}
                                B(7)=B(7)+1;
                            case{1025,1026,1027,1028,1033,1034,1035,1036}
                                B(8)=B(8)+1;
                            case{813,814,815,816,841,899,900,901,902,903,943,944,981,982,983,984,1004}
                                B(9)=B(9)+1;
                            case{405,406,407,408,410,411,412,413,471,483,484,516,517,754,755}
                                B(10)=B(10)+1;
                            case{102,103,109,142,145,152,153,191,192,197,198,235,236,300,301}
                                B(11)=B(11)+1;
                            case{288,292,294,295,296,297,365,366,367,368,369,376,400,444,447,589,590,733,734,735}
                                B(13)=B(13)+1;
                            case{242,245,303,305,475,477,540}
                                B(17)=B(17)+1;
                            case{203,205,247,249,250,254,304,306,307,474,478,487,488,537,631,632,758,760,826,827,857,858}
                                B(18)=B(18)+1;
                            case{258,259,263,264,265,266,315,316,414,415,417,418,421,423,425,426,427,428,479,480,566,568,699,700,818,819,820,843,844,848,849,986,987}
                                B(19)=B(19)+1;
                            case{8,13,15,16,18,19,23,26,29,46,53,123,124,135,167,171,172,217,218}
                                B(20)=B(20)+1;
                            case{77,78,79,97,275,276,286,352}
                                B(21)=B(21)+1;
                            case{442,443,445,446,458,460,587,588,634,764,765,766,767}
                                B(22)=B(22)+1;
                            case{706,707,762,763,829,830,934,935}
                                B(23)=B(23)+1;
                            case{469,470,941,942,1020,1021}
                                B(24)=B(24)+1;
                            case{226,227,234,237,281,282,380,381}
                                B(25)=B(25)+1;
                            case{432,433}
                                B(27)=B(27)+1;
                            case{94,95,169,170,345}
                                B(28)=B(28)+1;
                            case{670,671}
                                B(29)=B(29)+1;
                            case{429,430,431,481,482,660,661,663,665,669}
                                B(30)=B(30)+1;
                            case{756,757,823,824,856,859,863,905,909,915,916,920,923,932,994,995,998,1016,1017,1045,1049,1050,1053}
                                B(31)=B(31)+1;
                            case{395,396,610,612,945,946,947,948,1006,1008,1029,1030,1031,1032}
                                B(32)=B(32)+1;
                            case{838,839}
                                B(33)=B(33)+1;
                            case{177,178,179,180,279,280,283,285,287}
                                B(34)=B(34)+1;
                            case{140,141,272}
                                B(35)=B(35)+1;
                            case{75,76,129,130,136,166}
                                B(36)=B(36)+1;
                            case{208,209,213,214,267,268,269,270,271,323,324,328,329,330,331,332,333,334,420,422,558,560,627,628}
                                B(37)=B(37)+1;
                            case{33,34,35,36,38,41,44,45,59,60,101,221,222}
                                B(38)=B(38)+1;
                            case{702,703,704,705,708,709,865,906,907,908,910,917,918,958,959,967,970,974,975,989,990,991,1014,1015}
                                B(39)=B(39)+1;
                            case{774,778,781,785,869,872,877,878,881,885,927,928,992,993,999,1000,1061,1062,1082}
                                B(40)=B(40)+1;
                            case{672,673,719,720,721,722,723,768,770,771}
                                B(41)=B(41)+1;
                            case{674,675,727,729,741,742}
                                B(42)=B(42)+1;
                            case{738,739,740,792,793,797}
                                B(43)=B(43)+1;
                            case{689,690,691,748,749}
                                B(44)=B(44)+1;
                            case{602,603,653,655,656}
                                B(45)=B(45)+1;
                            case{607,608,658,696,697}
                                B(46)=B(46)+1;
                            case{146,147,185,188,298,299,302,392,393,399,401,402,409,466,467,468,511,604}
                                B(47)=B(47)+1;
                            otherwise
                                %no modifica nada
                          end;
                     end;
                end;
            end;
        end;
   end;
   
    %%%%%%%%%%%%%%%%%%%%%%%
    %iteracion de convergencia (flotacion de punto fijo)
    %%%%%%%%%%%%%%%%%%%%%%%
        
    if (dia <= 1500)
        max_error = 1;
    else
        max_error = 10;
    end;
    iter=0;
    error = 1000;
    C_k2 = C;
    while((iter < max_iter) && (error > max_error))
        C_k1 = C_k2;
               
       %Calcular dominio
       for k=2:kk-1
           for j=2:jj-1
               for i=2:ii-1
                if (C(io,jo,ko) < C_mig) %proliferacion
                    M(i,j,k) = 0;
                else %proliferacion y migracion
                    M(i,j,k) = M_optimizado(i,j,k) * (C_k1(i+1,j,k)+C_k1(i-1,j,k)+C_k1(i,j+1,k)+C_k1(i,j-1,k)+C_k1(i,j,k+1)+C_k1(i,j,k-1)-6*C_k1(i,j,k));
                end;
                P(i,j,k) = P_optimizado(i,j,k) * C_k1(i,j,k) * (1 - (C_k1(i,j,k) / C_max));
                C_k2(i,j,k) = C(i,j,k) + P(i,j,k) + M(i,j,k);
                if C_k2(i,j,k) > C_max
                    C_k2(i,j,k) = C_max;
                end;
                if C_k2(i,j,k) < 0.00001
                    C_k2(i,j,k) = 0;
                end;
               end;
           end;
       end;
  
       %% Calcular error y actualizar
       error = max(max(max(abs(C_k1-C_k2))));
       iter=iter+1; 
       if (mod(iter,10)==0)
           matriz_error=abs(C_k1-C_k2);
       end;
    end;
    
    %% Actualizar malla
    
    C = C_k2;
   
    %% Grabar datos
    if (mod(n,100)==0)
       save(['DATOS_1i_' num2str(dia) '.mat' ],'C')
    end
    
    %% Grabar info e informar por pantalla:
    if (C(io,jo,ko) > C_mig && migracion == 0)
            disp('comienza migracion');
            migracion = 1; 
    end;
    if error > 10
            disp(['error = ' num2str(error)]);
    end;
    if (mod(n,10)==0)
        disp(datestr(now,'HH:MM AM'));
        disp(['dia ' num2str(dia)]);
        info = fopen('Info_1i.txt','a+');
        fprintf(info,datestr(now,'HH:MM AM'));
        fprintf(info,[' Dia ' num2str(dia)]);
        fprintf(info,[' error = ' num2str(error)]);
        if migracion == 1
            fprintf(info,' en migracion');
        end;
        if (cantidad2 >= diagnostico)
            fprintf(info,' diagnosticado');
        end;
        if (cantidad3 >= 179594) %Para area letal (esfera de 70 mm de diametro) ponderada por areas vitales
            fprintf(info,' muerte del paciente');
        end;
        for a=1:47
            if B(a)>=1 %si esa area de Brodmann no es nula
                B_R(a)=B(a)*100/B_T(a); %la relativizo con respecto a la total del Talairach
                fprintf(info,[' Brodmann ' num2str(a) ' = ' num2str(eval(['B(' num2str(a) ')'])) ' porcentaje: ' num2str(eval(['B_R(' num2str(a) ')']))]);
            end;
        end;
        fprintf(info,'\n');
        status = fclose(info);
        dia = dia + 1;  
    end;
end

disp('Proceso Terminado!');

