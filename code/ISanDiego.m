% clear;
% clc;

%Llamada de datos

%Nombres de los 677 archivos
load 'C:\Users\Gabriela\Downloads\datoscal\filesname.mat' ;

% Numero de los dias que cuentan con datos de viento.
load 'C:\Users\Gabriela\Downloads\datoscal\DiasSd.mat' ;

% Matriz donde se guardan los datos
resultados= zeros(19,length(B));


%contador
n = 1;

 while n< length(B)+1

    
    % Obtenemos los datos del archivo
    file = filename(B(n));
    direccion = strcat('C:\Users\Gabriela\Downloads\datoscal\',file);
    load(direccion);

    % Creacion de los vectores
    time = time_day(:,1);                     % tiempo
    rglbt = GHI_day(:,1);                     % GHI (Radiacion global horizontal)
    rhum = data_day(:,3);                     % Humedad relativa
    temperatura = data_day(:,2);              %Temperatura

    
    % San Diego Coordenadas (32.879923 , -117.232848)

    Location.latitude = 32.879923;
    Location.longitude = -117.232848;
    Location.altitude = 0;

    %UTCoffset depende de la diferencia horaria -7 o -8 depende del horario de verano o invierno
    UTCoffset = -8; % Depende de la diferencia horaria -8 presente en todo el vector

    win_length = 20; % considera un intervalo de 20 
    sample_interval = 5; % datos cada  minutos en el vector

    ti = datetime(time,'ConvertFrom','datenum');
    ti.TimeZone = 'America/Los_Angeles';

    % Deteccion de tiempos cielos despejados
    [CD] = pvl_detect_clear_times(rglbt, time, UTCoffset, Location, win_length, sample_interval);
    
    sti =datetime(time,'ConvertFrom','datenum');
    sti.TimeZone = 'America/Los_Angeles';
    [pks,locs] = findpeaks(rglbt,sti);

    % Sacar los peaks
    [peaksghi,tiemppeak] = findpeaks(rglbt,ti);

    % Radiacion en los tiempos donde detecta cielo despejado
    radiacion = rglbt(CD);
    
    % Tiempos donde detecta cielo despejado
    tiempocd = datenum(sti(CD));
    DN = tiempocd.';
    Time = pvl_maketimestruct(DN, UTCoffset);
    
    % Vector de todos los indices de turbidez que se prueban desde 0 con aumento de 0.001 hasta 8
    vectorCS = 0:0.001:8;
    
    % Matriz con todos los cielos despejados son los datos que hay cada dia
    %estudiado
    VectorCD = zeros(length(radiacion),length(vectorCS));
    
    % Guarda la comparacion
    VectorComparacion = 0:length(vectorCS);
    
    % Contador
    indiceCS = 1;

    while indiceCS <= length(vectorCS)        
    
            % Entrega el Cielo Despejado para la turbidez estudiada
            [ClearSkyGHI]= pvl_clearsky_ineichen(Time, Location,vectorCS(indiceCS));

            VectorCD(:,indiceCS) = ClearSkyGHI;
    
            % Entrega un numero positivo que indica la cercania de la
            % radiacion de cielo despejado medida con la teorica en los momentos de cielo despejado. Menos es
            % mejor
            VectorComparacion(indiceCS) = abs(nanmean(ClearSkyGHI./radiacion-1));
    
            indiceCS = indiceCS + 1;
    end
    
    
    DN = time.';
    Time = pvl_maketimestruct(DN,UTCoffset);
    
    %Buscamos el indice del vector decomparacion menor 
    k = find(min(VectorComparacion.')==VectorComparacion.');

    % Turbidez
    turb =vectorCS(k);
    
    [ClearSkyGHI]= pvl_clearsky_ineichen(Time, Location,turb);

    %Indice de cielo despejado
    indice = rglbt./ClearSkyGHI;
    
    
    %Iniciamos las variables
    tfinalruptura = NaN;            % tiempo de ruptura final
    rfr = NaN;                      % radiación final ruptura

    
    tinicial = 1 ;                  % tiempo de amanecer
    
    tinicialruptura = NaN;          % tiempo de ruptura inicial
    rir= NaN;                       % radiación inicio ruptura
    kruptura = NaN;                     % indice inicial ruptura
 

    %Contador para determinar el amanecer
    f = 1 ;
    
    while f <(length(rglbt)+1)
    
        % El amanecer corresponde al tiempo anterior del registro de un
        % radiacion mayor a 5 w/m^2
        if rglbt(f) > 5     
            tinicial = ti(f-1);
            break;      
        end
            f=f+1;
    end
    
    
    % Contador para determinar el inicio de la ruptura
    p = f;
    
    while p <=(length(indice))
    
         if isnan(indice(p))==0  & isnan(find( locs==ti(p)))==0
                  
                    %Criterio para el tiempo inicial es que sea un peak y que su
                    %indice sea menor a 1.2 y que sea anterior a un cielo
                    %despejado
         
                    if  indice(p)< 1.2 & isnan(find(locs==ti(p)))==0  & datenum(ti(p))<=min(datenum(sti(CD)))
                        tinicialruptura = ti(p);
                        rir = rglbt(p);
                        kruptura = indice(p);
   
                       
                        break;

                    end
         end

         p=p+1;
    end
    
    
    % Inicio del contador para determinar el tiempo de finalizacion de la
    % ruptura

    q = p;      %parte en p, es decir, el final de la ruptura es posterior al inicio.
    
    
    while q <=(length(indice)-6) && isnan(datenum(tinicialruptura))==0 
    
        if isnan(indice(q))==0 
    
          
                %Criterio para el tiempo final es el primer tiempo que tenga un
                %cielo despejado de media hora.
               if isnan(find( sti(CD)==ti(q)))==0 &  isnan(find( sti(CD)==ti(q+1)))==0 & isnan(find( sti(CD)==ti(q+2)))==0 & isnan(find( sti(CD)==ti(q+3)))==0 & datenum(ti(p))<=min(datenum(sti(CD)))
                
                    tfinalruptura = ti(q);
                    rfr = rglbt(q);

                    break;
               end
    
    
        end
        q=q+1;
    end


    promkrup = NaN;              % Promedio del indice en la ruptura
    
    if isempty(indice)==0 &&  p <288 && q<288
        krango = indice(p:q);
        promkrup = nanmean(krango);
    
    end

    deltaruptura = NaN;

    if isnan(datenum(tfinalruptura))==0 & isnan(datenum(tinicialruptura))==0
        deltaruptura = tfinalruptura - tinicialruptura;
    end

    % Iniciacion de variables atmosfericas
    %radiacion global
    rglbprom = NaN;
    rglbmax = NaN;
    rglbmin = NaN;
    stdrglb = NaN;
    
    %temperatura
    tempromrup = NaN;
    tempminrup = NaN;
    tempmaxrup = NaN;
    stdtemp = NaN;
   
    %humedad relativa
    hrelpromrup = NaN;
    hrelmaxrup = NaN;
    hrelminrup = NaN;
    stdhrel = NaN;

    if p<289 && q<289 && p<q

        vectorrglbtruptura = rglbt(p:q);
        
        rglbprom = nanmean(vectorrglbtruptura);
        rglbmax = nanmax(vectorrglbtruptura);
        rglbmin = nanmin(vectorrglbtruptura);
        stdrglb = std(vectorrglbtruptura,'omitnan');
    
        vectortemruptura = temperatura(p:q);
        
        tempromrup = nanmean(vectortemruptura);
        tempminrup = nanmin(vectortemruptura);
        tempmaxrup = nanmax(vectortemruptura);
        stdtemp = std(vectortemruptura,'omitnan');
        
        
        vectorhrelruptura = rhum(p:q);
        
        hrelpromrup = nanmean(vectorhrelruptura);
        hrelmaxrup = nanmax(vectorhrelruptura);
        hrelminrup = nanmin(vectorhrelruptura);
        stdhrel = std(vectorhrelruptura,'omitnan');

    end

    resultados(1,n)= minute(tinicialruptura)+hour(tinicialruptura)*60; 
    resultados(2,n)= minute(tfinalruptura)+hour(tfinalruptura)*60; 

    resultados(3,n)=NaN;

    if isnan(deltaruptura)==0
        resultados(3,n)= minutes(deltaruptura);
    end

    resultados(4,n)= NaN;

    if isempty(turb)==0
       resultados(4,n)= turb;
    end

    resultados(5,n)= kruptura; 
    resultados(6,n)= minute(tinicial)+hour(tinicial)*60; 
    resultados(7,n)= rglbprom; 
    resultados(8,n)= rglbmax; 
    resultados(9,n)= rglbmin; 
    resultados(10,n)= stdrglb; 
    resultados(11,n)= promkrup; 
    resultados(12,n)= tempromrup; 
    resultados(13,n)= tempminrup;
    resultados(14,n)= tempmaxrup; 
    resultados(15,n)= stdtemp; 
    resultados(16,n)= hrelpromrup; 
    resultados(17,n)= hrelmaxrup; 
    resultados(18,n)= hrelminrup; 
    resultados(19,n)= stdhrel; 
    
    tiemposruptura = [tinicialruptura,tfinalruptura];
    radiacionrupt = [rir, rfr];

    diatitulo = " (" +num2str(day(time(1,1)))+"/" + num2str(month(time(1,1)))+ "/" + num2str(year(time(1,1))) + ")";
    diatitulo1 = num2str(day(time(1,1)))+"-" + num2str(month(time(1,1)))+ "-" + num2str(year(time(1,1))) + ".png";


    figure('Visible', 'off');
    hold on
    plot(ti, rglbt);
    plot(ti,ClearSkyGHI);
    hold off
    title(strcat(strcat('Irradianza global ',diatitulo),', San Diego'));
    xlabel('Hora (America/Los Angeles)');
    ylabel('GHI (W/m^2)');
    legend('GHI','GHI_{d}');


    nombredoc = strcat('img_Ant_',diatitulo1);
    imagen = gcf;
    exportgraphics(imagen,nombredoc);

    figure('Visible', 'off');
    hold on
    title(strcat(strcat('GHI y Tiempos de fragmentación ',diatitulo), ', San Diego'));
    plot(ti,rglbt,'-r');
    plot(tiemposruptura,radiacionrupt,'*b');
    xlabel('Time (America/Los Angeles)');
    ylabel('GHI (W/m^2)');
    legend('GHI','Tiempos de inicio'+ string(newline) + ' y término');
    hold off

    imagen = gcf;
    nombredoc = strcat('Tiempos_Ant_',diatitulo1);
    exportgraphics(imagen,nombredoc);


    figure('Visible', 'off');
    hold all
    plot(sti, rglbt, 'b-');
    plot(sti(CD), rglbt(CD), 'r.');
    xlabel('Hora (America/Los Angeles)');
    ylabel('GHI (W/m^2)');
    title(strcat(strcat('Momentos de Cielo Despejado ',diatitulo), ', San Diego'));
    %datetick('x','mm/dd','KeepTicks');
    legend('GHI','Momentos CD');

    imagen = gcf;
    nombredoc = strcat('Desp_Ant_', diatitulo1);
    exportgraphics(imagen,nombredoc);


    n = n+1 ;
 end

 %%
 load filenamewind.mat; % Nombre de los archivos de viento de cada dia

 resultadosV= zeros(13,length(filenamew));

 %%

%contador
n = 1;

while n< length(filenamew)+1
    
    file = filenamew(n);
    direccion = strcat('C:\Users\Gabriela\Documents\DatosVientoCalifornia\Cada 5 min\',file);
    T = readtable(direccion);

    time = T.Var1;
    velvt = T.wind_speed;
    dirvt = T.wind_dir;

    % Iniciamos las variables estudiadas

    % Vector de la velocidad en u (Predominante en San Diego)
    uA=-velvt.*sin(deg2rad(dirvt)); 
    % Vector de la velocidad en v
    vA=-velvt.*cos(deg2rad(dirvt));
    
    % Vector de Magnitud
    V = sqrt(uA.^2 + vA.^2);

    % Iniciamos las variables.
    inicioviento = 1 ;
    maginicioviento = 1;
    
    %contador
    contadorv = (resultados(6,n)+60*8)/5; %parte desde el amanecer
    
    % Viento suavizado de la magnitud
    Vientosuavizado = smooth(datenum(time),V,0.1,'rloess'); 
    
    while contadorv <=(length(V)-5)

            %Criterio para el tiempo inicial es que tenga un aumento por
            %los siguientes 30 minutos. 

            if  Vientosuavizado(contadorv)<Vientosuavizado(contadorv+1) && Vientosuavizado(contadorv+1)<Vientosuavizado(contadorv+2)  && Vientosuavizado(contadorv+2)<Vientosuavizado(contadorv+3)  && Vientosuavizado(contadorv+3)<Vientosuavizado(contadorv+4)

                  % Tiempo en el inicio de la brisa marina
                  inicioviento = time(contadorv);

                  % Magnitud en el inicio del viento
                  maginicioviento = V(contadorv); 

                  break;
            end


            contadorv=contadorv+1;
    end
    
    inicioviento = (hour(inicioviento))*60 + minute(inicioviento)-8*60;

    if contadorv>=(length(V))-5
        inicioviento = NaN;
        maginicioviento = NaN;

    end

    

    promvelv = nanmean(vA);
    promvelu = nanmean(uA);
    prommag = nanmean(V);
    maxmagnitudviento = nanmax(V);
    minmagnitudviento = nanmin(V);
 
    %Llamamos los indices de radiacion hacemos cambio a UTC sumandole los
    % minutos en 8 horas de diferencia.

    tiemporinicial= resultados(1,n)+60*8;
    tiemporfinal= resultados(2,n)+60*8;
 
    %contador
    p = 1;
    
    while p <length(V)

         % Calculamos los minutos para poder compararlos
         minIV = hour(time(p))*60 + minute(time(p));

         if tiemporinicial==minIV 

            break;        
         end
         p=p+1;
    end
    
    
    q = p; % parte en p, es decir, imponemos que el termino de la ruptura es poterior al inicio.
    
    
    while q <length(V)

        % Calculamos los minutos para poder compararlos
        minIF = hour(time(q))*60 + minute(time(q));
        
        if tiemporfinal==minIF
            break;      
        end

        q=q+1;
    end

    % Definimos los inidices en el caso que el tiempo de inicio o final
    % sean NaN
    
    maginicioruptura = NaN  ; 
    magfirup = NaN ;  

    promagruptura =  NaN ;
    maxmagruptura = NaN ;
    minmagruptura = NaN ;

    if  isnan(tiemporinicial)==0 && isnan(tiemporfinal)==0
        % Indices dependientes de la ruptura
    
   
    
        maginicioruptura = V(p) ; 
        magfirup = V(q); 
     
        % Vectores dependientes de la ruptura
        vectorruptura = time(p:q);
        Magrup = V(p:q);
    

        promagruptura =  nanmean(Magrup);
        maxmagruptura = nanmax(Magrup);
        minmagruptura = nanmin(Magrup);

    end

    % Se obtiene el tiempo de inicio del amanecer y se pasa a UTC
    tiempoF = resultados(7,n)+60*8;

    % Se aproxima el tiempo calculado de inicio del amanecer para evitar
    % los decimales

    if rem(tiempoF,10)==4 
           tiempoF = tiempoF+1;      
    end

    if rem(tiempoF,10)==9
           tiempoF = tiempoF+1;      
    end

    %Inicializacion de un contador para determinar el amanecer
    f = 1;

    while f <(length(time))

        % Se combierte el tiempo en minutos para realizar la comparacion
        TiempoIF = hour(time(f))*60 + minute(time(f));

        if tiempoF==TiempoIF
            break;      
        end
        f=f+1;
    end
    

    Vamanecer = V(f);
    
    % Copiar en la matriz de resultados

    resultadosV(1,n) = promvelv;
    resultadosV(2,n) = promvelu;
    resultadosV(3,n) = prommag;
    resultadosV(4,n) = maxmagnitudviento;
    resultadosV(5,n) = minmagnitudviento;
    resultadosV(6,n) = inicioviento;
    resultadosV(7,n) = maginicioviento;
    resultadosV(8,n) = Vamanecer;
    resultadosV(9,n) = promagruptura;
    resultadosV(10,n) = maxmagruptura;
    resultadosV(11,n) = minmagruptura;
    resultadosV(12,n) = maginicioruptura;
    resultadosV(13,n) = magfirup;
 

    
    n = n+1;
end

RSd1 = cat(1,resultados,resultadosV);

%%

% Importar datos de Radiosonda de NKX
load NKXRadiosonda.mat

% Matriz de resultados de Radiosonda.
resultadosR = zeros(15,1);
n= 1;
 while n < length(B)+1

    % Obtenemos los datos del archivo
    file = filename(B(n));
    direccion = strcat('C:\Users\Gabriela\Downloads\datoscal\',file);
    load(direccion);

    % Creacion de los vectores
    time = time_day(:,1);                     % tiempo

    %Se inicial las variables dependientes de la radiosonda
    n_clouds = NaN;
    LCL_srf = NaN;
    z_cloudbase= NaN;
    z_inv_base= NaN;
    z_inv_top= NaN;
    qT_BL= NaN;
    qT_jump= NaN;
    qT_3km= NaN;
    thetaL_BL= NaN;
    thetaL_jump= NaN;
    thetaL_3km= NaN;
    BLwnd_dir_avg= NaN;
    BLwnd_spd_avg= NaN;
    dq_decoupling= NaN;
    dtheta_decoupling= NaN;

    cont = 1 ;
    %Busqueda del dato de la radiosonda correspinduente al día
    while cont<1155
        
        if  day(NKXsoundingstable.date_i(cont,1))==day(time(23,1)) && month(NKXsoundingstable.date_i(cont,1))== month(time(23,1)) && year(NKXsoundingstable.date_i(cont,1))== year(time(23,1))
                
                n_clouds = NKXsoundingstable.n_clouds(cont,1);
                LCL_srf = NKXsoundingstable.LCL_srf(cont,1);
                z_cloudbase= NKXsoundingstable.z_cloudbase(cont,1);
                z_inv_base= NKXsoundingstable.z_inv_base(cont,1);
                z_inv_top= NKXsoundingstable.z_inv_top(cont,1);
                qT_BL= NKXsoundingstable.qT_BL(cont,1);
                qT_jump= NKXsoundingstable.qT_jump(cont,1);
                qT_3km= NKXsoundingstable.qT_3km(cont,1);
                thetaL_BL= NKXsoundingstable.thetaL_BL(cont,1);
                thetaL_jump= NKXsoundingstable.thetaL_jump(cont,1);
                thetaL_3km= NKXsoundingstable.thetaL_3km(cont,1);
                BLwnd_dir_avg= NKXsoundingstable.BLwnd_dir_avg(cont,1);
                BLwnd_spd_avg= NKXsoundingstable.BLwnd_spd_avg(cont,1);
                dq_decoupling= NKXsoundingstable.dq_decoupling(cont,1);
                dtheta_decoupling= NKXsoundingstable.dtheta_decoupling(cont,1);

        end

        cont = cont +1 ;

    end
 


    % Copia de los resultados en la matriz 
    resultadosR(1,n) = n_clouds;
    resultadosR(2,n) = LCL_srf;
    resultadosR(3,n) = z_cloudbase;
    resultadosR(4,n) = z_inv_base;
    resultadosR(5,n) = z_inv_top;
    resultadosR(6,n) = qT_BL;
    resultadosR(7,n) = qT_jump;
    resultadosR(8,n) = qT_3km;
    resultadosR(9,n) = thetaL_BL;
    resultadosR(10,n) = thetaL_jump;
    resultadosR(11,n) = thetaL_3km;
    resultadosR(12,n) = BLwnd_dir_avg;
    resultadosR(13,n) = BLwnd_spd_avg;
    resultadosR(14,n) = dq_decoupling;
    resultadosR(15,n) = dtheta_decoupling ;

    n = n+1;

 end


RSd = cat(1,RSd1,resultadosR);
