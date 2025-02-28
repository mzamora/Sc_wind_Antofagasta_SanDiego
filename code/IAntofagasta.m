
% clear all;
% clc;

load 'dAnt1.mat';

%Crear la matriz donde se guardan los resultados
resultados= zeros(64,length(dAnt1));

%%

%Se importan los datos desde la base de datos de Cerro Moreno.
load 'datacmorL1.mat';

time = datacmorL1.time;     % tiempo
temp = datacmorL1.temp;     % temperatura
hrel = datacmorL1.hrel;     % humedad relativa
rnet = datacmorL1.rnet;     % radiacion neta
velv = datacmorL1.velv;     % velocidad del viento
dirv = datacmorL1.dirv;     % direccion del viento
rglb =datacmorL1.rglb;      % radiacion global
pres =datacmorL1.pres;      % presion
tpot =datacmorL1.tpot;      % temperatura potencial
tdew =datacmorL1.tdew;      % temperatura de punto de rocio


t =datetime(time,'ConvertFrom','datenum',TimeZone='UTC');


%% Se importa los datos del Nefobasimetro

load 'Nubesbajas.mat'

%%

% Datos de la radiosonda.

load 'ALL.mat';


date = ALL_SOUNDINGS.date_i;         
n_clouds1 = ALL_SOUNDINGS.n_clouds;
LCL_srf1= ALL_SOUNDINGS.LCL_srf;
z_cloudbase1= ALL_SOUNDINGS.z_cloudbase;
z_inv_base1= ALL_SOUNDINGS.z_inv_base;
z_inv_top1= ALL_SOUNDINGS.z_inv_top;
qT_BL1= ALL_SOUNDINGS.qT_BL;
qT_jump1= ALL_SOUNDINGS.qT_jump;
qT_3km1= ALL_SOUNDINGS.qT_3km;
thetaL_BL1= ALL_SOUNDINGS.thetaL_BL;
thetaL_jump1= ALL_SOUNDINGS.thetaL_jump;
thetaL_3km1= ALL_SOUNDINGS.thetaL_3km;
BLwnd_dir_avg1= ALL_SOUNDINGS.BLwnd_dir_avg;
BLwnd_spd_avg1= ALL_SOUNDINGS.BLwnd_spd_avg;
dq_decoupling1= ALL_SOUNDINGS.dq_decoupling;
dtheta_decoupling1= ALL_SOUNDINGS.dtheta_decoupling;


date.TimeZone='UTC';

%%

%Contador
nummat= 1;

while nummat <= length(dAnt1)

    % Seleccion del dia buscado
    % inicio
    d1 = dAnt1(nummat);
    d1.TimeZone = 'UTC';

    % final
    d2 = d1+ day(1);
    dia = d1 + minutes(60*12);
    
    % Contadores de inicio y final del dia buscado
    inicio = 1;
    final = 1;
    
    i = 1  ;% El contador comienza en 1, ya que el orden de los dias buscados puede no ser cronologico.
            % Si la lista es cronologica se puede igualar i = final del
            % recorrido e inicial el contador antes de recorrer los dias.
    
    while i<=length(time)
        if (t(i)==d1) % Coincidencia entre los dias buscado y el contador.
            inicio=i;
        end
        if (t(i)==d2)
            final=i;
            break;
        end
        i=i+1;
    end

    % Se inicializan las variables nuevas, por dia.
    timet= time(inicio:final);   % tiempo
    rnett= rnet(inicio:final);   % radiacion neta        
    tempt= temp(inicio:final);   % temperatura
    velvt= velv(inicio:final);   % velocidad
    dirvt= dirv(inicio:final);   % direccion           
    rglbt= rglb(inicio:final);   % radiacion global
    prest= pres(inicio:final);   % presion
    tpott= tpot(inicio:final);   % temperatura potencial     
    tdewt= tdew(inicio:final);   % temperatura de punto de rocio
    hrelt= hrel(inicio:final);   % humedad relativa

    ti = datetime(timet,'ConvertFrom','datenum');

    
    % Cerro Moreno -23.466508498153548, -70.56416677506581
    Location.latitude = -23.466508498153548;
    Location.longitude = -70.56416677506581;
    Location.altitude = 115;
    UTCoffset = 0; % Porque los datos estan en UTC
    
    win_length = 20; % Consideramos 20 minutos de intervalo
    sample_interval = 5;
    
    
    % Demonstrate algorithm for unequally spaced data
    [CD] = pvl_detect_clear_times(rglbt, timet, UTCoffset, Location, win_length, sample_interval);
    
    sti =datetime(timet,'ConvertFrom','datenum');
    [pks,locs] = findpeaks(rglbt,sti);
    
    [peaksghi,tiemppeak] = findpeaks(rglbt,ti);
    radiacion = rglbt(CD);

    tiempocd = datenum(sti(CD));
    DN = tiempocd.';
    Time = pvl_maketimestruct(DN, 0);

    
    %vector de todos los indices de turbidez que se prueban 0.001
    vectorCS = 0:0.001:8;

    %Matriz con todos los cielos despejados/ 289 son los datos que hay cada dia
    %estudiado
    VectorCD = zeros(length(radiacion),length(vectorCS)-1);
    
    %Guarda la comparacion
    VectorComparacion = 0:length(vectorCS)-1;

    %contador
    indiceCS = 1;
    
    %Calcular los peaks de rglbt
    while indiceCS <= length(vectorCS)        
    
            [ClearSkyGHI]= pvl_clearsky_ineichen(Time, Location,vectorCS(indiceCS));
    
            VectorCD(:,indiceCS) = ClearSkyGHI;
    
            % Comparacion con la entre el cielo despejado y la radiacion
            % medida.
     
            VectorComparacion(indiceCS) = abs(nanmean(ClearSkyGHI./radiacion-1));
    
            indiceCS = indiceCS + 1;
    end
    
    
    DN = timet.';
    Time = pvl_maketimestruct(DN, 0);
    
    %Buscamos el indice del vector decomparacion menor 
    k = find(min(VectorComparacion.')==VectorComparacion.');

    % Se encuentra la turbulencia con la comparacion mas cercana a 0 ( igualdad con la radiacion)
    turb =vectorCS(k);
    
    %Crea el cielo despejado con la turbulencia 
    [ClearSkyGHI]= pvl_clearsky_ineichen(Time, Location,turb);
    
    % Calculo de la direccion del viento en u y v.
    uA=-velvt.*sin(deg2rad(dirvt));
    vA=-velvt.*cos(deg2rad(dirvt)); % Velocidad predominante en este caso
    
    V = sqrt(uA.^2 + vA.^2); % magnitud del viento

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

    Vamanecer = V(f);
    
    %contador 
    inicioviento = 1 ;
    maginicioviento = 1;
    velvinicioviento =1 ;

    %contador
    contadorv = f; %parte desde el amanecer
    
    % Suavizado del viento en direccion v (predominante en Antofagasta)
    Vientosuavizado = smooth( timet,V,0.1,'rloess');
    
    while contadorv <=(length(V)-5)
            %Criterio para el tiempo inicial es que tenga un aumento por
            %los siguientes 20 minutos. Despues de las 07:00 UTC -> 04:00 Hora Local aprox. 
           
            if  Vientosuavizado(contadorv)<Vientosuavizado(contadorv+1) & Vientosuavizado(contadorv+1)<Vientosuavizado(contadorv+2) & Vientosuavizado(contadorv+2)<Vientosuavizado(contadorv+3) & Vientosuavizado(contadorv+3)<Vientosuavizado(contadorv+4) 
                  inicioviento = ti(contadorv);
                  maginicioviento= V(contadorv);
                  velvinicioviento = vA(contadorv);
                  
                  break;
            end
    
    
            contadorv=contadorv+1;
    end
    
    if contadorv>(length(V)-5)
        inicioviento = NaN;
        maginicioviento = NaN;
        veluinicioviento = NaN;

    end

    %Calculo de variables
    promvelv = nanmean(vA);  
    promvelu = nanmean(uA);
    prommag = nanmean(V);
 

    maxmagnitudviento = nanmax(V);
    minmagnitudviento = nanmin(V);

    
    %Indice de cielo despejado (k)
    indice = rglbt./ClearSkyGHI;
    
    
    %Iniciamos las variables
    tfinalruptura = NaN; %tiempo de ruptura final
    rfr = NaN; %radiación final
    Vin = NaN;
    
   
    
    tinicialruptura = NaN; %tiempo de ruptura final
    rir= NaN; %radiación inicial
    kruptura = NaN; %indice inicial
    Vfin = NaN;
    
    %contador
    p = f;  % La inicio del proceso de ruptura es posterior al amanecer
    
    while p <=(length(indice)) 
    
         if isnan(indice(p))==0
                  
                    %Criterio para el tiempo inicial es que sea un peak, que su
                    %indice sea menor a 1.2 y que sean anterior al primer
                    %cielo despejado
                    if indice(p)< 1.2 & isnan(find( locs==ti(p)))==0 & datenum(ti(p))<=min(datenum(sti(CD)))
                      
                        tinicialruptura = ti(p);
                        rir = rglbt(p);
                        kruptura = indice(p);
   
                        break;
                    end
         end

         if p == (length(indice)) 
            p = NaN;
            break;
         end
         p=p+1;
    end
    
    
    %contador
    q = p; %parte en p
    
    
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


        
    % Indices dependientes de la ruptura
    
    maginicioruptura = NaN ; 
    magfirup = NaN; 
        
    promagruptura =  NaN;
    maxmagruptura = NaN;
    minmagruptura = NaN;

    %radiacion
    rglbprom = NaN;
    rglbmax = NaN;
    rglbmin = NaN;
    stdrglb = NaN;
    
    %Temperatura
    tempromrup = NaN;
    tempminrup = NaN;
    tempmaxrup = NaN;
    stdtemp = NaN;
    
    %Humedad Relativa
    hrelpromrup = NaN;
    hrelmaxrup = NaN;
    hrelminrup = NaN;
    stdhrel = NaN;
    
    % Presion
    prespromrup = NaN;
    presmaxrup = NaN;
    presminrup = NaN;
    stdpres = NaN;

    %Temperatura potencial
    tpotpromrup = NaN;
    tpotmaxrup = NaN;
    tpotminrup = NaN;
    stdtpot = NaN;
    
    % Temperatura de punto de rocio
    tdewpromrup = NaN;
    tdewminrup = NaN;  
    tdewmaxrup = NaN;
    stdtdew = NaN;

    % Si existe el preoceso de ruptura, inicio y final se calculan las
    % variables dentro del proceso
    if isnan(p)==0 && isnan(q)==0

        vectorruptura = ti(p:q);
        Magrup = V(p:q);

        velvinicioruptura = vA(p);
        velvfirup = vA(q);
            
        % Indices dependientes de la ruptura
        
        maginicioruptura = V(p) ; 
        magfirup = V(q); 
            
        promagruptura =  nanmean(Magrup);
        maxmagruptura = nanmax(Magrup);
        minmagruptura = nanmin(Magrup);

        % Calculo de las variables atmosfericas
        %Radiacion
        vectorrglbtruptura = rglbt(p:q);
    
        rglbprom = nanmean(vectorrglbtruptura);
        rglbmax = nanmax(vectorrglbtruptura);
        rglbmin = nanmin(vectorrglbtruptura);
        stdrglb = std(vectorrglbtruptura,'omitnan');
    
        %Temperatura
        vectortemruptura = tempt(p:q);
    
        tempromrup = nanmean(vectortemruptura);
        tempminrup = nanmin(vectortemruptura);
        tempmaxrup = nanmax(vectortemruptura);
        stdtemp = std(vectortemruptura,'omitnan');
    
        %Humedad Relativa
        vectorhrelruptura = hrelt(p:q);
    
        hrelpromrup = nanmean(vectorhrelruptura);
        hrelmaxrup = nanmax(vectorhrelruptura);
        hrelminrup = nanmin(vectorhrelruptura);
        stdhrel = std(vectorhrelruptura,'omitnan');
    
        % Presion
        vectorpresruptura = prest(p:q);
    
        prespromrup = nanmean(vectorpresruptura);
        presmaxrup = nanmax(vectorpresruptura);
        presminrup = nanmin(vectorpresruptura);
        stdpres = std(vectorpresruptura,'omitnan');
    
        vectortpottruptura = tpott(p:q);
    
        %Temperatura potencial
        tpotpromrup = nanmean(vectortpottruptura);
        tpotmaxrup = nanmax(vectortpottruptura);
        tpotminrup = nanmin(vectortpottruptura);
        stdtpot = std(vectortpottruptura,'omitnan');
    
        % Temperatura de punto de rocio
        vectortdewruptura = tdewt(p:q);
    
        tdewpromrup = nanmean(vectortdewruptura);
        tdewminrup = nanmin(vectortdewruptura);  
        tdewmaxrup = nanmax(vectortdewruptura);
        stdtdew = std(vectortdewruptura,'omitnan');

    end


    promkrup = NaN;
    Integralindice = NaN;

    
    if isempty(indice)==0 && isnan(p)==0 && isnan(q)==0
        krango = indice(p:q);
        promkrup = nanmean(krango);
        Integralindice = trapz(krango);
    
    end
    
    Vamanecer = V(f);
    velvamanecer = vA(f);
    
    deltaruptura = NaN;

    if isnan(datenum(tfinalruptura))==0 & isnan(datenum(tinicialruptura))==0
        deltaruptura = tfinalruptura - tinicialruptura;
    end   

    % Calculo de  datos de radiosonda.

    % S inician las variables
    nnub = NaN;
    altbase = NaN;
    altcapinv = NaN;
    alttopcapainv = NaN;
    lcl = NaN;
    qTBL = NaN ;
    qTj = NaN ;
    qT3 = NaN ;
    thetaBL = NaN ;
    thetaj = NaN ;
    theta3 = NaN ;
    wnd_dir =NaN;
    wnd_spd = NaN;
    
    %Ingresar el dia buscado

    primercontador = 1;

    while primercontador <= length(date)

        if dia==date(primercontador)


            n_clouds = n_clouds1(primercontador);
            LCL_srf = LCL_srf1(primercontador);
            z_cloudbase = z_cloudbase1(primercontador);
            z_inv_base = z_inv_base1(primercontador);
            z_inv_top = z_inv_top1(primercontador);
            qT_BL =  qT_BL1(primercontador);
            qT_jump = qT_jump1(primercontador);
            qT_3km = qT_3km1(primercontador);
            thetaL_BL = thetaL_BL1(primercontador);
            thetaL_jump = thetaL_jump1(primercontador);
            thetaL_3km =  thetaL_3km1(primercontador);
            BLwnd_dir_avg =  BLwnd_dir_avg1(primercontador);
            BLwnd_spd_avg =  BLwnd_spd_avg1(primercontador);
            dq_decoupling = dq_decoupling1(primercontador);
            dtheta_decoupling = dtheta_decoupling1(primercontador);
            
            break;
        end

        primercontador = primercontador +1;
    end

    
    %Iniciamos los indices del Nefo por defecto como NaN
    promedionb = NaN;
    mediananb = NaN;
    desvnb = NaN;
    maxnb = NaN;
    minnb = NaN;
    
    
    % Si existe un tiempo de inicio de ruptura
    if isnan(datenum(tinicialruptura))==0 % Si es diferente de NaN
        

        %Tiempo del Nefo a comparar
        tnef = datetime(Nubesbajas(1,:),'ConvertFrom','datenum');
        contadornefo = 1 ;
        while contadornefo <= length(Nubesbajas(1,:))
           
            if ((tnef(1,contadornefo)-tinicialruptura)<minutes(5) && minutes(-5)<(tnef(1,contadornefo)-tinicialruptura))
                    
                    
                        minnb = nanmean(Nubesbajas(6,(contadornefo-12):contadornefo));   
                        maxnb = nanmean(Nubesbajas(5,(contadornefo-12):contadornefo));            
                        desvnb =  nanmean(Nubesbajas(4,(contadornefo-12):contadornefo));
                        mediananb =  nanmean(Nubesbajas(3,(contadornefo-12):contadornefo));
                        promedionb =  nanmean(Nubesbajas(2,(contadornefo-12):contadornefo));   
                    
            end


    
            contadornefo = contadornefo+1;
        end

  end

    
    %Copiar en la matriz de resultados


    % Cambiando el UTC por hora local con hour(time(f))*60 +
    % minute(time(f)); restandole 3 horas
    resultados(1,nummat)= (hour(tinicialruptura)-3)*60 + minute(tinicialruptura);
    resultados(2,nummat)= (hour(tfinalruptura)-3)*60 + minute(tfinalruptura);

    resultados(3,nummat)=NaN;

    if isnan(deltaruptura)==0
        resultados(3,nummat)= minutes(deltaruptura);
    end

    resultados(4,nummat)= NaN;

    if isempty(turb)==0
       resultados(4,nummat)= turb;
    end

    resultados(5,nummat)=kruptura;
    resultados(6,nummat)=(hour(tinicial)-3)*60 + minute(tinicial);
    resultados(7,nummat)=rglbprom;
    resultados(8,nummat)=rglbmax;
    resultados(9,nummat)=rglbmin;
    resultados(10,nummat)=stdrglb;
    resultados(11,nummat)=promkrup;


    resultados(12,nummat)=tempromrup;
    resultados(13,nummat)=tempminrup;
    resultados(14,nummat)=tempmaxrup;
    resultados(15,nummat)=stdtemp;
    resultados(16,nummat)=hrelpromrup;
    resultados(17,nummat)=hrelmaxrup;
    resultados(18,nummat)=hrelminrup;
    resultados(19,nummat)=stdhrel;
    resultados(20,nummat)=prespromrup;
    resultados(21,nummat)=presmaxrup;
    resultados(22,nummat)=presminrup;
    resultados(23,nummat)=stdpres;
    resultados(24,nummat)=tpotpromrup;
    resultados(25,nummat)=tpotmaxrup;
    resultados(26,nummat)=tpotminrup;
    resultados(27,nummat)=stdtpot;
    resultados(28,nummat)=tdewpromrup;
    resultados(29,nummat)=tdewminrup;
    resultados(30,nummat)= tdewmaxrup;
    resultados(31,nummat)=stdtdew;

 

   
    resultados(32,nummat) = n_clouds;
    resultados(33,nummat) = LCL_srf;
    resultados(34,nummat) = z_cloudbase;
    resultados(35,nummat) = z_inv_base;
    resultados(36,nummat) = z_inv_top;
    resultados(37,nummat) = qT_BL;
    resultados(38,nummat) = qT_jump;
    resultados(39,nummat) = qT_3km;
    resultados(40,nummat) = thetaL_BL;
    resultados(41,nummat) = thetaL_jump;
    resultados(42,nummat) = thetaL_3km;
    resultados(43,nummat) = BLwnd_dir_avg;
    resultados(44,nummat) = BLwnd_spd_avg;
    resultados(45,nummat) = dq_decoupling;
    resultados(46,nummat) = dtheta_decoupling ;


    resultados(47,nummat) = promedionb;
    resultados(48,nummat) = mediananb;
    resultados(49,nummat) = desvnb;
    resultados(50,nummat) = maxnb;
    resultados(51,nummat) = minnb;

    resultados(52,nummat) = promvelv;
    resultados(53,nummat) = promvelu;
    resultados(54,nummat) =prommag;
    resultados(55,nummat) =maxmagnitudviento;
    resultados(56,nummat) =minmagnitudviento;
    resultados(57,nummat) = (hour(inicioviento)-3)*60 + minute(inicioviento);
    resultados(58,nummat) = maginicioviento;
    resultados(59,nummat) = Vamanecer;
    resultados(60,nummat) = promagruptura;
    resultados(61,nummat) = maxmagruptura;
    resultados(62,nummat) = minmagruptura;
    resultados(63,nummat) = maginicioruptura;
    resultados(64,nummat) = magfirup;


    diatitulo = " (" +num2str(day(timet(1,1)))+"/" + num2str(month(timet(1,1)))+ "/" + num2str(year(timet(1,1))) + ")";
    diatitulo1 = num2str(day(timet(1,1)))+"-" + num2str(month(timet(1,1)))+ "-" + num2str(year(timet(1,1))) + ".png";


    if (isnan(datenum(tfinalruptura))==0 && isnan(datenum(tinicialruptura))==0)

        tiemposruptura = [tinicialruptura,tfinalruptura];
        radiacionrupt = [rir, rfr];


            figure('Visible', 'off');
            hold on
            title(strcat(strcat('GHI y Tiempos de fragmentación ',diatitulo), ', Antofagasta'));
            plot(ti,rglbt,'-r');
            plot(tiemposruptura,radiacionrupt,'*b');
            xlabel('Time (UTC)');
            ylabel('GHI (W/m^2)');
            legend('GHI','Tiempos de inicio'+ string(newline) + ' y término','Location','northwest');
            hold off

            imagen = gcf;
            nombredoc = strcat('Tiempos_',diatitulo1);
            exportgraphics(imagen,nombredoc);

    end

    nombredoc = strcat(num2str(nummat),'_imgAnt.png');

    if isempty(turb)==0

        figure('Visible', 'off');
        hold on
        plot(ti, rglbt);
        plot(ti,ClearSkyGHI);
        hold off
        title(strcat(strcat('Irradianza global ',diatitulo),', Antofagasta'));
        xlabel('Hora (UTC)');
        ylabel('GHI (W/m^2)');
        legend('GHI','GHI_{d}','Location','northwest');

        nombredoc = strcat('img_',diatitulo1);
        imagen = gcf;
        exportgraphics(imagen,nombredoc);

    end


    figure('Visible', 'off');
    hold all
    plot(sti, rglbt, 'b-');
    plot(sti(CD), rglbt(CD), 'r.');
    xlabel('Hora (UTC)');
    ylabel('GHI (W/m^2)');
    title(strcat(strcat('Momentos de Cielo Despejado ',diatitulo), ', Antofagasta'));
    legend('GHI','Momentos CD','Location','northwest');

    imagen = gcf;
    nombredoc = strcat('Desp_Ant_', diatitulo1);
    exportgraphics(imagen,nombredoc);



    nummat= nummat + 1;
end