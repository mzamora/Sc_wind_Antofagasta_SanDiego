# -*- coding: utf-8 -*-
"""
@author: G Pallauta
"""


import pandas as pd
import numpy as np
import pyEDM
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.feature_selection import mutual_info_regression
from joblib import Parallel, delayed


# Load dataset
data = pd.read_csv(r"C:\\Users\\Gabriela\\Downloads\\Proyecto\\Datos_Combinados.csv")

data['Tiempo'] = pd.to_datetime(data['Tiempo'], errors='coerce')

data = data.dropna(subset=['Tiempo'])
data.set_index('Tiempo', inplace=True)
data = data.tz_localize('UTC')
data['Fecha'] = data.index.normalize()

# Filter data after 14:00 UTC (sunrise in San Diego)
data = data[data.index.hour >= 14]
data = data[~data.index.duplicated(keep='first')]

#Interpolate missing data in columns 'VelocidadViento' and 'Radiacion'
cols_interpolar = ['VelocidadViento', 'Radiacion']
for col in cols_interpolar:
    if col in data.columns:
        data[col] = data[col].interpolate(method='cubic')

# Remove specific dates
fechas_a_eliminar = pd.to_datetime(['2017-09-30', '2017-08-31']).tz_localize('UTC')
data = data[~data['Fecha'].isin(fechas_a_eliminar)]

# Group data by date
vectores_por_dia = {
    fecha: grupo.resample('T').mean()
    for fecha, grupo in data.groupby('Fecha')
}



#%% Plotting GHI and Wind Speed for each day

for fecha, datos in vectores_por_dia.items():
 
    fig, ax1 = plt.subplots(figsize=(12, 6))
    
    # Plot GHI
    ax1.plot(datos.index, datos['Radiacion'], label='Radiación (W/m²)', color='red', linewidth=2)
    ax1.set_xlabel('Time (UTC)', fontsize=20)
    ax1.set_ylabel('GHI (W/m²)', fontsize=20, color='red')
    ax1.tick_params(axis='y', labelcolor='red')
    
    # Plot Wind Speed
    ax2 = ax1.twinx()
    ax2.plot(datos.index, datos['VelocidadViento'], label='Velocidad del Viento (m/s)', color='blue', linewidth=2)
    ax2.set_ylabel('Wind Speed (m/s)', fontsize=20, color='blue')
    ax2.tick_params(axis='y', labelcolor='blue')
    
 
    plt.title(f"GHI and Wind Speed ({fecha.strftime('%d-%m-%Y')}), La Jolla", fontsize=16)
    fig.tight_layout()
    plt.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    plt.show()


#%% Mutual Information Calculation

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.feature_selection import mutual_info_regression
from joblib import Parallel, delayed

def calcular_mutual_info(serie, max_lag=10, n_jobs=-1):
    
    serie = serie[~np.isnan(serie)]  

    def calcular_mi(lag):
        X = serie[:-lag].reshape(-1, 1)
        y = serie[lag:]
        min_len = min(len(X), len(y))
        return mutual_info_regression(X[:min_len], y[:min_len])[0]

    ami_vals = Parallel(n_jobs=n_jobs)(delayed(calcular_mi)(lag) for lag in range(1, max_lag))
    return ami_vals

def calcular_tau_ami(serie, max_lag=10):
    ami_vals = calcular_mutual_info(serie, max_lag)

    
    tau_optimo = next((i + 1 for i in range(1, len(ami_vals) - 1) if ami_vals[i] < ami_vals[i - 1] and ami_vals[i] < ami_vals[i + 1]), 1)

    # Function to find optimal tau
    plt.figure(figsize=(8, 4))
    plt.plot(range(1, max_lag), ami_vals, marker='o', linestyle='-')
    plt.xlabel('Retraso (τ)')
    plt.ylabel('Información Mutua')
    plt.title('Selección de τ usando Información Mutua')
    plt.legend()
    plt.show()

    return tau_optimo


print(f"Optimal tau (GHI): {calcular_tau_ami(data['Radiacion'].values)}")
print(f"Optimal tau (Wind Speed): {calcular_tau_ami(data['VelocidadViento'].values)}")

#%%

# Sliding Window for Analysis
window_size =599  
overlap = 0       

resultados_por_dia = {}  
resultados_E = []
resultados_E2 = []


for fecha, datos in vectores_por_dia.items():
    print(f"Procesando el día: {fecha}")
    
    if len(datos) < window_size:
        print(f"Día {fecha}: datos insuficientes para formar ventanas.")
        continue

    
    columnas_requeridas = ['VelocidadViento', 'Radiacion']
    try:
        datos_numericos = datos.loc[:, columnas_requeridas]
    except KeyError as e:
        print(f"Error: Columnas requeridas no están presentes en los datos del día {fecha}: {e}")
        continue

    resultados_del_dia = []

    step = window_size - overlap  
    for i in range(0, len(datos_numericos) - window_size + 1, step):
        print(f"Ventana de {i} a {i + window_size} ({len(datos_numericos)})")  
        
        ventana = datos_numericos.iloc[i:i + window_size]
        t = np.arange(0, len(ventana), 1)
        data = pd.DataFrame({'Time': t, 'Radiacion': ventana['Radiacion'].values, 'VelocidadViento': ventana['VelocidadViento'].values})
    
        #  CCM 
        try:
            result = pyEDM.CCM(
                dataFrame=data,
                E=2, 
                columns="Radiacion",  
                target="VelocidadViento", 
                libSizes="20 550 50",  
                sample=10,  
                showPlot=False,  
                tau=1  
            )
            resultados_del_dia.append(result)
            print(f"Resultado de CCM (Ventana {i}): {result}")
            
            optimal_embedding = pyEDM.EmbedDimension(dataFrame=data, lib="1 150", pred="151 400", columns='Radiacion', target='VelocidadViento', showPlot=False)
            resultados_E.append(optimal_embedding)  
            
            optimal_embedding2 = pyEDM.EmbedDimension(dataFrame=data, lib="1 150", pred="151 400", columns='VelocidadViento', target='Radiacion', showPlot=False)
            resultados_E2.append(optimal_embedding2)  
        except Exception as e:
            print(f"Error en CCM para ventana {i}: {e}")
    
    resultados_por_dia[fecha] = resultados_del_dia  


#%%

suma_rho_por_E = {i: 0 for i in range(2, 11)}
suma_rho_por_E2 = {i: 0 for i in range(2, 11)}

for df in resultados_E:
    for i in range(2, 11): #Parte desde E = 2
        suma_rho_por_E[i] += df[df['E'] == i]['rho'].values[0]
        
    
for df in resultados_E2:
    for i in range(2, 11): 
        suma_rho_por_E2[i] += df[df['E'] == i]['rho'].values[0]
        
max_E = max(suma_rho_por_E, key=suma_rho_por_E.get)
max_rho = suma_rho_por_E[max_E]/48

max_E2 = max(suma_rho_por_E2, key=suma_rho_por_E.get)
max_rho2 = suma_rho_por_E2[max_E]/48


E_values = list(range(2, 11))  # E [1,10]
rho_E_values = [suma_rho_por_E.get(E, 0) for E in E_values]  
rho_E2_values = [suma_rho_por_E2.get(E, 0) for E in E_values]  

# Crear el gráfico de barras
plt.figure(figsize=(10, 6))

# Plot de barras para resultados_E
plt.bar(E_values, rho_E_values, width=0.4, label='GHI-Wind Speed', align='center', alpha=0.7)

# Plot de barras para resultados_E2
plt.bar(E_values, rho_E2_values, width=0.4, label='Wind Speed-GHI', align='edge', alpha=0.7)
plt.xlabel("Dimension (E)")
plt.ylabel(r'$\rho$')
plt.title("Determination of Embedding Dimension E")
plt.xticks(E_values)
plt.legend()
plt.show()


#%%


lib_sizes = None
corr_radiacion_u = []  # Correlations for GHI -> Wind Speed
corr_u_radiacion = []  # Correlations for Wind Speed -> GHI

# Process results
for fecha, resultados in resultados_por_dia.items():
    for result in resultados:
        lib_sizes.append(result['LibSize'])
        corr_radiacion_u.append(result['Radiacion:VelocidadViento'])
        corr_u_radiacion.append(result['VelocidadViento:Radiacion'])

# Convert to numpy arrays
corr_radiacion_u = np.array(corr_radiacion_u)
corr_u_radiacion = np.array(corr_u_radiacion)


# Compute mean and standard error of the mean (SEM)
avg_radiacion_u = np.mean(corr_radiacion_u, axis=0)
sem_radiacion_u = np.std(corr_radiacion_u, axis=0) / np.sqrt(corr_radiacion_u.shape[0])
avg_u_radiacion = np.mean(corr_u_radiacion, axis=0)
sem_u_radiacion = np.std(corr_u_radiacion, axis=0) / np.sqrt(corr_u_radiacion.shape[0])


# Plot results
plt.figure(figsize=(16, 10))

# GHI -> Wind speed
plt.errorbar(lib_sizes[0], avg_radiacion_u, yerr=sem_radiacion_u, fmt='o', label="Wind Speed → GHI", color="red", capsize=4)
plt.plot(lib_sizes[0], avg_radiacion_u, color="darkred", linestyle='--')  

# Wind Speed -> GHI
plt.errorbar(lib_sizes[0], avg_u_radiacion, yerr=sem_u_radiacion, fmt='o', label="GHI → Wind Speed", color="blue", capsize=4)
plt.plot(lib_sizes[0], avg_u_radiacion, color="darkblue", linestyle='--')  


plt.legend(fontsize=25, loc='upper left')


plt.xlabel("Library Size (L)", fontsize=30)
plt.ylabel("Correlation (ρ)", fontsize=30)
plt.title(r"CCM ($E=2$, $\tau=1$), San Diego", fontsize=30)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(0.4, 0.72)  
plt.xlim(-20, max(lib_sizes[0])+20)  

plt.show()
