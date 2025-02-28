# -*- coding: utf-8 -*-

import pyEDM
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.feature_selection import mutual_info_regression

# Load dataset
BASE = pd.read_csv(r"C:\\Users\\Gabriela\\Downloads\\Proyecto\\BaseCCM20112024.csv", delimiter=';', parse_dates=['Tiempo'])

# Convert columns to one-dimensional arrays
Radiacion = BASE['Radiacion Solar']
Velocidad = BASE['Velocidad']

# Generate time index
t = np.arange(0, len(Velocidad), 1)

# Create DataFrame with unidimensional columns
data = pd.DataFrame({'Time': t,'Radiacion': Radiacion, 'Velocidad': Velocidad})

# Interpolate missing values in each column
data['Radiacion'] = data['Radiacion'].interpolate(method='linear')
data['Velocidad'] = data['Velocidad'].interpolate(method='linear')

#%%
from sklearn.feature_selection import mutual_info_regression

# Compute tau for GHI
lags_R, mi_R = mutual_information(data['Radiacion'].to_numpy(), max_lag=10)

# Plot mutual information for GHI
plt.plot(lags_R, mi_R, marker='o', label='Mutual Information (GHI)')
plt.xlabel(r'$\tau$')
plt.ylabel('Mutual Information')
plt.legend()
plt.title(r'Selection of $\tau$ based on Mutual Information (GHI)')
plt.show()

# Optimal tau
tau_optimo = np.argmax(mi_R) + 1
print(f"Optimal Tau for GHI: {tau_optimo}")

# Compute tau for Wind Speed
lags_Ve, mi_Ve = mutual_information(data['Velocidad'].to_numpy(), max_lag=10)

# Plot mutual information for Wind Speed
plt.plot(lags_R, mi_Ve, marker='o', label='Mutual Information (Wind Speed)')
plt.xlabel(r'$\tau$')
plt.ylabel('Mutual Information')
plt.legend()
plt.title(r'Selection of $\tau$ based on Mutual Information (Wind Speed)')
plt.show()

# Optimal tau
tau_optimo = np.argmax(mi_Ve) + 1
print(f"Optimal Tau for Wind Speed: {tau_optimo}")

#%%

# It is recommended to vary lib="1 2000" and pred="2001 5000" to evaluate different parts of the cycle.

optimal_emdeddingr = pyEDM.EmbedDimension(dataFrame=data, lib="1 2000", pred="2001 5000", columns='Radiacion', target='Radiacion')
print(optimal_emdeddingr)

optimal_emdeddingvel = pyEDM.EmbedDimension(dataFrame=data, lib="1 2000", pred="2001 5000", columns='Velocidad', target='Velocidad')
print(optimal_emdeddingvel)

import pandas as pd
import seaborn as sns

# Convert results to DataFrame 
df_radiacion = pd.DataFrame(optimal_emdeddingr)
df_velocidad = pd.DataFrame(optimal_emdeddingvel)

# Configure Seaborn style
sns.set(style="whitegrid")


fig, ax = plt.subplots(figsize=(8, 5))

# Plot rho vs E curves
sns.lineplot(data=df_radiacion, x="E", y="rho", marker="o", label="Radiacion", ax=ax)
sns.lineplot(data=df_velocidad, x="E", y="rho", marker="s", label="Velocidad", ax=ax, linestyle="dashed")

# Highlight optimal E values (maximum rho)
opt_E_rad = df_radiacion.loc[df_radiacion["rho"].idxmax(), "E"]
opt_rho_rad = df_radiacion["rho"].max()
opt_E_vel = df_velocidad.loc[df_velocidad["rho"].idxmax(), "E"]
opt_rho_vel = df_velocidad["rho"].max()

ax.axvline(opt_E_rad, linestyle="--", color="blue", alpha=0.5)
ax.axvline(opt_E_vel, linestyle="--", color="red", alpha=0.5)
ax.text(opt_E_rad, opt_rho_rad, f"Optimal E (GHI)", color="blue", fontsize=10, verticalalignment='bottom')
ax.text(opt_E_vel, opt_rho_vel, f"Optimal E (Wind Speed)", color="red", fontsize=10, verticalalignment='bottom')
ax.set_title("Determination of Embedding Dimension E")
ax.set_xlabel("Dimension (E)")
ax.set_ylabel(r'$\rho$')
ax.legend()
plt.show()

#%%
from pyEDM import PredictInterval, PredictNonlinear

PredictInterval(dataFrame=data, lib="1 1000", pred="301 600", columns="Radiacion", E=5)
PredictInterval(dataFrame=data, lib="1 1000", pred="301 600", columns="Velocidad", E=5)
PredictNonlinear(dataFrame=data, lib="1 1000", pred="1001 700", columns="Radiacion", E=5)
PredictNonlinear(dataFrame=data, lib="1 1000", pred="1001 700", columns="Velocidad", E=5)

#%%
from pyEDM import CCM

window_size = 8000  # Window size
overlap = 1000  # Overlap between windows
results = []

# Perform CCM analysis in sliding windows
for i in range(0, len(data) - window_size + 1, window_size - overlap):
    window_data = data.iloc[i:i + window_size]
    result = CCM(dataFrame=window_data, E=5, columns="Radiacion", target="Velocidad", libSizes="50 2500 250", sample=10, showPlot=True, tau=1)
    results.append(result)

#%%
import numpy as np
import matplotlib.pyplot as plt

lib_sizes = None
corr_radiacion_u = []  # Correlations for GHI -> Wind Speed
corr_u_radiacion = []  # Correlations for Wind Speed -> GHI

# Process results
for result in results:
    if lib_sizes is None:
        lib_sizes = result['LibSize']
    corr_radiacion_u.append(result['Radiacion:Velocidad'])
    corr_u_radiacion.append(result['Velocidad:Radiacion'])

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
plt.errorbar(lib_sizes, avg_radiacion_u, yerr=sem_radiacion_u, fmt='o', label="Wind Speed → GHI", color="red", capsize=4)
plt.errorbar(lib_sizes, avg_u_radiacion, yerr=sem_u_radiacion, fmt='o', label="GHI → Wind Speed", color="blue", capsize=4)
plt.xlabel("Library Size (L)", fontsize=30)
plt.ylabel("Correlation (ρ)", fontsize=30)
plt.title(r"CCM ($E=5$, $\tau=1$), Antofagasta", fontsize=30)
plt.show()
