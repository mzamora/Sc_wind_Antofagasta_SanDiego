# Based on original Matlab code by M. Zamora and X. Zhong @UCSD 
# Translated to python using Chat GPT, minor adaptions by M. Zamora @UCh

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from calendar import monthrange
from matplotlib import pyplot as plt  # Uncomment if plotting is needed

# Initialize variables
ind = 0
Cp = 1004
Rd = 287.05
Lv = 2.5e6

# other functions

def get_sounding_var(time, time_end, path_name, list_of_vars):
    """
    Extracts sounding variables from a file within a given time range.

    Parameters:
        time (datetime): Start time for extracting data.
        time_end (datetime): End time for extracting data.
        path_name (str): Path to the file containing the sounding data.
        list_of_vars (list of str): List of variable names to extract.

    Returns:
        dict: A dictionary where keys are variable names, and values are arrays of data.
    """

    var_list = ['PRES', 'HGHT', 'TEMP', 'DWPT', 'RELH', 'MIXR', 'WDIR', 'WSPD', 'THTA', 'THTE', 'THTV']

    # Read file using pandas
    with open(path_name, 'r') as file:
        raw_lines = file.readlines()

    # Convert lines to list of lists (split by commas)
    raw_data = [line.strip().split(',') for line in raw_lines]

    # Find the indices corresponding to the start and end times
    time_str = time.strftime('%Y%m%d_%H')
    time_end_str = (time + timedelta(hours=12)).strftime('%Y%m%d_%H')
    index_start = next((i for i, line in enumerate(raw_lines) if time_str in line), None) + 2

    if abs((time - time_end).total_seconds()) < 60:
        index_end = len(raw_data)
    else:
        index_end = next((i for i, line in enumerate(raw_data) if time_end_str in line), None)
        if index_end is None:
            time_next_day = (time + timedelta(hours=24)).strftime('%Y%m%d_%H')
            try:
                index_end = next((i for i, line in enumerate(raw_data) if time_next_day in line), None) - 1
            except Exception:
                index_end=-1

    # Extract and process variables
    result = {}
    for var in list_of_vars:
        if var not in var_list:
            print(f"Wrong input variable name: {var}")
            #return result

        var_index = var_list.index(var)
        var_data = [line[var_index] for line in raw_data[index_start:index_end] if line[var_index] != '']

        # Convert to numeric data
        var_data = np.array([float(value) for value in var_data])

        # Apply conversions if needed
        if var in ['TEMP', 'DWPT']:
            var_data += 273.15
        elif var == 'WSPD':
            var_data *= 0.514  # Convert knots to m/s
        result[var] = var_data

    # Remove duplicates and clean up
    if 'TEMP' in result:
        tp = result['TEMP']
        duplicate_indices = [i for i in range(len(tp) - 1) if tp[i] == tp[i + 1]]
        for key in result:
            result[key] = np.delete(result[key], duplicate_indices)

    # Remove rows with NaN values
    nan_mask = np.isnan(list(result.values())[0])
    for key in result:
        result[key] = result[key][~nan_mask]

    return result


def cal_lcl_bolton(TSK, PSFC, HGT, Tdew):
    """
    Calculate LCL temperature, pressure, and height based on Bolton's formula.

    Parameters:
        TSK (numpy array): Surface temperature in Kelvin.
        PSFC (numpy array): Surface pressure in Pa.
        QVG (numpy array): Mixing ratio in kg/kg.
        QSG (numpy array or None): Saturation mixing ratio in kg/kg. If None, it's calculated.
        HGT (numpy array): Surface height in meters.
        Tdew (numpy array): Dew point temperature in Kelvin.

    Returns:
        T_LCL2 (numpy array): LCL temperature in Kelvin.
        P_LCL2 (numpy array): LCL pressure in Pa.
        H_LCL2 (numpy array): LCL height in meters.
    """

    # Constants
    Rd = 287.04  # Specific gas constant for dry air (J/(kg-K))
    g = 9.80665  # Gravitational constant (m/s^2)

    # Initialize output arrays
    T_LCL2 = np.zeros_like(TSK)
    P_LCL2 = np.zeros_like(PSFC)
    H_LCL2 = np.zeros_like(HGT)

    # Loop through the grid
    for i in range(len(TSK)):
        # Calculate LCL temperature using Bolton's formula
        T_LCL2[i] = 1 / (1 / (Tdew[i] - 56) + np.log(TSK[i] / Tdew[i]) / 800) + 56
        # Calculate LCL pressure using the Poisson equation
        P_LCL2[i] = PSFC[i] * (T_LCL2[i] / TSK[i]) ** 3.5
        # Calculate LCL height using the hypsometric equation
        H_LCL2[i] = Rd * (TSK[i] + T_LCL2[i]) / 2 / g * np.log(PSFC[i] / P_LCL2[i]) + HGT[i]
    return T_LCL2, P_LCL2, H_LCL2




def tmp_inversion_strength_cal(T, hght):
    """
    Calculate temperature inversion strength.

    Parameters:
        T (numpy array): Temperature profile (K), 1D array.
        hght (numpy array): Height profile (m), 1D array.

    Returns:
        DT_max (float): Maximum temperature inversion strength (K).
        numofInv (int): Number of inversions.
        hght_top (float): Height of inversion top (km).
        hght_base (float): Height of inversion base (km).
        eta_top (int): Index of inversion top.
        eta_base (int): Index of inversion base.
    """

    hght=hght/1000 #m to km

    # Constants
    H_Crit = 3  # Maximum height to consider (km)
    DT_Crit=1 # minimum inversion jump (K)

    # Initialize outputs
    DT_max = 0
    numofInv = 0
    hght_top = np.nan
    hght_base = np.nan
    eta_top = np.nan
    eta_base = np.nan

    # Restrict height within H_Crit km
    mask = hght < H_Crit
    T = T[mask]
    hght = hght[mask]

    # Calculate temperature gradient
    dT = np.sign(np.diff(T))

    # Find indices of positive temperature gradients
    indices = np.where(dT > 0)[0]

    if len(indices) == 0:
        return DT_max, numofInv, hght_top, hght_base, eta_top, eta_base

    # Identify consecutive segments of positive gradients
    dI = np.diff(indices)
    segments = np.split(indices, np.where(dI > 1)[0] + 1)

    # Initialize inversion characteristics
    DT = []
    DT_DZ = []
    eta_top_temp = []
    eta_base_temp = []

    for segment in segments:
        base = segment[0]
        top = segment[-1] + 1  # Include the last point of the segment
        DT.append(T[top] - T[base])
        DT_DZ.append((T[top] - T[base]) / (hght[top] - hght[base]))
        eta_base_temp.append(base)
        eta_top_temp.append(top)

    # Determine the strongest inversion
    DT_max = max(DT)
    numofInv = len(segments)

    index_max = DT.index(DT_max)
    eta_top = eta_top_temp[index_max]
    eta_base = eta_base_temp[index_max]
    hght_top = hght[eta_top]
    hght_base = hght[eta_base]

    # Check inversion strength
    if DT_max < DT_Crit: # original was 3 K
        return 0, 0, np.nan, np.nan, np.nan, np.nan

    # Adjust inversion top height if temperature gradient < 5
    try:
        dTdz = np.diff(T) / np.diff(hght)
        first_low_gradient = np.where(dTdz[eta_base:eta_top] < 5)[0][0]
        adjusted_top = eta_base + first_low_gradient
        adjusted_hght_top = hght[adjusted_top]

        if adjusted_hght_top < hght_top:
            hght_top = adjusted_hght_top
            eta_top = adjusted_top
            DT_max = T[eta_top] - T[eta_base]
            if DT_max < DT_Crit:
                return 0, 0, np.nan, np.nan, np.nan, np.nan
    except IndexError:
        pass

    return DT_max, numofInv, hght_top, hght_base, eta_top, eta_base

#%% Main code

# Initialize lists to store results
date_i = []
n_clouds = []
LCL_srf = []
z_cloudbase = []
dLCL_decoupling = []
LCL_CB = []
dtv1 = []
dtv2 = []
fog = []
z_inv_base = []
z_inv_top = []
sfc_inversion = []
top_inversion = []
qT_BL = []
qT_3km = []
qT_jump = []
thetaL_BL = []
thetaL_3km = []
thetaL_jump = []
q_bot = []
q_up = []
theta_bot = []
theta_up = []
dq_decoupling = []
dtheta_decoupling = []
BLwnd_spd_avg = []
BLwnd_dir_avg = []

# Define the years and months to process
for yy in range(2016, 2018):
    for mm in range(1,13):
        dmax=monthrange(yy,mm)[1]
        for dd in range(1,dmax+1):
            for hh in [0, 12]:
                # Construct the path name for the data file
                pathName = f'NKX/csv/72293_{yy}_{mm:02d}_0100_{dmax:02d}12.csv'
                # List of parameters available
                VarList = ['PRES', 'HGHT', 'TEMP', 'DWPT', 'RELH', 'MIXR', 'WDIR', 'WSPD', 'THTA', 'THTE', 'THTV']
                date0 = datetime(yy, mm, dd, hh, 0, 0)
                
                # Get sounding variables (You need to implement this function)
                try:
                    daytable=get_sounding_var(date0, date0 + timedelta(hours=12), pathName, VarList)
                except Exception:
                    #print('Day missing')
                    continue
                p=daytable[VarList[0]]
                z=daytable[VarList[1]]
                tp=daytable[VarList[2]]
                td=daytable[VarList[3]]
                RH=daytable[VarList[4]]
                w=daytable[VarList[5]]
                winddir=daytable[VarList[6]]
                windspeed=daytable[VarList[7]]
                theta=daytable[VarList[8]]
                theta_e=daytable[VarList[9]]
                theta_v=daytable[VarList[10]]
                                
                if z.size == 0:
                    continue
    
                date_i.append(date0)
    
                # Filter sounding data
                zup = 3000  # Up to 3km
                # Wipe out NaN values and data above desired height or weird values
                f = np.isnan(w) | (z > zup) | (w > 50)
                p = p[~f]
                td = td[~f]
                tp = tp[~f]
                w = w[~f]
                theta = theta[~f]
                theta_e = theta_e[~f]
                theta_v = theta_v[~f]
                windspeed = windspeed[~f]
                winddir = winddir[~f]
                z = z[~f]
                RH = RH[~f]
    
                #fig,ax=plt.subplots(1,3,figsize=[9,4])
                #ax[0].plot(theta,z)
                #ax[1].plot(w,z)
                #ax[2].plot(RH,z)
    
                # Find LCL
                try:
                    # You need to implement Cal_LCL_Bolton function
                    _, _, LCL_bolton = cal_lcl_bolton(tp, p, z, td)
                    LCL_srf.append(LCL_bolton[0])
                except Exception:
                    LCL_srf.append(np.nan)
    
                # Cloud base and number of cloud segments from RH profile (last cloud)
                #ax[2].plot([95,95],[z[0],z[-1]],'r:')
                #ax[2].plot([0,100],[LCL_bolton[0],LCL_bolton[0]],'r--')
                cloudpoints = np.where(RH > 95)[0]
                if len(cloudpoints)>0:
                    n_clouds.append(1)
                    zb_index = cloudpoints[0]
                    z_cloudbase.append(z[cloudpoints[0]])
                    for ic in range(len(cloudpoints) - 1):
                        if cloudpoints[ic + 1] != (cloudpoints[ic] + 1):
                            n_clouds[-1] += 1
                            z_cloudbase[-1] = z[cloudpoints[ic + 1]]
                            zb_index = cloudpoints[ic + 1]
                    dLCL_decoupling.append(z_cloudbase[-1] - LCL_bolton[0])
                    try:
                        LCLCB = LCL_bolton[np.where(LCL_bolton < z_cloudbase[-1])[0][-1]]
                        LCL_CB.append(LCLCB)
                    except Exception:
                        LCL_CB.append(np.nan)
                    dtv1.append(theta_v[zb_index] - theta_v[0])
                    dtv2.append(theta_v[zb_index] - theta_v[1])
                else:
                    n_clouds.append(np.nan)
                    z_cloudbase.append(np.nan)
                    dLCL_decoupling.append(np.nan)
                    LCL_CB.append(np.nan)
                    zb_index = np.nan
                    dtv1.append(np.nan)
                    dtv2.append(np.nan)
    
                # Detect if there's fog present
                try:
                    if RH[0] > 95:
                        fog.append(1)
                    else:
                        fog.append(0)
                except Exception:
                    fog.append(np.nan)
    
                # Find inversion
                try:
                    DT_max, numofInv, hght_top, hght_base, _, _ = tmp_inversion_strength_cal(tp, z)
                    _, _, inv_w_top, inv_w_base, _, _ = tmp_inversion_strength_cal(-w, z)
                    assert not np.isnan(hght_base)
                    z_inv_base.append(hght_base * 1000)  # Inversion base height (IBH)
                    z_inv_top.append(min(hght_top, inv_w_top) * 1000)  # Inversion top height
                    zi_index1 = np.where(z == round(z_inv_base[-1]))[0][0]  # Find index for IBH
                    if zi_index1 == 0:
                        #print('Inversion at sfc')
                        sfc_inversion.append(1)
                        top_inversion.append(0)
                    elif z_inv_top[-1] == z[-1]:
                        #print('Inversion at top')
                        top_inversion.append(1)
                        sfc_inversion.append(1)
                    else:
                        sfc_inversion.append(0)
                        top_inversion.append(0)
                    zi_index2 = np.where(z == round(z_inv_top[-1]))[0][0]
                except Exception:
                    z_inv_base.append(np.nan)
                    z_inv_top.append(np.nan)
                    sfc_inversion.append(np.nan)
                    top_inversion.append(np.nan)
                    zi_index1 = np.nan
                    zi_index2 = np.nan
    
                # Idealized well-mixed profile
                try:
                    qT_BL.append(np.trapz(w[:zi_index1+1], z[:zi_index1+1]) / (z_inv_base[-1] - z[0]))
                    qT_3km.append(np.trapz(w[zi_index2:], z[zi_index2:]) / (z[-1] - z_inv_top[-1]))
                    qT_jump.append(qT_3km[-1] - qT_BL[-1])
                except Exception:
                    qT_BL.append(np.nan)
                    qT_3km.append(np.nan)
                    qT_jump.append(np.nan)
                    
                try:
                    if not np.isnan(zb_index) and zb_index > 0:
                        tlbl=np.trapz(theta[:zb_index + 1], z[:zb_index + 1]) / (z[zb_index] - z[0])
                        thetaL_BL.append(tlbl)
                    else:
                        thetaL_BL.append(np.trapz(theta[:zi_index1+1], z[:zi_index1+1]) / (z_inv_base[-1] - z[0]))
                except Exception:
                    thetaL_BL.append(np.nan)
                try:        
                    poly_coeffs = np.polyfit(z[zi_index2:], theta[zi_index2:], 1)
                    thetaL_3km.append(np.polyval(poly_coeffs, 3000))
                    thetaL_jump.append(np.polyval(poly_coeffs, z_inv_base[-1]) - thetaL_BL[-1])
                except Exception:
                    thetaL_3km.append(np.nan)
                    thetaL_jump.append(np.nan)
    
                # Delta values for decoupling as in Jones
                if not np.isnan(zi_index1) and zi_index1 > 5:
                    try:
                        n2 = round(zi_index1 / 3)
                        q_bot.append(np.trapz(w[:n2+1], z[:n2+1]) / (z[n2] - z[0]))
                        q_up.append(np.trapz(w[zi_index1 - n2:zi_index1], z[zi_index1 - n2:zi_index1]) / (z[zi_index1-1] - z[zi_index1 - n2]))
                        theta_bot.append(np.trapz(theta[:n2+1], z[:n2+1]) / (z[n2] - z[0]))
                        theta_up.append(np.trapz(theta[zi_index1 - n2:zi_index1], z[zi_index1 - n2:zi_index1]) / (z[zi_index1-1] - z[zi_index1 - n2]))
                        dq_decoupling.append(q_bot[-1] - q_up[-1])
                        dtheta_decoupling.append(theta_bot[-1] - theta_up[-1])
                    except Exception:
                        q_bot.append(np.nan)
                        q_up.append(np.nan)
                        theta_bot.append(np.nan)
                        theta_up.append(np.nan)
                        dq_decoupling.append(np.nan)
                        dtheta_decoupling.append(np.nan)
                else:
                    q_bot.append(np.nan)
                    q_up.append(np.nan)
                    theta_bot.append(np.nan)
                    theta_up.append(np.nan)
                    dq_decoupling.append(np.nan)
                    dtheta_decoupling.append(np.nan)
    
                # BL wind
                try:
                    BLwnd_spd_avg.append(np.trapz(windspeed[:zi_index1], z[:zi_index1]) / (z[zi_index1-1] - z[0]))
                    u = windspeed * np.sin(np.deg2rad(winddir))
                    BLu = np.trapz(u[:zi_index1], z[:zi_index1]) / (z[zi_index1-1] - z[0])
                    v = windspeed * np.cos(np.deg2rad(winddir))
                    BLv = np.trapz(v[:zi_index1], z[:zi_index1]) / (z[zi_index1-1] - z[0])
                    BLwnd_dir_avg.append(np.mod(np.arctan2(BLu, BLv) * 180 / np.pi, 360))
                except Exception:
                    BLwnd_spd_avg.append(np.nan)
                    BLwnd_dir_avg.append(np.nan)
    
                ind += 1
                
#%% Post processing                

# Create table
wellmixed_dtv = np.array(dtv1) < 0.25
decoupled_dtv = np.array(dtv1) > 1
wellmixed_dLCL = np.array(dLCL_decoupling) < 150
decoupled_dLCL = np.array(dLCL_decoupling) > 150

# Create a pandas DataFrame
ALL_SOUNDINGS = pd.DataFrame({
'date_i': date_i,
'n_clouds': n_clouds,
'LCL_srf': LCL_srf,
'z_cloudbase': z_cloudbase,
'dLCL_decoupling': dLCL_decoupling,
'z_inv_base': z_inv_base,
'z_inv_top': z_inv_top,
'sfc_inversion': sfc_inversion,
'top_inversion': top_inversion,
'qT_BL': qT_BL,
'qT_jump': qT_jump,
'qT_3km': qT_3km,
'thetaL_BL': thetaL_BL,
'thetaL_jump': thetaL_jump,
'thetaL_3km': thetaL_3km,
'BLwnd_dir_avg': BLwnd_dir_avg,
'BLwnd_spd_avg': BLwnd_spd_avg,
'dq_decoupling': dq_decoupling,
'dtheta_decoupling': dtheta_decoupling,
'wellmixed_dtv': wellmixed_dtv,
'decoupled_dtv': decoupled_dtv,
'wellmixed_dLCL': wellmixed_dLCL,
'decoupled_dLCL': decoupled_dLCL
})

# Add metadata (similar to MATLAB's table properties)
ALL_SOUNDINGS.attrs['Description'] = 'Sounding parameters for NKX, May-Sept 2014-2017'
ALL_SOUNDINGS.attrs['VariableUnits'] = ['UTC time', '', 'm', 'm', 'm', 'm', 'm', '', '', 'g/kg', 'g/kg', 'g/kg', 'K', 'K', 'K', '°', 'm/s', 'g/kg', 'K', 'K', 'K', 'm', 'm']
ALL_SOUNDINGS.attrs['VariableDescriptions'] = {
'date_i': 'Sounding Time',
'n_clouds': 'Number of cloud layer from the RH profile (RH>95%). Cases with 0 clouds are not wanted',
'LCL_srf': 'LCL computed using the Bolton formula',
'z_cloudbase': 'Cloud base of the last layer of cloud based on the RH profile',
'dLCL_decoupling': 'Difference between LCL and the cloud base. Measure of decoupling following Jones',
'z_inv_base': 'Inversion base height',
'z_inv_top': 'Inversion top height',
'sfc_inversion': '1 if the inversion found is at the surface (not wanted)',
'top_inversion': '1 if the inversion found is at the top of the profile (not wanted)',
'qT_BL': 'Total water mixing ratio in the boundary layer. Well mixed assumption',
'qT_jump': 'Total water mixing ratio jump at the top of the boundary layer. Well mixed assumption',
'qT_3km': 'Total water mixing ratio at 3km. Well mixed assumption',
'thetaL_BL': 'Liquid water potential temperature in the boundary layer. Well mixed assumption',
'thetaL_jump': 'Liquid water potential temperature jump at the top of the boundary layer. Well mixed assumption',
'thetaL_3km': 'Liquid water potential temperature at 3km. Well mixed assumption',
'BLwnd_dir_avg': 'Direction of the boundary layer wind velocity average',
'BLwnd_spd_avg': 'Speed of the boundary layer wind velocity average',
'dq_decoupling': 'Δq between bottom and top of the BL. Measure for decoupling following Ghate',
'dtheta_decoupling': 'Δθv between bottom and top of the BL.´ Measure for decoupling following Ghate',
'wellmixed_dtv': '1: it is probably well-mixed since Δθv<0.25',
'decoupled_dtv': '1: it is probably decoupled since Δθv>1',
'wellmixed_dLCL': '1: it is probably well-mixed since ΔLCL<150',
'decoupled_dLCL': '1: it is probably decoupled since ΔLCL>150'
}

# Save the DataFrame to a CSV file
ALL_SOUNDINGS.to_csv('NKX_soundings_table_v1K.csv', index=False)


#%% Match SD event days
import csv

with open('../../Sc_Antofa_SD/DiasSD.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))
    
dts=[datetime(int(dd[0]),int(dd[1]),int(dd[2]),12,0,0) for dd in data]
iid=np.zeros(len(dts))
for dt in dts:
    try: 
        iid[dts.index(dt)]=np.where(dt==ALL_SOUNDINGS.date_i)[0][0]
    except Exception:
        iid[dts.index(dt)]=np.nan
        
sd_iid=iid[~np.isnan(iid)]
sd_days=ALL_SOUNDINGS.iloc[sd_iid]
