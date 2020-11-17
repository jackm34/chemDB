import numpy as np
import pandas as pd

chem_df = pd.read_excel("Z:\\Python Projects\\chemDB\\chemlog.xlsx", engine="openpyxl")
print(f"First material on list: {chem_df.iloc[0].name}, type = {type(chem_df.iloc[0].name)}")
chem_df.head()


def andrews(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W):
    Ac3_celcius = 910 - (203*C**(1/2)) + 44.7*Si - 15.2*Ni + 31.5*Mo + 104.4*C + 13.1*W
    Ac3_fahrenheit = Ac3_celcius *(9/5) +32
    return Ac3_fahrenheit

def hougardy(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W):
    Ac3_celcius = 902 - 255*C +19*Si - 11*Mn - 5*Cr + 13*Mo - 20*Ni +55*V
    Ac3_fahrenheit = Ac3_celcius *(9/5) +32
    return Ac3_fahrenheit

def kasatkin(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W):
    Ac3_celcius = 912 - 370*C - 27.4*Mn + 27.3*Si - 6.35*Cr - 32.7*Ni + 95.2*V + 190*Ti + 72*Al + 65.6*Nb + 5.57*W + 332*S + 276*P + 485*N - 900*B + 16.2*C*Mn + 32.3*C*Si + 15.4*C*Cr + 48*C*Ni + 4.32*Si*Cr-17.3*Si*Mo + 18.6*Si*Ni + 4.8*Mn*Ni + 40.5*Mo*V + 174*C**2 + 2.46*Mn**2 - 6.86*Si**2 + .322*Cr**2 + 9.9*Mo**2 + 1.24*Ni**2 - 60.2*V**2
    Ac3_fahrenheit = Ac3_celcius *(9/5) + 32
    return Ac3_fahrenheit

def TDS(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W):
    Ac3_celcius = 973-224.5*C*(1/2)-17*Mn+34*Si-14*Ni+21.6*Mo+41.8*V-20*Cu
    Ac3_fahrenheit = Ac3_celcius *(9/5) + 32
    return Ac3_fahrenheit


def PS(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W):
    MS_celcius = 489.9 - 316.7*C - 33.4*Mn - 27.8*C - 16/7*Ni -11.1*(Si+Mo+W)
    MS_fahrenheit = MS_celcius *(9/5) +32
    return MS_fahrenheit

def GS(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W):
    MS_celcius = 537.8 - 361.1*C - 38.9*(Mn+Cr)-19.4*Ni-27.8*Mo
    MS_fahrenheit = MS_celcius *(9/5) +32
    return MS_fahrenheit

def MS_L_andrews(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W):
    MS_celcius = 539 - 423*C - 30.4*Mn - 17.7*Ni - 12.1*Cr - 7.5*Mo
    MS_fahrenheit = MS_celcius *(9/5) +32
    return MS_fahrenheit

def MS_NL_andrews(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W):
    MS_celcius = 512 - 453*C - 16.9*Ni + 15*Cr - 9.5*Mo +217*C**2 - 71.5*C*Mn - 67.6*C*Cr
    MS_fahrenheit = MS_celcius *(9/5) +32
    return MS_fahrenheit

def Wang(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W):
    MS_celcius = 545 - 470.4*C - 2.96*Si - 27.7*Mn - 21.5*Cr + 38.9*Mo
    MS_fahrenheit = MS_celcius *(9/5) +32
    return MS_fahrenheit



col_C = chem_df.loc[:, "C_min":"C_max"]
col_Mn = chem_df.loc[:, "Mn_min":"Mn_max"]
col_P = chem_df.loc[:, "P_min":"P_max"]
col_S = chem_df.loc[:, "S_min":"S_max"]
col_Si = chem_df.loc[:, "Si_min":"Si_max"]
col_Ni = chem_df.loc[:, "Ni_min":"Ni_max"]
col_Cr = chem_df.loc[:, "Cr_min":"Cr_max"]
col_Mo = chem_df.loc[:, "Mo_min":"Mo_max"]
col_Cu = chem_df.loc[:, "Cu_min":"Cu_max"]
col_Sn = chem_df.loc[:, "Sn_min":"Sn_max"]
col_V = chem_df.loc[:, "V_min":"V_max"]
col_Al = chem_df.loc[:, "Al_min":"Al_max"]
col_N = chem_df.loc[:, "N_min":"N_max"]
col_B = chem_df.loc[:, "B_min":"B_max"]
col_Ti = chem_df.loc[:, "Ti_min":"Ti_max"]
col_Nb = chem_df.loc[:, "Nb_min":"Nb_max"]
col_W = chem_df.loc[:, "W_min":"W_max"]

chem_df['C_ave'] = col_C.mean(axis=1)
chem_df['Mn_ave'] = col_Mn.mean(axis=1)
chem_df['P_ave'] = col_P.mean(axis=1)
chem_df['S_ave'] = col_S.mean(axis=1)
chem_df['Si_ave'] = col_Si.mean(axis=1)
chem_df['Ni_ave'] = col_Ni.mean(axis=1)
chem_df['Cr_ave'] = col_Cr.mean(axis=1)
chem_df['Mo_ave'] = col_Mo.mean(axis=1)
chem_df['Cu_ave'] = col_Cu.mean(axis=1)
chem_df['Sn_ave'] = col_Sn.mean(axis=1)
chem_df['V_ave'] = col_V.mean(axis=1)
chem_df['Al_ave'] = col_Al.mean(axis=1)
chem_df['N_ave'] = col_N.mean(axis=1)
chem_df['B_ave'] = col_B.mean(axis=1)
chem_df['Ti_ave'] = col_Ti.mean(axis=1)
chem_df['Nb_ave'] = col_Nb.mean(axis=1)
chem_df['W_ave'] = col_W.mean(axis=1)

for index, material in enumerate(chem_df['material'].tolist()):
    C = chem_df['C_ave'].iloc[index]
    Mn = chem_df['Mn_ave'].iloc[index]
    P = chem_df['P_ave'].iloc[index]
    S = chem_df['S_ave'].iloc[index]
    Si = chem_df['Si_ave'].iloc[index]
    Ni = chem_df['Ni_ave'].iloc[index]
    Cr = chem_df['Cr_ave'].iloc[index]
    Mo = chem_df['Mo_ave'].iloc[index]
    Cu = chem_df['Cu_ave'].iloc[index]
    Sn = chem_df['Sn_ave'].iloc[index]
    V = chem_df['V_ave'].iloc[index]
    Al = chem_df['Al_ave'].iloc[index]
    N = chem_df['N_ave'].iloc[index]
    B = chem_df['B_ave'].iloc[index]
    Ti = chem_df['Ti_ave'].iloc[index]
    Nb = chem_df['Nb_ave'].iloc[index]
    W = chem_df['W_ave'].iloc[index]

    chem_df.loc[index, 'andrews'] = andrews(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W)
    chem_df.loc[index, 'hougardy'] = hougardy(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W)
    chem_df.loc[index, 'kasatkin'] = kasatkin(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W)
    chem_df.loc[index, 'TDS'] = TDS(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W)

    chem_df.loc[index, 'PS'] = PS(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W)
    chem_df.loc[index, 'GS'] = GS(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W)
    chem_df.loc[index, 'MS_L_andrews'] = MS_L_andrews(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W)
    chem_df.loc[index, 'MS_NL_andrews'] = MS_NL_andrews(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W)
    chem_df.loc[index, 'Wang'] = Wang(C, Mn, P, S, Si, Ni, Cr, Mo, Cu, Sn, V, Al, N, B, Ti, Nb, W)

print(chem_df)
