import pandas as pd
import random


random.seed(1)
par = "parameters.csv"
par_cts = "parameters_constants.csv"

df_pars = pd.read_csv(par)
df_pars_cts = pd.read_csv(par_cts)

alpha = df_pars_cts.alpha[0]
k = df_pars_cts.k[0]
tau = df_pars_cts.tau[0]

pi = list(df_pars.pi)
di = list(df_pars.di)
Di = list(df_pars.Di)
ci2 = list(df_pars.ci2)
epi = [round(val, 2) for val in list(df_pars.epi)]
efi = [round(val, 2) for val in list(df_pars.efi)]
evi = [round(val, 2) for val in list(df_pars.evi)]
en = [0 for _ in range(len(efi))]
ei = [[random.choice([epi[i], efi[i], en[i], evi[i], en[i]]) for j in range(5)] for i in range(len(en))]

