import os 
import pandas as pd
import json
import math
from sksurv.metrics import concordance_index_censored

#database_uri = os.getenv('DATABASE_URI','no_database_uri')
#data = pd.read_csv(database_uri, index_col=0)

data = pd.read_csv("C:\\Users\\P70070487\\OneDrive - Maastro - Clinic\\Radiomics\\Clipped\\2022_hn1.csv")
data['event_overall_survival'] = data['event_overall_survival'].astype('bool')
data["lp"] = ""
lp_list=[]

mean_dict = {}
val_dict = {}
lp_val = {}

#Read input coefficients from text file
with open('input.txt') as json_file:
    list_betas = json.load(json_file)
    betas = list_betas[0]

#calculate linear predictor
lp = 0
for i, j in data.iterrows():
    for key in betas:
        val_dict[key] = j[key]
        lp_val[key] = (val_dict[key]*betas[key])
    lp = sum(lp_val.values())
    exp_lp = math.exp(lp)
    lp_list.append(exp_lp)
    
data['lp']=lp_list

#calculate concordance index
result = concordance_index_censored(data["event_overall_survival"], data["overall_survival_in_days"], data["lp"])
#print(result)
cindex = result[0]
            
with open('output.txt', 'w') as f:
    f.write(json.dumps({'Concordance Index':cindex}))



