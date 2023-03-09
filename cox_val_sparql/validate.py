import os 
import pandas as pd
import numpy as np
import json
import math
from sksurv.metrics import concordance_index_censored
import requests
from io import StringIO
import re

# read database uri
database_uri = os.getenv('DATABASE_URI','no_database_uri')
args = open(database_uri).read().strip()

def run_sparql_for_cox (args):
    #print(args)
    url1 = "http://gateway.docker.internal:7200/rest/repositories"
    try:
        url = "http://"+args+":7200/rest/repositories"
        node = requests.get(url, headers={"Accept": "application/json"})
        nodes = list(node.json())
    except:
        node = requests.get(url1, headers={"Accept": "application/json"})
        nodes = list(node.json())
    for i in nodes:
        txt = i['id']
        x = re.search("userRepo", txt)
        if x:
            repo = txt
    query = """
        PREFIX dbo: <http://um-cds/ontologies/databaseontology/>
        PREFIX roo: <http://www.cancerdata.org/roo/>
        PREFIX ncit: <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        SELECT ?data ?Tstage ?Treatment ?overall_survival_in_days ?HPV ?Nstage ?event_overall_survival
        WHERE { 
    	?tablerow dbo:has_column ?dataset.
    	?dataset rdf:type ncit:C25474.  
    	?dataset dbo:has_cell ?datavalue.
    	?datavalue dbo:has_value ?data.
        ?tablerow roo:P100231 ?chemov. 
        ?tablerow roo:P100029 ?neoplasm.
        ?neoplasm roo:P100202 ?tumourv.
        ?neoplasm roo:P100244 ?tstagev.
        ?tablerow roo:P100254 ?survivalv.

        ?chemov dbo:has_cell ?chemocell.
        ?tumourv dbo:has_cell ?tumourcell.
        ?tstagev dbo:has_cell ?tcell.
        ?survivalv dbo:has_cell ?scell.

        ?tumourcell a ?t.
        ?tcell a ?tS.
        ?scell a ?s.
        ?chemocell a ?c.

        FILTER regex(str(?t), ("http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C12762"))
        FILTER regex(str(?tS), ("http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C48737|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C48719|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C48720|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C48724|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C48728|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C48732"))
        FILTER regex(str(?s), ("http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C28554|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C37987"))
        FILTER regex(str(?c), ("http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C94626|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C15313"))
        BIND(strafter(str(?t), "http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#") AS ?TumourLocation)
        BIND(strafter(str(?tS), "http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#") AS ?Tstage)
        BIND(strafter(str(?s), "http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#") AS ?event_overall_survival)
        BIND(strafter(str(?c), "http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#") AS ?Treatment)
        ?tablerow roo:has ?survivald.
        ?survivald dbo:has_cell ?sdcell.
        ?sdcell roo:P100042 ?overall_survival_in_days.
    OPTIONAL {
        ?tablerow roo:P100022 ?hpvv.
        ?hpvv dbo:has_cell ?hpvcell.
        ?hpvcell a ?h.
        FILTER regex(str(?h), ("http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C128839|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C131488|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C10000"))
        BIND(strafter(str(?h), "http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#") AS ?HPV)
       }
    OPTIONAL {
        ?tablerow roo:P100029 ?neoplasm.
        ?neoplasm roo:P100242 ?nstagev.
        ?nstagev dbo:has_cell ?cell.
        ?cell a ?n.
        FILTER regex(str(?n), ("http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C48705|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C48706|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C48786"))
        BIND(strafter(str(?n), "http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#") AS ?Nstage)
    }
    }
    """
    codedict = {
        "C16576": "Female", "C20197": "Male", "C27966": "Stage I", "C28054": "Stage II",
        "C27970": "Stage III", "C27971": "Stage IV", "C12762": "Oropharynx", "C12420": "Larynx",
        "C12246": "Hypopharynx", "C12423": "Nasopharynx", "C00000": "Unknown", "C48737": "Tx", "C48719": "T0", "C48720": "T1", "C48724": "T2",
        "C48728": "T3", "C48732": "T4", "C48705": "N0", "C48706": "N1", "C48786": "N2", "C48714": "N3",
        "C48704": "Mx", "C48699": "M0", "C48700": "M1", "C28554": "1", "C37987": "0", "C128839": "HPV_Positive",
        "C131488": "HPV_Negative", "C10000": "Unknown", "C94626": "ChemoRadiotherapy", "C15313": "Radiotherapy"
    }

    endpoint = "http://"+args+":7200/repositories/" + repo
    annotationResponse = requests.post(endpoint,data="query=" + query,
                                       headers={
                                           "Content-Type": "application/x-www-form-urlencoded",
                                           # "Accept": "application/json"
                                       })
    output = annotationResponse.text
    #return output 
    TList = ["Tstage_Tx", "Tstage_T0", "Tstage_T1", "Tstage_T2", "Tstage_T3", "Tstage_T4"]
    NList = ["Nstage_Nx", "Nstage_N0", "Nstage_N1", "Nstage_N2"]
    hnscc = pd.read_csv(StringIO(output))
    for col in hnscc.columns:
        hnscc[col] = hnscc[col].map(codedict).fillna(hnscc[col])
    hnscc.loc[hnscc["overall_survival_in_days"] >= 1826, "overall_survival_in_days"] = 1826
    hnscc.loc[hnscc["overall_survival_in_days"] >= 1826, "event_overall_survival"] = 0
    hnscc = pd.get_dummies(hnscc, columns=['Treatment', 'Tstage', 'HPV', 'Nstage'])
    for i in TList:
        if i not in hnscc:
            hnscc[i] = 0 
    for j in NList:
        if j not in hnscc:
            hnscc[j] = 0         
    hnscc['T2orLower'] = (hnscc['Tstage_T0'] | hnscc['Tstage_T1'] | hnscc['Tstage_T2'])
    hnscc['T3orHigher'] = (hnscc['Tstage_T3'] | hnscc['Tstage_T4'])   
    hnscc["T2orLower"] = hnscc["T2orLower"].astype(int)
    hnscc["T3orHigher"] = hnscc["T3orHigher"].astype(int)
    hnscc['N1orLower'] = (hnscc['Nstage_N0'] | hnscc['Nstage_N1'])
    hnscc['N2orHigher'] = (hnscc['Nstage_N2'])   
    hnscc["N1orLower"] = hnscc["N1orLower"].astype(int)
    hnscc["N2orHigher"] = hnscc["N2orHigher"].astype(int)    
    hnscc.to_csv('df.csv', index = False) 
    #print(hnscc)
    #print(hnscc.iloc[:, 9:])
    #return hnscc

run_sparql_for_cox(args)

data = pd.read_csv("df.csv")
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

cindex = result[0]
            
with open('output.txt', 'w') as f:
    f.write(json.dumps({'Concordance Index':cindex}))



