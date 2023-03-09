import os 
import pandas as pd
import numpy as np
import json
import math
from sksurv.metrics import concordance_index_censored
import requests
from io import StringIO
import re
import subprocess

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
    PREFIX ro: <http://www.radiomics.org/RO/>
    PREFIX ncit: <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    SELECT DISTINCT ?ID ?firstorderenergy ?orig_glrlm ?volume ?sphericity ?wav_glrlm ?event_overall_survival ?overall_survival_in_days
    WHERE
     {         
        ?tablerow roo:P100061 ?patient.
        ?patient dbo:has_cell ?patientcell.
        ?patientcell dbo:has_value ?ID.     
        ?tablerow ro:P00088 ?energy. 
        ?energy a ?evar.
        ?evar owl:equivalentClass ro:N8CA.
        ?energy dbo:has_cell ?energycell.
        ?energycell dbo:has_value ?firstorderenergy.
        ?tablerow ro:P00088 ?orig_g.
        ?tablerow a ?data.
        ?data owl:equivalentClass ncit:C16960.
        ?orig_g a ?ogvar.
        ?ogvar owl:equivalentClass ro:R5YN.
        ?orig_g dbo:has_cell ?orig_gcell.
        ?orig_gcell dbo:has_value ?orig_glrlm.
        ?tablerow ro:P00088 ?vol.
        ?vol a ?volvar.
        ?volvar owl:equivalentClass ro:RNU0.
        ?vol dbo:has_cell ?volcell.
        ?volcell dbo:has_value ?volume.
        ?tablerow ro:P00088 ?sphere.
        ?sphere a ?svar.
        ?svar owl:equivalentClass ro:QCFX.
        ?sphere dbo:has_cell ?spherecell.
        ?spherecell dbo:has_value ?sphericity.
        ?tablerow ro:P00039 ?waveletfilter. 
        ?waveletfilter a ?filter.
        ?filter owl:equivalentClass ro:4BEG.
        ?waveletfilter ro:P00088 ?wav_g.
        ?wav_g dbo:has_cell ?wav_gcell.
        ?wav_gcell dbo:has_value ?wav_glrlm.
    
        ?table roo:P100061 ?patient1.
        ?patient1 dbo:has_cell ?patientcell1.
        ?patientcell1 dbo:has_value ?ID.
        ?table roo:P100029 ?neoplasm.
        ?neoplasm roo:P100202 ?tumour.
        ?tumour dbo:has_cell ?tcell.
        ?tcell a ?t.
        FILTER regex(str(?t), ("http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C12762"))
        ?table roo:P100254 ?survivalv.
        ?survivalv dbo:has_cell ?scell.
        ?scell a ?s.
        FILTER regex(str(?s), ("http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C28554|http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C37987"))
        BIND(strafter(str(?s), "http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#") AS ?event_overall_survival)
        ?table roo:has ?survivald.
        ?survivald dbo:has_cell ?sdcell.
        ?sdcell roo:P100042 ?overall_survival_in_days.
    }
    """
    codedict = {
        "C12762": "Oropharynx", "C12420": "Larynx",
        "C12246": "Hypopharynx", "C12423": "Nasopharynx", 
        "C00000": "Unknown", "C28554": "1", "C37987": "0"
    }
    endpoint = "http://"+args+":7200/repositories/" + repo
    annotationResponse = requests.post(endpoint,data="query=" + query,
                                       headers={
                                           "Content-Type": "application/x-www-form-urlencoded",
                                           # "Accept": "application/json"
                                       })
    output = annotationResponse.text
    #return output 
    hnscc = pd.read_csv(StringIO(output))
    for col in hnscc.columns:
        hnscc[col] = hnscc[col].map(codedict).fillna(hnscc[col])
    hnscc.loc[hnscc["overall_survival_in_days"] >= 1826, "overall_survival_in_days"] = 1826
    hnscc.loc[hnscc["overall_survival_in_days"] >= 1826, "event_overall_survival"] = 0  
    hnscc.to_csv('df.csv', index = False) 
    #print(hnscc)
    #print(hnscc.iloc[:, 9:])
    #return hnscc

run_sparql_for_cox(args)

data = pd.read_csv("df.csv")
data['event_overall_survival_bool'] = ""
data['event_overall_survival_bool'] = data['event_overall_survival'].astype('bool')
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
result = concordance_index_censored(data["event_overall_survival_bool"], data["overall_survival_in_days"], data["lp"])
cindex = result[0]
            
with open('output.txt', 'w') as f:
    f.write(json.dumps({'Concordance Index':cindex}))