# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'

# %%
import json
import urllib.request
import pandas as pd
import requests
import os

class Experiment:
    def __init__(self, id, biosample, assay, file_name, target):
        self.id = id
        self.biosample = biosample
        self.assay = assay
        self.file_name = file_name
        self.target = target
def print_line(experiment):
    #return experiment.id+","+experiment.biosample+","+experiment.assay+","+experiment.file_name+","+experiment.target
    return experiment.id.split("/")[-2]+","+experiment.assay+","+experiment.file_name.split("/")[-2]

#url = "https://www.encodeproject.org/search/?type=Annotation&software_used.software.name=avocado&searchTerm=H1-heSC&biosample_ontology.term_name=H1&biosample_ontology.term_name=mesenchymal+stem+cell&annotation_type=imputation"
url = "https://www.encodeproject.org/search/?type=Annotation&software_used.software.name=avocado&searchTerm=H1-heSC&biosample_ontology.term_name=H1&annotation_type=imputation"
get_all = True
if(get_all):
    url += "&limit=all"
url += "&format=json"

folder = "encode_data"

if not os.path.exists(folder):
    os.makedirs(folder)

data = urllib.request.urlopen(url).read().decode()

obj = json.loads(data)


# %%
all_experiments = dict()
for item in obj["@graph"]:
    url = 'https://www.encodeproject.org'+item['@id']+'?format=json'
    item_data = urllib.request.urlopen(url).read().decode() 
    item_obj = json.loads(item_data)

    if("nuclear" in item_obj["description"]): continue
    if("cytosolic" in item_obj["description"]): continue

    if(len(item_obj['original_files']) < 1): # or original files
        file_name = ""
    elif(len(item_obj['original_files']) < 2):
        file_name = item_obj['original_files'][0]
    else:
        file_name = "MULTIPLE_FILES"

    if(len(item_obj['targets']) < 1):
        target = ""
    elif(len(item_obj['targets']) < 2):
        target = item_obj['targets'][0]["label"]
    else:
        target = "MULTIPLE_TARGETS"
    
    if("_RNA-seq" in item_obj['description']):
        target = "RNA-Seq "
        if("minus" in item_obj["description"]):
            target += "minus "
        if("plus" in item_obj["description"]):
            target += "plus "

        if("small" in item_obj["description"]):
            target += "small"
        if("total" in item_obj["description"]):
            target += "total"
        if("polyA_depleted" in item_obj["description"]):
            target += "polyA_depleted"
        elif("polyA" in item_obj["description"]):
            target += "polyA"
    if(target == ""):
        target = item_obj['assay_term_name'] + " "
        if("minus" in item_obj["description"]):
            target += "minus "
        if("plus" in item_obj["description"]):
            target += "plus "
    experiment = Experiment(item['@id'], item_obj['biosample_ontology']['term_name'], item_obj['assay_term_name'], file_name, target)
    if not(experiment.biosample in all_experiments):
        all_experiments[experiment.biosample] = dict()
    all_experiments[experiment.biosample][target] = experiment


# %%
df = pd.DataFrame.from_dict(all_experiments, orient='index').T
df_csv = df.applymap(lambda experiment: print_line(experiment))
for col in df_csv.columns:
    df_csv[[col+" ID",col+' assay',col+' file']] = df_csv[col].str.split(',',expand=True)
    df_csv = df_csv.drop([col], axis=1)
print(df_csv)


# %%
df_csv.to_csv(folder+"/data.csv")
print("Report generated. Starting to download.")


# %%
for sample in all_experiments:
    for target in all_experiments[sample]:
        experiment = all_experiments[sample][target]
        fileName = experiment.file_name.split("/")[-2]
        url = "https://www.encodeproject.org/files/"+fileName+"/@@download/"+fileName+".bigWig"
        r = requests.get(url, allow_redirects=True)

        open(folder+"/"+fileName+".bigWig", 'wb').write(r.content)


# %%



