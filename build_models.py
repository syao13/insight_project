# data analysis and wrangling
import os
import sys
import time
import json

import multiprocessing
import tempfile
import pandas as pd
import numpy as np
import random as rnd
import pickle

# machine learning
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.metrics import explained_variance_score
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import median_absolute_error
from sklearn.model_selection import train_test_split

# import data
data = pd.read_csv("sanger1018_brainarray_ensemblgene_rma.txt", sep='\t')
cellline = pd.read_excel("Cell_Lines_Details.xlsx")
dose = pd.read_excel("v17.3_fitted_dose_response.xlsx")

### Test some different models using 5-fold cross validation
## one drug at a time
drug_ids = dose.DRUG_ID.value_counts().index.tolist()
drug_counts = dose.DRUG_ID.value_counts()
#drug_id = 1494

def test_models(drug_id):
    onedrug_name = dose.loc[dose['DRUG_ID'] == drug_id]['DRUG_NAME'].tolist()[0]
    onedrug_dose = dose.loc[dose.DRUG_ID == drug_id]

    # select all cell lines that were tested on the drug
    # select and sort rnaseq data by cell line order
    onedrug_dose = dose.loc[dose.DRUG_ID == drug_id]
    onedrug_ind = [str(x) for x in set(onedrug_dose.COSMIC_ID) if str(x) in data.columns and x in cellline['COSMIC identifier'].tolist()]

    onedrug_cellline = cellline[cellline['COSMIC identifier'].isin(onedrug_ind)]
    onedrug_data = data[['ensembl_gene'] + [i for i in onedrug_cellline['COSMIC identifier'].astype(str).tolist()]]
    onedrug_dose = onedrug_dose[onedrug_dose['COSMIC_ID'].isin(onedrug_ind)]
    onedrug_dose['sort'] = pd.Categorical(
        onedrug_dose['COSMIC_ID'].astype(str).tolist(),
        categories=onedrug_data.columns.tolist(),
        ordered=True
    )
    onedrug_dose = onedrug_dose.sort_values('sort')


    temp = onedrug_cellline['Cancer Type\n(matching TCGA label)'].astype(str) + onedrug_cellline['Screen Medium'].astype(str)
    stratified_category = temp.replace(temp.value_counts().index[temp.value_counts() == 1], ['nanR'] * np.sum(temp.value_counts() == 1))
    
    ## First random forest
    X = onedrug_data.drop(['ensembl_gene'], axis=1).T
    y = np.array(onedrug_dose['LN_IC50'].tolist())
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1, stratify=stratified_category) # test size 20%

    model = RandomForestRegressor(max_depth=2, random_state=0, n_estimators=100)
    model.fit(X_train, y_train)
    
    ## save important genes
    importance_idx = np.argsort(model.feature_importances_)
    important_genes = {onedrug_name: onedrug_data['ensembl_gene'][importance_idx[-50:]].tolist()}
    filepath = './conf/'+ onedrug_name + '.json'
    with open(filepath, 'w') as outfile:
        json.dump(important_genes, outfile)
    
    ## second random forest
    X_train_subset = X_train.iloc[:, importance_idx[-50:]]
    X_test_subset = X_test.iloc[:, importance_idx[-50:]]
    model = RandomForestRegressor(max_depth=5, random_state=0, n_estimators=200)
    model.fit(X_train_subset, y_train)
    
    ## save the model
    filepath = './models/'+ onedrug_name + '_model.sav'
    pickle.dump(model, open(filepath, 'wb'))

    baseline = np.sum((y_test - np.median(y_test)) ** 2) / y_test.shape[0]
    results = [drug_id, onedrug_name, drug_counts[drug_id], baseline]

    y_pred_train = model.predict(X_train_subset)
    results = results + [mean_squared_error(y_pred_train, y_train),
                         mean_absolute_error(y_pred_train, y_train),
                         median_absolute_error(y_pred_train, y_train),
                         r2_score(y_pred_train, y_train),
                         explained_variance_score(y_pred_train, y_train)]
    
    y_pred = model.predict(X_test_subset)
    results = results + [mean_squared_error(y_pred, y_test),
                         mean_absolute_error(y_pred, y_test),
                         median_absolute_error(y_pred, y_test),
                         r2_score(y_pred, y_test),
                         explained_variance_score(y_pred, y_test)]

    return results

def processFunction(drug_id):
    results = test_models(drug_id)
    
    if results:
        item = ','.join([str(j) for j in results])
        globalTempFile.write("%s\n" % item)
    return globalTempFile.name

def openTemporaryFile(temp):
    dirname = os.getcwd()
    global globalTempFile
    globalTempFile = tempfile.NamedTemporaryFile(mode='w', buffering=1, dir=dirname, prefix="tempModel_", delete=False)
    time.sleep(0.01) # sleep 0.01 seconds to prevent the same worker from calling this twice.
    return globalTempFile.name

def closeTemporaryFile(temp):
    filename = globalTempFile.name
    globalTempFile.close()
    # Sleep 0.01 seconds to prevent the same worker from calling this twice.
    # May need to increase this to 1 second or even 10 seconds to make this work.
    time.sleep(0.01)
    return filename


if __name__ == '__main__':
    resultFile = sys.argv[1]
    ids = drug_ids[0:4]

    with multiprocessing.Pool() as pool:
        result_filenames = pool.map(openTemporaryFile, (1 for i in range(multiprocessing.cpu_count())))
        pool.map(processFunction, ids)
        closed_filenames = pool.map(closeTemporaryFile, (1 for i in range(multiprocessing.cpu_count())))
    
    unclosed_filenames = set(result_filenames) - set(closed_filenames)
    if unclosed_filenames: # check whether all result files were closed.  Increase sleep time in closeTemporary File to prevent this.
        print("Unclosed Files: ", unclosed_filenames)
    
    with open(resultFile, "w") as outfile:
        for f in result_filenames:
            if f:
                with open(f, "r") as infile:
                    outfile.write(infile.read())
                os.remove(f) # remove the file

