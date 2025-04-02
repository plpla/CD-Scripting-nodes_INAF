import json
import csv
import os
import sys
from CdScriptingNodeHelper import ScriptingResponse, print_error

import pandas as pd
import numpy as np

# Author: Pier-Luc Plante
# e-mail: plpla2 @ ulaval.ca

#
# Scripting Node for Compound Discoverer v3.4 that calculates QC metrics and add them to the Compounds table.
# CDScriptingNodeHelper was provided by CD and is available at mycompounddiscoverer.com
#

print('CD Scripting Node')
            
# start in development mode where nodeargs are given explicitely rather than reading it as command line argument
if sys.argv[1] == '-devel':    
    print(f'Development mode: Current Dir is {os.getcwd()}')
    nodeargs_path = 'node_args.json'        
else:
    nodeargs_path = sys.argv[1]        

# parse node args from Compound Discoverer and extract location of Compounds and Study Information table
try:
    with open(nodeargs_path, "r") as rf:
        nodeargs = json.load(rf) 
        compounds_path = ""
        studyInformation_path = ""
        response_path = nodeargs["ExpectedResponsePath"]
        tables = nodeargs["Tables"]
        for table in tables:
            if table["TableName"] == "Compounds":
                compounds_path = table['DataFile']                                    
                if table['DataFormat'] != 'CSV':
                    print_error(f"Unknown Data Format {table['DataFormat']}")    
                    exit(1)
            elif table["TableName"] == "Study Information":
                studyInformation_path = table['DataFile']
                if table['DataFormat'] != 'CSV':
                    print_error(f"Unknown Data Format {table['DataFormat']}")
                    exit(1)        
except Exception as e: 
    print_error(f"Could not read Compound Discoverer node args\n{str(e)}")
    exit(1)

if compounds_path == "" or studyInformation_path == "":
    print_error('Features file not defined in node args.')
    exit(1)

try:
    studyInformation_table = pd.read_csv(studyInformation_path, header=0, sep="\t")
    QCsamples_fileID = studyInformation_table["Study File ID"][studyInformation_table["Sample Type"] == "QualityControl"]
    compounds_table = pd.read_csv(compounds_path, header=0 , sep="\t")
    norm_area_QCs_columns = [col for col in compounds_table.columns if col.startswith("Norm Area") and col.split()[-1] in list(QCsamples_fileID)]
    raw_area_QCs_columns = [col for col in compounds_table.columns if col.startswith("Area") and col.split()[-1] in list(QCsamples_fileID)]

    if len(raw_area_QCs_columns) > 0:
        selected_raw_data = compounds_table[raw_area_QCs_columns]
        raw_qc_rsd = selected_raw_data.std(axis=1) / selected_raw_data.mean(axis=1) * 100
    if len(norm_area_QCs_columns) > 0:
        selected_normalized_data = compounds_table[norm_area_QCs_columns]
        normalized_qc_rsd = selected_normalized_data.std(axis=1) / selected_normalized_data.mean(axis=1) * 100
        n_usable_QC = selected_normalized_data.notna().sum(axis=1)


except Exception as e:
    print_error('Could not process data')
    print_error(e)        
    exit(1)
    

# Write to file
output_df = pd.DataFrame({"Compounds ID" : compounds_table["Compounds ID"], 
                        "RSD QC Areas [%]" : raw_qc_rsd, 
                        "RSD Corr. QC Areas [%]" : normalized_qc_rsd,
                        "# Usable QC" : n_usable_QC})

outfilename = 'CompoundsWithQCMetrics.txt'
(workdir, _ ) = os.path.split(response_path)

#outfile_path = os.path.join(workdir, outfilename) # This does not work on my current setup: it adds \\ to a path containing / as separator)
outfile_path = workdir + f"/{outfilename}"
with open(outfile_path, mode='w') as compounds_csv:
    output_df.to_csv(compounds_csv, sep="\t", index=False)
            
# entries for new column in Features table
response = ScriptingResponse()
response.add_table('Compounds', outfile_path)
response.add_column('Compounds', 'Compounds ID', 'Int', 'ID')
response.add_column('Compounds', 'RSD QC Areas [%]', 'Float')
response.add_column('Compounds', 'RSD Corr. QC Areas [%]', 'Float')
response.add_column('Compounds', '# Usable QC', 'Int')

# save to disk
response.save(response_path)
