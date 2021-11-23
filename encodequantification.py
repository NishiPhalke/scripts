#!/usr/bin/env python

import requests
import numpy as np
from joblib import Parallel, delayed

response = requests.get('https://www.encodeproject.org/search/?type=Experiment&searchTerm=RNA-seq&assay_title=polyA%20plus%20RNA-seq&format=json&assay_title=total%20RNA-seq&limit=all')
data = response.json()
experimentAccessions = list(map(lambda x: x['accession'], data['@graph']))

def isPipeLine(file):
    if file.get('analysis_step_version') and file.get('analysis_step_version').get('analysis_step') and file.get('analysis_step_version').get('analysis_step').get('pipelines'):
        if file.get('analysis_step_version').get('analysis_step').get('pipelines')[0]['title'] == "RNA-seq pipeline (Reddy GGR)":
            return True
        else:
            return False        
    else:
        return False


def filterFiles(file):    
    if file.get('cloud_metadata') and file['status'] == "released" and file['file_type']=="tsv" and (file['output_type']=="gene quantifications" or file['output_type']=="transcript quantifications") and (not isPipeLine(file)):        
        if file['output_type']=="transcript quantifications" and 'kallisto' in file.get('aliases')[0]:            
            return False
        else:    
            if file.get('assembly') == "GRCh38" and file.get('genome_annotation') == "V29":
                return True
            elif file.get('assembly') == "hg19" and file.get('genome_annotation') == "V19":
                return True
            elif file.get('assembly') == "mm10" and file.get('genome_annotation') == "M21":
                return True
            else:
                return False    
    else:
        return False



def getData(expAccession):
    experiment =  requests.get('https://www.encodeproject.org/experiments/'+expAccession+'/?format=json')
    experimentData = experiment.json()
    files  = experimentData['files']        
    filteredFiles = list(filter(filterFiles, files))
    q = []
    for expFile in filteredFiles:
        q.append({
            "experimentAccession": expAccession,
            "fileAccession": expFile['accession'],
            "file_type": expFile['file_type'],
            "assembly": expFile['assembly'],
            "output_type": expFile['output_type'],
            "genome_annotation": expFile['genome_annotation'],
            "url": expFile.get('cloud_metadata').get("url")
        })
    return q

results = Parallel(n_jobs=5)(delayed(getData)(expAccession) for expAccession in experimentAccessions)    

fullList = results

geneQuantassemblyList =  {}
transcriptQuantassemblyList =  {}
for lst in fullList:
    for mainList in lst:
        if mainList.get('output_type')=="gene quantifications":
            if geneQuantassemblyList.get(mainList.get('assembly')):          
                geneQuantassemblyList[mainList.get('assembly')].append(mainList)
            else:
                geneQuantassemblyList[mainList.get('assembly')] = [mainList]
        else:
            if transcriptQuantassemblyList.get(mainList.get('assembly')):          
                transcriptQuantassemblyList[mainList.get('assembly')].append(mainList)
            else:
                transcriptQuantassemblyList[mainList.get('assembly')] = [mainList]

assemblyList = ['GRCh38','mm10','hg19']

for asm in assemblyList:
    genesList = geneQuantassemblyList.get(asm)
    for g in genesList:
        with open('/home/niship/encodegenequant/encode_'+asm+'_genes.txt', 'a') as gfl:
            gfl.write(g.get('experimentAccession')+'\t'+g.get('fileAccession')+'\t'+g.get('url')+'\n')


for asm in assemblyList:
    transcriptList = transcriptQuantassemblyList.get(asm)
    for t in transcriptList:
        with open('/home/niship/encodegenequant/encode_'+asm+'_transcript.txt', 'a') as tfl:
            tfl.write(t.get('experimentAccession')+'\t'+t.get('fileAccession')+'\t'+t.get('url')+'\n')

def runnumpy(assemblyList,quanttype,metric):
    for asm in assemblyList:
        allmetricvals = []
        with open('/home/niship/encodegenequant/encode_'+asm+'_'+quanttype+'.txt', 'r') as r:            
            for fline in r:
                url = fline.strip().split('\t')[2]            
                tsvContent = requests.get(url).text
                metricarr = []
                lines = tsvContent.split('\n')
                header = lines[0]
                headervals = list(map(lambda x: x.lower(),header.split('\t')))
                if metric.lower() in headervals:
                    metricindex = headervals.index(metric.lower())
                    if len(lines) > 1: 
                        for l in lines[1:]:
                            vals = l.split('\t')
                            if len(vals) > 1:							
                                metricarr.append(np.float32(vals[metricindex]))
                print(len(metricarr),url,len(lines),fline.strip().split('\t')[0],fline.strip().split('\t')[1])
                allmetricvals.append(metricarr)					
        allmetricvals_tr = np.transpose(allmetricvals)
        with open('/home/niship/encodegenequant/encode_'+asm+'_'+metric+'_'+quanttype+'.npy', 'wb') as f:
            np.save(f, np.array(allmetricvals_tr,order='C'), allow_pickle = True)
        print('done creating numpy matrix for ', asm , quanttype, metric)

#Run
runnumpy(['mm10','hg19','GRCh38'],'genes','tpm')        
runnumpy(['mm10','hg19','GRCh38'],'genes','expected_count')      

runnumpy(['GRCh38','mm10', 'hg19'],'transcript','tpm')      
runnumpy(['mm10','hg19','GRCh38'],'transcript','expected_count')        