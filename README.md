# dragenGermline

A workflow for calling SNVs on fastq inputs in germline mode

## Overview

## Dependencies

* [gsi hg38 modules : hg38-dbsnp 138](https://gitlab.oicr.on.ca/ResearchIT/modulator)


## Usage

### Cromwell
```
java -jar cromwell.jar run dragenGermline.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`sampleInputs`|Array[InputGroup]|Input structure with tumor fastq files and read group strings
`outputFileNamePrefix`|String|Prefix for output files
`reference`|String|The genome reference build. For example: hg19, hg38, mm10


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extractInfoLine.timeout`|Int|4|Timeout for the job
`extractInfoLine.jobMemory`|Int|4|Job allocated RAM
`composeList.jobMemory`|Int|4|Job allocated RAM
`composeList.timeout`|Int|4|Timeout for the job
`runDragenGermline.enableDupMarking`|Boolean|true|Flag for duplicate marking, true by  default
`runDragenGermline.enableTargeted`|Boolean|true|Flag for enabling calling on targets like HBA, GBA etc. clusters
`runDragenGermline.additionalParameters`|String?|None|Additional dragen parameters
`runDragenGermline.timeout`|Int|96|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`unfilteredVcf`|File|SNV calls before applying any filters|vidarr_label: unfilteredVcf
`filteredVcf`|File|SNV calls with filter information attached|vidarr_label: filteredVcf
`ploidyVcf`|File?|Ploidy vcf file|vidarr_label: ploidyVcf
`targetedVcf`|File?|Targeted vcf file|vidarr_label: targetedVcf


## Commands
This section lists command(s) run by dragenGermline workflow
 
* Running dragenGermline
 
dragenGermline is a workflow which launches DRAGEN SNV calling pipeline. It creates
input list based on information passed by the respective olive and then aligns
all reads using input fastq files, calling SNVs in Germline mode after that. 
It applies a number of filters and adds annotations from dbSNP database, if available 
 
### Extracting information from RG line
 
```
     python3<<CODE
     import json
     import re
     jsonInput = "~{write_json(fastqInput)}"
     with open(jsonInput, "r") as ji:
         inputData = json.load(ji)
     ji.close()
 
     try:
         myPattern = r'\S+?\:\S+'
         rgs = re.findall(myPattern, inputData['readGroup'])
         for rgroup in rgs:
             if rgroup.startswith("ID:"):
                 RGID = rgroup.split(":")[1]
                 Lane = rgroup.split("_")[-2]
             if rgroup.startswith("SM:"):
                 RGSM = rgroup.split(":")[1]
             if rgroup.startswith("LB:"):
                 RGLB = rgroup.split(":")[1]
         fastqR1 = inputData['fastqR1']
         fastqR2 = inputData['fastqR2']
         myResult = ",".join([RGID, RGSM, RGLB, Lane, fastqR1, fastqR2])
         print(myResult)
     except:
         print("Error parsing string")
     CODE 
```
### Composing input lists
 
```
    python3<<CODE
    l = "~{sep=' ' inputLines}"
    inLines = l.split()
    linesToPrint = ["RGID,RGSM,RGLB,Lane,Read1File,Read2File\n"]
    for inputString in inLines:
        inputString.rstrip()
        if not inputString.startswith("Error"):
            linesToPrint.append(inputString + "\n")
 
    with open("~{outputFileName}", "w") as tl:
        tl.writelines(linesToPrint)
    tl.close() 
    CODE
```
### Running SNV caller
 
```
       dragen -f -r ~{refDir} \
       --fastq-list ~{sampleFastqList} \
       --enable-duplicate-marking ~{enableDupMarking} \
       --enable-variant-caller true \
       --enable-targeted ~{enableTargeted} \
       --dbsnp ~{dbSNP} \
       --output-directory . \
       --output-file-prefix ~{outputFileNamePrefix} ~{additionalParameters}
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
