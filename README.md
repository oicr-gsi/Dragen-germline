# dragenGermline

A workflow for calling SNVs on fastq inputs in germline mode

## Overview

## Dependencies

* [gsi hg38 modules : hg38-dbsnp 151](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi modules : dragen-scripts 0.1](https://gitlab.oicr.on.ca/ResearchIT/modulator)


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
`extractInfoLine.parsingScript`|String|"$DRAGEN_SCRIPTS_ROOT/bin/composeList.py"|Script for parsing inputs into a line
`extractInfoLine.timeout`|Int|4|Timeout for the job
`extractInfoLine.jobMemory`|Int|4|Job allocated RAM
`extractInfoLine.modules`|String|"dragen-scripts/0.1"|dependency modules
`composeList.listWritingScript`|String|"$DRAGEN_SCRIPTS_ROOT/bin/writeFile.py"|Script for writing out list of inputs
`composeList.jobMemory`|Int|4|Job allocated RAM
`composeList.timeout`|Int|4|Timeout for the job
`composeList.modules`|String|"dragen-scripts/0.1"|dependency modules
`runDragenGermline.enableDupMarking`|Boolean|true|Flag for duplicate marking, true by  default
`runDragenGermline.enableTargeted`|Boolean|true|Flag for enabling calling on targets like HBA, GBA etc. clusters
`runDragenGermline.additionalParameters`|String?|None|Additional dragen parameters
`runDragenGermline.timeout`|Int|96|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`unfilteredVcf`|File|SNV calls before applying any filters|vidarr_label: unfilteredVcf
`unfilteredIdx`|File|Index of SNV calls before applying any filters|vidarr_label: unfilteredIdx
`filteredVcf`|File|SNV calls with filter information attached|vidarr_label: filteredVcf
`filteredIdx`|File|Index of SNV calls with filter information attached|vidarr_label: filteredIdx
`ploidyVcf`|File?|Ploidy vcf file|vidarr_label: ploidyVcf
`ploidyIdx`|File?|Index of Ploidy vcf file|vidarr_label: ploidyIdx
`targetedVcf`|File?|Targeted vcf file|vidarr_label: targetedVcf
`targetedIdx`|File?|Index of Targeted vcf file|vidarr_label: targetedIdx


## Commands
This section lists command(s) run by dragenGermline workflow
 
* Running dragenGermline
 
dragenGermline is a workflow which launches DRAGEN SNV calling pipeline. It creates
input list based on information passed by the respective olive and then aligns
all reads using input fastq files, calling SNVs in Germline mode after that. 
It applies a number of filters and adds annotations from dbSNP database, if available 
 
### Extracting information from RG line

```
     python3 ~{parsingScript} -i ~{write_json(fastqInput)}
```

### Composing list of inputs for dragen
 
```
    python3 ~{listWritingScript} -o ~{outputFileName} -l "~{sep=';' inputLines}"
```

### Running dragen SNV caller in germline mode
 
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
