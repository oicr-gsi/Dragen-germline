version 1.0

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroup
}

struct GenomeResources {
    String dbSNP
    String referenceDirectory
}

workflow dragenGermline {
    input {
        Array[InputGroup] sampleInputs
        String outputFileNamePrefix
        String reference
    }

    parameter_meta {
        sampleInputs: "Input structure with tumor fastq files and read group strings"
        outputFileNamePrefix: "Prefix for output files"
        reference: "The genome reference build. For example: hg19, hg38, mm10"
    }

    Map[String,GenomeResources] dragen_resources_by_genome = { 
    "hg38": {
      "dbSNP": "/.mounts/labs/gsiprojects/gsi/Dragen/reference/dbSNP.151/common_all_dbSNP151_hg38p7_sorted.vcf.gz",
      "referenceDirectory": "/.mounts/labs/gsiprojects/gsi/Dragen/reference/hg38fa.p12_v4.5"
      }
    }

    String dragen_ref = dragen_resources_by_genome [ reference ].referenceDirectory
    String dragen_dbsnp = dragen_resources_by_genome [ reference ].dbSNP

    meta {
        author: "Peter Ruzanov, Lawrence Heisler"
        email: "pruzanov@oicr.on.ca, lheisler@oicr.on.ca"
        description: "A workflow for calling SNVs and CNVs on fastq inputs in germline mode"
        dependencies: [
        {
          name: "gsi hg38 modules : hg38-dbsnp/138",
          url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
        }]
        output_meta: {
          mappingMetrics: {
             description: "Mapping metrics",
             vidarr_label: "mappingMetrics"
           },
          fastqcMetrics: {
             description: "FastQC metrics",
             vidarr_label: "fastqcMetrics"
           },
          coverageMetrics: {
             description: "Coverage metrics",
             vidarr_label: "coverageMetrics"
           },
          vcVcf: {
             description: "SNV calls",
             vidarr_label: "variantCalls"
           },
          vcMetrics: {
             description: "SNV call metrics",
             vidarr_label: "variantCallMetrics"
           },
          hardfilteredVcf: {
             description: "HardFiltered variant calls",
             vidarr_label: "hardfilteredCalls"
           },
          targetedVcf: {
             description: "Targeted variant calls",
             vidarr_label: "targetedCalls"
           },
          ploidyVcf: {
             description: "Ploidy metrics",
             vidarr_label: "ploidyMetrics"
           },
          cnvVcf: {
            description: "Copy Number calls",
            vidarr_label: "cnvCalls"
          },
          cnvMetrics: {
            description: "Copy Number metrics",
            vidarr_label: "cnvCallMetrics"
          },
          svVcf: {
            description: "Structural Variant calls",
            vidarr_label: "svCalls"
          },
          svMetrics: {
            description: "Structural Variant metrics",
            vidarr_label: "svCallMetrics"
          }



        }
    }

    # Compose input lines using read group data
    scatter(s in sampleInputs) {
      call extractInfoLine {
        input:
        fastqInput = object{fastqR1: s.fastqR1, fastqR2: s.fastqR2, readGroup: s.readGroup}
      }
    }

    # file writing jobs for Normal and tumor
    call composeList {
      input:
        inputLines = extractInfoLine.outputLine,
        outputFileName = "sample_inputs.csv"
    } 

    call runDragenGermline {
      input:
        sampleFastqList = composeList.inputList,
        refDir = dragen_ref,
        dbSNP = dragen_dbsnp,
        outputFileNamePrefix = outputFileNamePrefix
    } 

    output {
        File mappingMetrics = runDragenGermline.mappingMetrics
        File fastqcMetrics = runDragenGermline.fastqcMetrics
        File coverageMetrics = runDragenGermline.coverageMetrics
        File? vcVcf = runDragenGermline.vcVcf
        File? vcMetrics = runDragenGermline.vcMetrics
        File? hardfilteredVcf =runDragenGermline.hardfilteredVcf
        File? targetedVcf = runDragenGermline.targetedVcf
        File? ploidyVcf = runDragenGermline.ploidyVcf
        File? cnvVcf = runDragenGermline.cnvVcf 
        File? cnvMetrics = runDragenGermline.cnvMetrics
        File? svVcf = runDragenGermline.svVcf
    }
}

# =====================================================================
# A scripted extraction of info from RG line to dragen-compliant string
# =====================================================================
task extractInfoLine {
   input {
       InputGroup fastqInput
       Int timeout = 4
       Int jobMemory = 4
   }

   parameter_meta {
     fastqInput: "InputGroup struct entry with fastq files"
     timeout: "Timeout for the job"
     jobMemory: "Job allocated RAM"
   }

   command <<<
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
   >>>

   runtime {
     timeout: "~{timeout}"
     memory:  "~{jobMemory} GB"
   }

   output {
     String outputLine = read_string(stdout())
   }

   meta {
     output_meta: {
       outputLine: "Output line to use in a list of fastq files in dragen-compliant format"
     }
   }
}

# =====================================================================
#  Compose gragen-comliant list of inputs to use with dragen snv caller
# =====================================================================
task composeList {
   input  {
      Array[String] inputLines
      String outputFileName
      Int jobMemory = 4
      Int timeout = 4
   }

   parameter_meta {
     inputLines: "Array of input lines to print"
     outputFileName: "Name of an output file, list of inputs"
     jobMemory: "Job allocated RAM"
     timeout: "Timeout for the job"
   }

   command<<<
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
   >>>
   

   runtime {
      timeout: "~{timeout}"
      memory:  "~{jobMemory} GB"
   }

   output {
     File inputList = "~{outputFileName}"
   }

   meta {
     output_meta: {
       inputList: "Output file to use with dragen SNV caller"
     }
   }
}
# ================================================================
# Main task for generating SNV calls in germline mode (DRAGEN mode)
#
# we need CSV files with a header and data lines organized as:
#
# RGID Read Group
# RGSM Sample ID
# RGLB Library
# Lane Flow cell lane
# Read1File - Full path to a valid FASTQ input file
# Read2File - Full path to a valid FASTQ input file. Required for paired-end input. If not using paired-end input, leave empty.
# Each FASTQ file can only be referenced once in the CSV list.
# All values in the Read2File column must be reference valid files or must all be empty.
# See https://support-docs.illumina.com/SW/dragen_v42/Content/SW/DRAGEN/TargetedCalling.htm for targeted calling option
# =====================================================================================================================

task runDragenGermline {
    input {
        File sampleFastqList
        Boolean enableDupMarking = true
        Boolean enableTargeted = true
        Boolean enableVariantCaller = true
        Boolean enableCnv = false
        Boolean enableCnvSelfNormalization = false
        Boolean enableSv = false
        String refDir
        String? additionalParameters
        String outputFileNamePrefix
        String dbSNP
        Int timeout = 96
        Int jobMemory = 96
    }

    parameter_meta {
        sampleFastqList: "List of tumor fastq files, required input"
        enableDupMarking: "Flag for duplicate marking, true by  default"
        enableTargeted: "Flag for enabling calling on targets like HBA, GBA etc. clusters"
        enableVariantCaller: "Flag for enabling variant calling"
        enableCnv: "Flag for enabling CNV calling"
        enableCnvSelfNormalization : "Flag to enable selfnormalization, required for uniformity of coverage metric"
        enableSv : "Flag to enable SV calling"
        refDir: "The reference genome directory"
        additionalParameters: "Additional dragen parameters"
        dbSNP: "Path to the dbSNP reference file"
        outputFileNamePrefix: "Output file name prefix"
        timeout: "Hours before task timeout"
    }
    

    command <<<
      ### add dragen 4.5.4 to the path
      export PATH=$PATH:/opt/dragen/4.5.4/bin/

      dragen -f -r ~{refDir} \
      --fastq-list ~{sampleFastqList} \
      --enable-duplicate-marking ~{enableDupMarking} \
      --enable-variant-caller ~{enableVariantCaller} \
      --enable-targeted ~{enableTargeted} \
      --enable-cnv ~{enableCnv} \
      --cnv-enable-self-normalization ~{enableCnvSelfNormalization} \
      --enable-sv ~{enableSv} \
      --dbsnp ~{dbSNP} \
      --validate-pangenome-reference=false \
      --output-directory . \
      --output-file-prefix ~{outputFileNamePrefix} ~{additionalParameters}
    >>>
    runtime {
        backend: "DRAGEN" 
        timeout: "~{timeout}"
        memory: "~{jobMemory} GB"
    }
    
    output {
        File fastqcMetrics = "~{outputFileNamePrefix}.fastqc_metrics.csv"
        File mappingMetrics = "~{outputFileNamePrefix}.mapping_metrics.csv"
        File coverageMetrics = "~{outputFileNamePrefix}.wgs_coverage_metrics.csv"
        File? vcVcf = "~{outputFileNamePrefix}.vcf.gz"
        File? vcMetrics = "~{outputFileNamePrefix}.vc_metrics.csv"
        File? hardfilteredVcf = "~{outputFileNamePrefix}.hard-filtered.vcf.gz"
        File? targetedVcf = "~{outputFileNamePrefix}.targeted.vcf.gz"
        File? ploidyVcf = "~{outputFileNamePrefix}.ploidy.vcf.gz"
        File? cnvVcf = "~{outputFileNamePrefix}.cnv.vcf.gz"
        File? cnvMetrics = "~{outputFileNamePrefix}.cnv_metrics.csv"
        File? svVcf = "~{outputFileNamePrefix}.sv.vcf.gz"
        File? svMetrics = "~{outputFileNamePrefix}.sv_metrics.csv"
    }

    meta {
        output_meta: {
            fastqcMetrics: "fastqc metrics, required",
            mappingMetrics: "mapping metrics, required",
            coverageMetrics: "coverage metrics, required",
            ploidyVcf: "Ploidy vcf",
            vcVcf: "unfiltered vcf from variant calling, requires enableVariantCaller=true",
            vcMetrics: "metrics from variant calling, requires enableVariantCaller=true",
            hardfilteredVcf: "hard-filtered vcf file with variants with filter info attached, requires enableVariantCaller=true",
            targetedVcf: "vcf file from targeted, requires enableTargeted=true",
            cnvVcf: "vcf file with copy number variation, requires enableCnv=true",
            cnvMetrics: "metrics from copy number variation, requires enableCnv=true",
            svVcf: "vcf file with structural variants, requires enableSv=true",
            svMetrics: "metrics from structual variation, requires enableSv=true"
        }
    }

}


