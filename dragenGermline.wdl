version 1.0

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroup
}

struct GenomeResources {
    String dbSNP
    String referenceDirectory
    String dragenVersion
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
      "referenceDirectory": "/.mounts/labs/gsiprojects/gsi/Dragen/reference/hg38fa.p12/",
      "dragenVersion": "4.2.4"
      },
    "hg38_noAlt": {
      "dbSNP": "/.mounts/labs/gsiprojects/gsi/Dragen/reference/dbSNP.151/common_all_dbSNP151_hg38p7_noalt_sorted.vcf.gz",
      "referenceDirectory": "/.mounts/labs/gsiprojects/gsi/Dragen/reference/hg38_noAlt-p12/",
      "dragenVersion": "4.2.4"
      }
    }

    String dragen_ref = dragen_resources_by_genome [ reference ].referenceDirectory
    String dragen_dbsnp = dragen_resources_by_genome [ reference ].dbSNP
    String dragen_version = dragen_resources_by_genome [ reference ].dragenVersion

    meta {
        author: "Peter Ruzanov"
        email: "pruzanov@oicr.on.ca"
        description: "A workflow for calling SNVs on fastq inputs in germline mode"
        dependencies: [
        {
          name: "gsi hg38 modules : hg38-dbsnp/151",
          url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
        },
        {
          name: "gsi modules : dragen-scripts/0.3",
          url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
        }]
        output_meta: {
          unfilteredVcf: {
            description: "SNV calls before applying any filters",
            vidarr_label: "unfilteredVcf"
          },
          unfilteredIdx: {
            description: "Index of SNV calls before applying any filters",
            vidarr_label: "unfilteredIdx"
          },
          filteredVcf: {
            description: "SNV calls with filter information attached",
            vidarr_label: "filteredVcf"
          },
          filteredIdx: {
            description: "Index of SNV calls with filter information attached",
            vidarr_label: "filteredIdx"
          },
          targetedVcf: {
            description: "Targeted vcf file",
            vidarr_label: "targetedVcf"
          },
          targetedIdx: {
            description: "Index of Targeted vcf file",
            vidarr_label: "targetedIdx"
          },
          ploidyVcf: {
            description: "Ploidy vcf file",
            vidarr_label: "ploidyVcf"
          },
          ploidyIdx: {
            description: "Index of Ploidy vcf file",
            vidarr_label: "ploidyIdx"
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
        outputFileNamePrefix = outputFileNamePrefix,
        dragenVersion = dragen_version
    } 

    output {
        File unfilteredVcf = runDragenGermline.outputVcf
        File unfilteredIdx = runDragenGermline.outputIdx
        File filteredVcf = runDragenGermline.hardfilteredVcf
        File filteredIdx = runDragenGermline.hardfilteredIdx
        File? ploidyVcf = runDragenGermline.ploidyVcf
        File? ploidyIdx = runDragenGermline.ploidyIdx
        File? targetedVcf = runDragenGermline.targetedVcf
        File? targetedIdx = runDragenGermline.targetedIdx
    }
}

# =====================================================================
# A scripted extraction of info from RG line to dragen-compliant string
# =====================================================================
task extractInfoLine {
   input {
       InputGroup fastqInput
       String parsingScript = "$DRAGEN_SCRIPTS_ROOT/bin/composeList.py"
       Int timeout = 4
       Int jobMemory = 4
       String modules = "dragen-scripts/0.3"
   }

   parameter_meta {
     fastqInput: "InputGroup struct entry with fastq files"
     parsingScript: "Script for parsing inputs into a line"
     timeout: "Timeout for the job"
     jobMemory: "Job allocated RAM"
     modules: "dependency modules"
   }

   command <<<
    python3 ~{parsingScript} -i ~{write_json(fastqInput)}
   >>>

   runtime {
     timeout: "~{timeout}"
     modules: "~{modules}"
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
      String listWritingScript = "$DRAGEN_SCRIPTS_ROOT/bin/writeFile.py"
      String outputFileName
      Int jobMemory = 4
      Int timeout = 4
      String modules = "dragen-scripts/0.3"
   }

   parameter_meta {
     inputLines: "Array of input lines to print"
     listWritingScript: "Script for writing out list of inputs"
     outputFileName: "Name of an output file, list of inputs"
     jobMemory: "Job allocated RAM"
     timeout: "Timeout for the job"
     modules: "dependency modules"
   }

   command<<<
   python3 ~{listWritingScript} -o ~{outputFileName} -l "~{sep=';' inputLines}"
   >>>
   

   runtime {
      timeout: "~{timeout}"
      modules: "~{modules}"
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
        String refDir
        String dragenVersion
        String? additionalParameters
        String outputFileNamePrefix
        String dbSNP
        Int timeout = 96
    }

    parameter_meta {
        sampleFastqList: "List of tumor fastq files, required input"
        enableDupMarking: "Flag for duplicate marking, true by  default"
        enableTargeted: "Flag for enabling calling on targets like HBA, GBA etc. clusters"
        refDir: "The reference genome directoty"
        dragenVersion: "Expected version of dragen software on the DRAGEN node"
        additionalParameters: "Additional dragen parameters"
        dbSNP: "Path to the dbSNP reference file"
        outputFileNamePrefix: "Output file name prefix"
        timeout: "Hours before task timeout"
    }
    
    String resultVcf = "~{outputFileNamePrefix}.vcf.gz"
    String hardfilteredVcfName = "~{outputFileNamePrefix}.hard-filtered.vcf.gz"
    String ploidyVcfName = "~{outputFileNamePrefix}.ploidy.vcf.gz"
    String targetedVcfName = "~{outputFileNamePrefix}.targeted.vcf.gz"

    command <<<
      dragen -f -r ~{refDir} \
      --fastq-list ~{sampleFastqList} \
      --enable-duplicate-marking ~{enableDupMarking} \
      --enable-variant-caller true \
      --enable-targeted ~{enableTargeted} \
      --dbsnp ~{dbSNP} \
      --output-directory . \
      --output-file-prefix ~{outputFileNamePrefix} ~{additionalParameters}
    >>>
    runtime {
        backend: "DRAGEN"
        dragen_version: "~{dragenVersion}"
        timeout: "~{timeout}"
    }
    
    output {
        File outputVcf = "~{resultVcf}"
        File outputIdx = "~{resultVcf}.tbi"
        File hardfilteredVcf = "~{hardfilteredVcfName}"
        File hardfilteredIdx = "~{hardfilteredVcfName}.tbi"
        File? targetedVcf = "~{targetedVcfName}"
        File? targetedIdx = "~{targetedVcfName}.tbi"
        File? ploidyVcf = "~{ploidyVcfName}"
        File? ploidyIdx = "~{ploidyVcfName}.tbi"
   }

    meta {
        output_meta: {
            outputVcf: "output unfiltered vcf with SNV calls",
            outputIdx: "Index of unfiltered vcf with SNV calls",
            hardfilteredVcf: "Hard-filtered vcf file with variants with filter info attached",
            hardfilteredIdx: "Index of hard-filtered vcf file with variants with filter info attached",
            targetedVcf: "Targeted vcf",
            targetedIdx: "Index of targeted vcf",
            ploidyVcf: "Ploidy vcf",
            ploidyIdx: "Index of ploidy vcf"
        }
    }

}


