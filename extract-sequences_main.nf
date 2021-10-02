#!/usr/bin/env nextflow

/*================================================================
The Aguilar Lab presents...

- A pipeline to extract and create miRNA and 3'UTR consensus sequences for analysis
   with targetscan and miRmap.

==================================================================
Version: 0.1
Project repository:
==================================================================
Authors:

- Bioinformatics Design
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)



- Bioinformatics Development
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)


- Nextflow Port
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)


=============================
Pipeline Processes In Brief:

Pre-processing:
_pre1_split_chromosomes
_pre2_make_consensus_sequence

Core-processing:
_001A_extract_mirna_FASTA_ref
_002A_extract_utr_FASTA_ref
_001B_extract_mirna_FASTA_mut
_002B_extract_utr_FASTA_mut

Pos-processing:
_01A_merge_mirref_fastas
_02A_merge_utrref_fastas
_01B_merge_mirmut_fastas
_02B_merge_utrmut_fastas

Analysis:

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
	The miRNA and 3'UTR consensus sequences extractor pipeline
  v${version}
  ==========================================

	Usage:

	nextflow run ${pipeline_name}.nf --mirnabed <path to input 1> --utrbed <path to input 2>
  --vcf <path to input 3> --fasta <path to input 4>   [--output_dir path to results ]

	  --mirnabed	<- miRNA bed file;

	  --utrbed	<- UTR bed file;

    --vcf <- VCF file;

    --fasta_dir <- Directory whith the FASTA files;

	  --output_dir     <- directory where results, intermediate and log files will be stored;
	      default: same dir where --query_fasta resides

	  -resume	   <- Use cached results if the executed project has been run before;
	      default: not activated
	      This native NF option checks if anything has changed from a previous pipeline execution.
	      Then, it resumes the run from the last successful stage.
	      i.e. If for some reason your previous run got interrupted,
	      running the -resume option will take it from the last successful pipeline stage
	      instead of starting over
	      Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "nf-extract-sequences_main.nf"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.mirnabed = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.utrbed = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.vcf = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.fasta_dir = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "${pipeline_name} v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at MAY 2021
*/
nextflow_required_version = '20.01.0'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
*/

/* Check if the input directory is provided
    if it was not provided, it keeps the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
if ( !params.mirnabed | !params.utrbed | !params.vcf | !params.fasta_dir ) {
  log.error " Please provide the --mirnabed AND --utrbed AND --vcf --fasta \n\n" +
  " For more information, execute: nextflow run extract-sequences.nf --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = file(params.mirnabed).getParent() //!! maybe creates bug, should check

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable pipeline_name defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*
Useful functions definition
*/

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The nf-miRNA-compare
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['Input miRNA bed']			= params.mirnabed
pipelinesummary['Input 3UTR bed']			= params.utrbed
pipelinesummary['Input VCF']			= params.vcf
pipelinesummary['Input FASTA Dir']			= params.fasta_dir
pipelinesummary['Results Dir']		= results_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/

/* enable DSL2*/
nextflow.enable.dsl=2

/*
	READ GENERAL INPUTS
*/

/* _pre1_split_chromosomes */
/* Load VCF file into channel */
Channel
	.fromPath( "${params.vcf}" )
	.set{ vcf_input }

/* Process _pre2_make_consensus_sequence
/* Load FASTA files into channel */
Channel
	.fromPath( "${params.fasta_dir}*.fa" )
	.toList()
	.set{ fasta_input }

/* Process _001A_extract_mirna_FASTA_ref */
/* Load mirna bed file into channel */
Channel
.fromPath( "${params.mirnabed}" )
.set{ mirnabed_input}

/*_002A_extract_utr_FASTA_ref
/* Load utr bed file into channel */
Channel
.fromPath( "${params.utrbed}" )
.set{ utrbed_input}


/*
	READ mk-files
*/

/* _pre1_split_chromosomes */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-split-chromosomes/*")
	.toList()
	.set{ mk_files }

/* Process _pre2_make_consensus_sequence
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/make-consensus-sequence/*")
	.toList()
	.set { mkfiles_pre2 }

/* Process _001A_extract_mirna_FASTA_ref */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-extract-mirna-FASTA-reference/*")
	.toList()
	.set{ mkfiles_001A }

/*_001B_extract_mirna_FASTA_mut
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-extract-mirna-FASTA-consensus/*")
	.toList()
	.set{ mkfiles_001B }

/*_002A_extract_utr_FASTA_ref
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-extract-utr-FASTA-reference/*")
	.toList()
	.set{ mkfiles_002A }

/*_002B_extract_utr_FASTA_mut
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-extract-utr-FASTA-consensus/*")
	.toList()
	.set{ mkfiles_002B }


/*	 * Import modules */
include {  _pre1_split_chromosomes;
					_pre2_make_consensus_sequence;
					_001A_extract_mirna_FASTA_ref;
					_001B_extract_mirna_FASTA_mut;
					_002A_extract_utr_FASTA_ref;
					_002B_extract_utr_FASTA_mut;
					_01A_merge_mirref_fastas;
					_02A_merge_utrref_fastas;
					_01B_merge_mirmut_fastas;
					_02B_merge_utrmut_fastas
					} from './extract-sequences_modules.nf'

/*  main pipeline logic */
workflow  {
						// PRE 1: Split VCF file
						_pre1_split_chromosomes(
									vcf_input,
									mk_files)
									
						// PŔE 2: Make consensus FASTA sequence
						_pre2_make_consensus_sequence(
										_pre1_split_chromosomes.out,
										fasta_input,
										mkfiles_pre2)
										
						// COre 1A: Extract miRNA FASTA reference sequences
						_001A_extract_mirna_FASTA_ref(
										fasta_input,
										mirnabed_input,
										mkfiles_001A)
										
						// Core 1B: Extract miRNA FASTA mutate sequences
						_001B_extract_mirna_FASTA_mut(
										_pre2_make_consensus_sequence.out,
										mirnabed_input,
										mkfiles_001B)
										
						// Core 2A: Extract UTR FASTA reference sequences
						_002A_extract_utr_FASTA_ref(
										fasta_input,
										utrbed_input,
										mkfiles_002A)
										
						// Core 2B: Extract UTR FASTA mutate sequences
						_002B_extract_utr_FASTA_mut(
										_pre2_make_consensus_sequence.out,
										utrbed_input,
										mkfiles_002B)
										
						// Pos 1A: Merge miRNA reference utr_FASTAs
						_01A_merge_mirref_fastas(_001A_extract_mirna_FASTA_ref.out)
						
						// Pos 2A: Merge UTR reference FASTAs
						_02A_merge_utrref_fastas(_002A_extract_utr_FASTA_ref.out)
						
						// Pos !B: Merge miRNA mutate FASTAs
						_01B_merge_mirmut_fastas(_001B_extract_mirna_FASTA_mut.out)
						
						// POs 2B: Merge UTR mutate FASTAs
						_02B_merge_utrmut_fastas(_002B_extract_utr_FASTA_mut.out)

}
