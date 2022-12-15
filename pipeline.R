
# Read command line arguments through optparse package ####


library(optparse)

option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "Input file path relative to working directory,
              excluding read and file extension (e.g. /path/to/sample1 instead of /path/to/sample1_R1.fastq",
              metavar = "character"),
  make_option(c("-d", "--work_dir"), type = "character", default = NULL,
              help = "Path to working directory", metavar = "character"),
  make_option(c("--rdna_index"), type = "character",
              default = "/home/cluster/o_heiningc/data/genome_indices/mouse_rDNA/bowtie2_index/mm_rdna_bt2",
              help = "Path to rDNA index to remove rDNA from reads", metavar = "character"),
  make_option(c("-t", "--threads"), type = "interger", default = NULL,
              help = "Number of threads for pipeline (SLURM!)", metavar = "integer"),
  make_option(c("--adapter_r1"), type = "character", default = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
              help = "Sequence for read1 adapter to be trimmed", metavar = "character"),
  make_option(c("--adapter_r2"), type = "character", default = "GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT",
              help = "Sequence for read2 adapter to be trimmed", metavar = "character"),
  make_option(c("--genome_index"), type = "character",
              default = "home/gpfs/o_heiningc/data/genome_indices/mm39/STAR_index/",
              help = "Path to genome index for alignment", metavar = "character"),
  make_option(c("--chrom_sizes"), type = "character",
              default = "home/cluster/o_heiningc/data/genome_indices/GRCm39/STAR_index/chrNameLength.txt",
              help = "Path to file containing chromosome sizes", metavar = "character")
)


opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);


if(is.null(opt$input_file)) {
  print_help(opt_parser)
  stop("Please supply path to input file relative to working directory.n")
}

if(is.null(opt$work_dir)) {
  print_help(opt_parser)
  stop("Please supply path to working directory.n")
}

if(is.null(opt$threads)) {
  print_help(opt_parser)
  stop("Please supply number of threads to be used for pipeline.n")
}


# Set basic parameters ####

work_dir <- opt$work_dir

threads <- registerDoParallel(cores = (as.integer(opt$threads)))

read_1 <- paste0(opt$input_file, "_R1.fastq")
read_2 <- paste0(opt$input_file, "_R2.fastq")


setwd(opt$work_dir)

# Load package for execute function or define if package is not installed ####


if(require("execute")) {
  library(execute)
} else {

  execute <- function(x, outputfile = NA, intern = FALSE, quitOnError = FALSE) {
    if(!is.na(outputfile) && file.exists(outputfile)) {
      cat("Output for step exists, skipping this step\n");
      return("")
    }
    cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
    if(res[1] >= 1) {
      cat("Error external process did not finish \n\n");
      if(quitOnError) q("no")
    }
  }
}


# Step 1: Removing reads without UMIs with cutadapt ####

cutadapt_out <- paste0(work_dir, "/data/UMIonly")
if(!file.exists(cutadapt_out)) { dir.create(cutadapt_out, recursive = TRUE) }

cutadapt_log <- paste0(work_dir, "/analysis/logs/cutadapt")
if(!file.exists(cutadapt_log)) { dir.create(cutadapt_log, recursive = TRUE) }

cutadapt_threads <- as.integer(threads / 3)


cutadapt_cmd <- paste0("cutadapt -g ^\"N{6}A\" -G ^\"N{6}C\" --discard-untrimmed --action=retain --cores=",
                       cutadapt_threads, " --pair-filter=both ", read_1, " ", read_2,
                       " -o ", cutadapt_out, "/", basename(read_1),
                       " -p ", cutadapt_out, "/", basename(read_2),
                       " 2>&1 ", cutadapt_log, "/", basename(opt$input_file), "_cutadapt.log")


# cutadapt \
# -g ^"N{6}A" \
# -G ^"N{6}C" \
# --discard-untrimmed \
# --action=retain \
# --cores=$(echo ${THREADS}/3 | bc) \
# --pair-filter=both ${PAIR}_R1.fastq ${PAIR}_R2.fastq \
# -o data/UMIonly/$(basename ${PAIR})_R1.fastq \
# -p data/UMIonly/$(basename ${PAIR})_R2.fastq 2>&1 analysis/logs/cutadapt/$(basename ${PAIR})_cutadapt.log





















