
# Read command line arguments through optparse package ####

suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "Input file path relative to working directory,
              excluding read and file extension (e.g. /path/to/sample1 instead of /path/to/sample1_R1.fastq",
              metavar = "character"),
  make_option(c("-d", "--work_dir"), type = "character", default = NULL,
              help = "Path to working directory", metavar = "character"),
  make_option(c("--rdna_index"), type = "character",
              default = "/home/cluster/o_heiningc/data/genome_indices/mouse_rDNA/bowtie2_index/mm_rdna_bt2",
              help = "Path to rDNA index to remove rDNA from reads [default %default]", metavar = "character"),
  make_option(c("-t", "--threads"), type = "interger", default = NULL,
              help = "Number of threads for pipeline (SLURM!)", metavar = "integer"),
  make_option(c("--adapter_r1"), type = "character", default = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
              help = "Sequence for read1 adapter to be trimmed [default %default]", metavar = "character"),
  make_option(c("--adapter_r2"), type = "character", default = "GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT",
              help = "Sequence for read2 adapter to be trimmed [default %default]", metavar = "character"),
  make_option(c("--genome_index"), type = "character",
              default = "home/gpfs/o_heiningc/data/genome_indices/mm39/STAR_index/",
              help = "Path to genome index for alignment [default %default]", metavar = "character"),
  make_option(c("--chrom_sizes"), type = "character",
              default = "home/cluster/o_heiningc/data/genome_indices/GRCm39/STAR_index/chrNameLength.txt",
              help = "Path to file containing chromosome sizes [default %default]", metavar = "character"),
  make_option(c("-v", "--verbose"), action = "store", default = TRUE,
              help = "Have the program print information [default %default]"),
  make_option(c("-q", "--quiet"), action = "store_false", dest = "verbose",
              help = "Turn off verbose")
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


if(opt$verbose) {
  
  # Print the input file
  cat("Input_file:\n")
  cat(opt$input_file)
  
  # Print the working directory path
  cat("\n\nWorking directory:\n")
  cat(opt$work_dir)
  
  # Print rDNA index path
  cat("\n\nPath to rDNA index:\n")
  cat(opt$rdna_index)
  
  # Print number of threads for pipeline
  cat("\n\nNumber of threads:\n")
  cat(opt$threads)
  
  # Print read1 adapter to be trimmed
  cat("\n\nRead1 adapter to trim:\n")
  cat(opt$adapter_r1)
  
  # Print read2 adapter to be trimmed
  cat("\n\nRead2 adapter to trim:\n")
  cat(opt$adapter_r2)
  
  # Print genome index path
  cat("\n\nPath to genome index:\n")
  cat(opt$genome_index)
  
  # Print chromosome sizes path
  cat("\n\nPath to chrom.sizes:\n")
  cat(opt$chrom_sizes)
  
}


# Set basic parameters ####

work_dir <- opt$work_dir

logs_dir <- paste0(work_dir, "/analysis/logs")
if(!file.exists(logs_dir)) { dir.create(logs_dir, recursive = TRUE) }

read_1 <- paste0(opt$input_file, "_R1.fastq")
read_2 <- paste0(opt$input_file, "_R2.fastq")


threads <- as.integer(opt$threads)


setwd(work_dir)

# Define execute function ####

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


# Step 1: Removing reads without UMIs with cutadapt ####

cutadapt_out <- paste0(work_dir, "/data/UMIonly")
if(!file.exists(cutadapt_out)) { dir.create(cutadapt_out, recursive = TRUE) }

cutadapt_log <- paste0(logs_dir, "/cutadapt")
if(!file.exists(cutadapt_log)) { dir.create(cutadapt_log, recursive = TRUE) }

cutadapt_threads <- as.integer(threads / 3)


cutadapt_cmd <- paste0("cutadapt -g ^\"N{6}A\" -G ^\"N{6}C\" --discard-untrimmed --action=retain --cores=",
                       cutadapt_threads, " --pair-filter=both ", read_1, " ", read_2,
                       " -o ", cutadapt_out, "/", basename(read_1), " -p ", cutadapt_out, "/", basename(read_2),
                       " 2>&1 ", cutadapt_log, "/", basename(opt$input_file), "_cutadapt.log")

execute(cutadapt_cmd, outputfile = paste0(cutadapt_out, "/", basename(read_1)))


# Step 2: Processing reads with fastp and bowtie2 ####

processing_out <- paste0(work_dir, "/data/processed_fastq")
if(!file.exists(processing_out)) { dir.create(processing_out, recursive = TRUE) }

fastp_log <- paste0(logs_dir, "/fastp")
if(!file.exists(fastp_log)) { dir.create(fastp_log, recursive = TRUE) }

rrna_log <- paste0(logs_dir, "/rRNA")
if(!file.exists(rrna_log)) { dir.create(rrna_log, recursive = TRUE) }

processing_threads <- as.integer((threads / 3) * 2)


processing_cmd <- paste0("(fastp -i ", cutadapt_out, "/", basename(read_1), " -I ", cutadapt_out, "/", basename(read_2),
                         " --umi --umi_loc=per_read --umi_len=6 --trim_poly_g --trim_poly_x --adapter_sequence ", opt$adapter_r1,
                         " --adapter_sequence_r2 ", opt$adapter_r2, " --stdout --html ", fastp_log, "/", basename(opt$input_file),
                         "_fastp.html --json ", fastp_log, "/", basename(opt$input_file), "_fastp.json --thread ", processing_threads,
                         " -c 2>&1 ", fastp_log, "/", basename(opt$input_file), "_fastp.log) | (bowtie2 --fast-local --un-conc ",
                         processing_out, "/", basename(opt$input_file), ".fastq --interleaved - -x ", opt$rdna_index,
                         " --threads ", processing_threads, " 2>&1 ", rrna_log, "/", basename(opt$input_file), "_rRNA_bowtie.log) > /dev/null")


execute(processing_cmd, outputfile = paste0(processing_out, "/", basename(read_1)))


# Step 2.1: Renaming output from bowtie2 ####

old_names <- list.files(path = processing_out)

new_names <- gsub("\\.R", "\\_R", old_names)

file.rename(paste0(processing_out, old_names), paste0(processing_out, new_names))


# Step 3: Aligning reads to mm39 using STAR aligner ####

sorted_out <- paste0(work_dir, "/data/sorted_fastq")
if(!file.exists(sorted_out)) { dir.create(sorted_out, recursive = TRUE) }

align_out <- paste0(work_dir, "/data/BAM")
if(!file.exists(align_out)) { dir.create(align_out, recursive = TRUE) }

sorted_files <- c(paste0(sorted_out, "/", basename(opt$input_file), "_R1_sorted.fastq"),
                  paste0(sorted_out, "/", basename(opt$input_file), "_R2_sorted.fastq"))

sort1_cmd <- paste0("fastq-sort -n -S 4G \"", processing_out, "/", basename(read_1), "\" > ", sorted_files[1])

sort2_cmd <- paste0("fastq-sort -n -S 4G \"", processing_out, "/", basename(read_2), "\" > ", sorted_files[2])

execute(sort1_cmd, sorted_files[1])

execute(sort2_cmd, sorted_files[2])


aligned_file <- c(paste0(align_out, "/", basename(opt$input_file), "_Aligned.out.sam"))


align_cmd <- paste0("STAR --runThreadN ", threads, " --genomeDir ", opt&genome_index, 
                    " --readFilesIn ", sorted_files[1], " ", sorted_files[2], " outFileNamePrefix \"", align_out, 
                    "/", opt$input_file, "_\"")

execute(align_cmd, aligned_file)


# Step 3.1: Sort and index with samtools ####

bam_out <- paste0(align_out, "/", basename(opt$input_file), ".BAM")

samtools_threads <- as.integer((threads / 3) * 2)

samviewsort_cmd <- paste0("samtools view -bS \"", aligned_file[1], "\" | samtools sort -@ ",
                          samtools_threads, " -o \"", bam_out, "\"")

execute(samviewsort_cmd, bam_out)


samindex_cmd <- paste0("samtools index \"", bam_out, "\"")

execute(samindex_cmd, outputfile = paste0(bam_out, ".bai"))


# Step 4: Deduplicating BAM files based on UMIs ####

deduped_out <- paste0(work_dir, "/data/deduped_BAM")
if(!file.exists(deduped_out)) { dir.create(deduped_out, recursive = TRUE) }

deduped_log <- paste0(logs_dir, "/dedup")
if(!file.exists(deduped_log)) { dir.create(deduped_log, recursive = TRUE) }


deduped_file <- paste0(deduped_out, "/", basename(opt$input_file), "_deduped.BAM")


deduped_cmd <- paste0("(umi_tools dedup -I \"", bam_out, "\" --paired --umi-separator=\":\" --output-stats=\"",
                      deduped_log, "/", basename(opt$input_file), "\" -S \"", deduped_file, "\") 2>&1 \"",
                      deduped_log, "/", basename(opt$input_file), "_deduped.log && samtools index \"", deduped_file, "\"")

execute(deduped_cmd, deduped_file)


# Step 5: Generating bigwig files for visualization ####

bigwig_out <- paste0(work_dir, "/data/bw")
if(!file.exists(bigwig_out)) { dir.create(bigwig_out, recursive = TRUE) }

bigwig_log <- paste0(logs_dir, "bigwig")
if(!file.exists(bigwig_log)) { dir.create(bigwig_log, recursive = TRUE) }


fwd_out <- paste0(bigwig_out, "/", basename(opt$input_file), "_fwd.bw")

rev_out <- paste0(bigwig_out, "/", basename(opt$input_file), "_rev.bw")

inv_out <- paste0(bigwig_out, "/", basename(opt$input_file), "_rev_inv.bw")


fwd_cmd <- paste0("bamCoverage --bam \"", deduped_file, "\" --skipNonCoveredRegions --outFileName \"", fwd_out,
                  "\" --binSize 1 --numberOfProcessors ", threads, " --normalizeUsing None --Offset 1 --samFlagInclude 82")

rev_cmd <- paste0("bamCoverage --bam \"", deduped_file, "\" --skipNonCoveredRegions --outFileName \"", rev_out,
                  "\" --binSize 1 --numberOfProcessors ", threads, " --normalizeUsing None --Offset 1 --samFlagInclude 98")


execute(fwd_cmd, fwd_out)

execute(rev_cmd, rev_out)


# Step 5.1: Invert values of rev strand bigwig for visualization ####

tmp_dir <- paste0(work_dir, "/tmp")
if(!file.exists(tmp_dir)) { dir.create(tmp_dir, recursive = TRUE) }


inv_cmd <- paste0("bigWigToBedGraph ", rev_out, " \"", tmp_dir, "/", basename(opt$input_file), ".bedgraph\" && cat \"",
                  tmp_dir, "/", basename(opt$input_file), ".bedgraph\" | awk 'BEGIN{OFS=\"t\"} {print $1,$2,$3,-1*$4}' | sort -S 2G -k1,1, -k2,2n > \"",
                  tmp_dir, "/", basename(opt$input_file), "_inv.bedgraph\" && bedGraphToBigWig \"", tmp_dir, "/", basename(opt$input_file),
                  "_inv.bedgraph\" ", opt$chrom_sizes, " \"", inv_out, "\"")

execute(inv_cmd, inv_out)


