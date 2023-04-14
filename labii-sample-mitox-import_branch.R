# labii-sample-import-mitox.R
#
# 2020-15-20
#
# This script outputs two tsv files, one for patients, one for samples, that can be imported into Labii.
# The inputs are the RecordID values. Five samples are created for each RecordID
#
# opt$record_id <- 4002
# opt$next_box_loc <- "C3"

if (!require("optparse")) {
  install.packages("optparse")
  library(optparse)
}

option_list <- list(
  make_option(c("-r", "--record_id"), type = "character", default = NULL,
              help='record ids, format is id1,id2,id3', metavar = "character"),
  make_option(c("-o", "--out_patients"), type = "character", default = "out-patients.tsv",
              help="output file name for the patients [default= %default]",
              metavar="character"),
  make_option(c("-s", "--out_samples"), type = "character", default="out-samples.tsv",
              help="output file name for the samples [default= %default]",
              metavar="character"),
  make_option(c("-n", "--next_box_loc"), type = "character", default = NULL,
              help = "next free location in the freezer box",
              metavar = "character"),
  make_option(c("-b", "--freezer_box"), type = "integer", default = "1",
              help = "freezer box number",
              metavar = "character"),
  make_option(c("-p", "--project"), type = "character", default = "mitox",
              help = "freezer box number",
              metavar = "character"),
  make_option(c("-d", "--date"), type = "character", default = Sys.Date(),
              help = "date samples were processed",
              metavar = "character"),
  make_option(c("-t", "--time"), type = "character", default = format(Sys.time(), "%H:%M"),
              help = "time samples were processed",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$record_id)){
  print_help(opt_parser)
  stop("At least one record id must be supplied", call.=FALSE)
}

# Convert record_id
try({
  opt$record_id <- as.character(
    unlist(
      strsplit(opt$record_id, split = ",")
    )
  )
})

# ### Patient table ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# ###
#
# p.fields <- c("name", "record_id", "project")
#
# pts <- matrix(nrow = length(opt$record_id), ncol = length(p.fields))
#
# pts <- as.data.frame(pts)
#
# names(pts) <- p.fields
#
# pts$record_id <- opt$record_id
#
# pts$name <- pts$record_id
#
# write.table(pts,
#             file = paste(opt$date, opt$out_patients, sep = "_"),
#             sep = "\t", row.names = FALSE, quote = FALSE)

### Sample table ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###

s.fields <- c("name", "patient", "record_id", "freezer_box", "box_position", "date_processed",
              "time_processed", "project", "mass")

samp <- matrix(nrow = length(opt$record_id) * 2, ncol = length(s.fields))

samp <- as.data.frame(samp)

names(samp) <- s.fields

samp$record_id <- rep(opt$record_id, each = 2)
samp$patient <- samp$record_id

samp$freezer_box <- opt$freezer_box

# box position
#
box_map <- paste(rep(LETTERS[1:9], each = 9), rep(seq(1:9), 9), sep = "")
start_loc <- grep(opt$next_box_loc, box_map)
samp$box_position <- box_map[start_loc:(start_loc + (nrow(samp) - 1))]

# Date
samp$date_processed <- opt$date
samp$date_processed <- as.Date(samp$date_processed, format = "%F")

# Time
samp$time_processed <- opt$time

# Project
samp$project <- opt$project

# Sample ID
#
samp$sample_id <- paste(samp$record_id, format(samp$date_processed, "%y%m%d"), seq(1:2), sep = ".")

# Name (label in Labii)
samp$name <- paste(samp$project, samp$record_id, paste(samp$freezer_box,
                                                       samp$box_position, sep = "."), samp$date_processed, sep = " ")
# Write output
write.table(samp,
            file = paste(as.character(samp$date_processed[1]), opt$out_samples, sep = "_"),
            sep = "\t", row.names = FALSE, quote = FALSE)
