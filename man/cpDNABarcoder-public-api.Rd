
\name{cpDNABarcoder-public-api}
\alias{cpDNABarcoder-public-api}
\alias{read_cpGenome}
\alias{extract_marker}
\alias{extract_barcode}
\alias{extract_barcodes}
\alias{batch_extract_barcode}
\alias{batch_extract_barcodes}
\alias{qc_barcode}
\alias{summarize_barcodes}
\alias{report_qc}
\alias{report_summary}
\alias{as_fasta}
\alias{as_dataframe}
\alias{write_fasta}
\alias{write_csv}
\alias{write_barcode}
\alias{write_report}
\alias{get_metadata}
\alias{get_sequence}
\alias{get_feature_table}
\alias{get_accession}
\alias{get_version}
\alias{get_definition}
\alias{get_organism}
\alias{get_genome_length}
\alias{get_marker}
\alias{get_coordinates}
\alias{get_strand}
\alias{get_feature_type}

\title{Public API for cpDNABarcoder}

\description{
Public user-facing functions for reading chloroplast GenBank genomes,
extracting DNA barcode markers, computing quality-control metrics,
summarizing results, exporting outputs, and accessing object fields.
}

\usage{
read_cpGenome(file)

extract_marker(genome, marker, options = PIPELINE_OPTIONS)

extract_barcode(file, marker, output = NULL, format = "fasta",
  options = PIPELINE_OPTIONS)

extract_barcodes(file, markers, output_dir = NULL, format = "fasta",
  options = PIPELINE_OPTIONS)

batch_extract_barcode(files, marker, output_dir = NULL, format = "fasta",
  options = PIPELINE_OPTIONS)

batch_extract_barcodes(files, markers, output_dir = NULL, format = "fasta",
  options = PIPELINE_OPTIONS)

qc_barcode(barcode)

summarize_barcodes(barcodes)

report_qc(qc)

report_summary(summary)

as_fasta(barcode, width = DEFAULT_FASTA_WIDTH)

as_dataframe(barcode)

write_fasta(barcode, file, width = 80)

write_csv(barcode, file)

write_barcode(barcode, file, format = "fasta", ...)

write_report(report, file, format = NULL)

get_metadata(x)

get_sequence(x)

get_feature_table(x)

get_accession(x)

get_version(x)

get_definition(x)

get_organism(x)

get_genome_length(x)

get_marker(x)

get_coordinates(x)

get_strand(x)

get_feature_type(x)
}

\arguments{
\item{file}{Path to an input GenBank file or output file, depending on the function.}

\item{files}{Character vector of GenBank file paths.}

\item{genome}{A \code{cpGenome} object.}

\item{barcode}{A \code{cpBarcode} object.}

\item{barcodes}{A list of \code{cpBarcode} objects.}

\item{marker}{A single marker name, such as \code{"RBCL"}, \code{"MATK"}, or \code{"PSBA_TRNH"}.}

\item{markers}{Character vector of marker names.}

\item{output}{Optional output file path.}

\item{output_dir}{Optional output directory.}

\item{format}{Output format. Supported values include \code{"fasta"}, \code{"csv"}, \code{"txt"}, \code{"md"}, and \code{"html"}, depending on the function.}

\item{options}{Pipeline options controlling extraction behavior.}

\item{qc}{A \code{cpBarcodeQC} object.}

\item{summary}{A barcode summary object.}

\item{report}{A \code{cpBarcodeReport} object.}

\item{width}{Maximum FASTA sequence line width.}

\item{x}{An object with an appropriate S3 accessor method.}

\item{...}{Additional arguments passed to lower-level exporters.}
}

\value{
Depending on the function, returns a \code{cpGenome}, \code{cpBarcode},
\code{cpBarcodeQC}, summary object, report object, character vector,
data frame, named list, accessor value, or invisibly returns the written
output file path.
}

\details{
The package currently focuses on chloroplast GenBank genome records and
plastid barcode markers. Supported marker workflows include gene markers
such as \code{RBCL} and \code{MATK}, and intergenic markers such as
\code{PSBA_TRNH}.
}

\seealso{
\code{\link{extract_barcode}},
\code{\link{extract_barcodes}},
\code{\link{batch_extract_barcode}},
\code{\link{qc_barcode}},
\code{\link{write_barcode}}
}

\keyword{package}

