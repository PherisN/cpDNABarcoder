make_test_extract_genome <- function(){
  
  metadata <- list(
    accession = "TEST001",
    version = "TEST001.1",
    organism = "Test organism",
    definition = "Test chloroplast genome",
    genome_length = 2000L
  )
  
  sequence <- paste(
    rep(
      "ACGT",
      500
    ),
    collapse = ""
  )
  
  feature_table <- data.frame(
    type = c(
      "gene",
      "gene",
      "gene",
      "gene"
    ),
    start = c(
      101L,
      301L,
      501L,
      601L
    ),
    end = c(
      120L,
      330L,
      520L,
      610L
    ),
    strand = c(
      STRAND_FORWARD,
      STRAND_REVERSE,
      STRAND_FORWARD,
      STRAND_FORWARD
    ),
    location_type = c(
      LOCATION_SIMPLE,
      LOCATION_SIMPLE,
      LOCATION_SIMPLE,
      LOCATION_SIMPLE
    ),
    location = c(
      "101..120",
      "complement(301..330)",
      "501..520",
      "601..610"
    ),
    gene = c(
      "rbcL",
      "matK",
      "psbA",
      "trnH-GUG"
    ),
    product = c(
      "ribulose-bisphosphate carboxylase large chain",
      "maturase K",
      "photosystem II protein D1",
      "tRNA-His"
    ),
    locus_tag = c(
      "LOC001",
      "LOC002",
      "LOC003",
      "LOC004"
    ),
    protein_id = c(
      "PROT001",
      "PROT002",
      "PROT003",
      NA_character_
    ),
    stringsAsFactors = FALSE
  )
  
  new_cpGenome(
    metadata = metadata,
    sequence = sequence,
    feature_table = feature_table
  )
  
}


testthat::test_that(
  "compute_interval computes a non-wrapping interval between features",
  {
    
    genome <- make_test_extract_genome()
    
    feature_table <- get_feature_table(
      genome
    )
    
    interval <- compute_interval(
      as.list(
        feature_table[3, ]
      ),
      as.list(
        feature_table[4, ]
      )
    )
    
    testthat::expect_true(
      inherits(
        interval,
        "Interval"
      )
    )
    
    testthat::expect_equal(
      get_interval_start(interval),
      521L
    )
    
    testthat::expect_equal(
      get_interval_end(interval),
      600L
    )
    
    testthat::expect_equal(
      get_interval_width(interval),
      80L
    )
    
    testthat::expect_false(
      interval_wraps_origin(interval)
    )
    
  }
)


testthat::test_that(
  "compute_interval orders features by genomic position",
  {
    
    feature1 <- list(
      start = 100L,
      end = 200L,
      gene = "geneA"
    )
    
    feature2 <- list(
      start = 50L,
      end = 80L,
      gene = "geneB"
    )
    
    interval <- compute_interval(
      feature1,
      feature2
    )
    
    testthat::expect_false(
      interval_wraps_origin(interval)
    )
    
    testthat::expect_equal(
      get_interval_start(interval),
      81L
    )
    
    testthat::expect_equal(
      get_interval_end(interval),
      99L
    )
    
    testthat::expect_equal(
      get_interval_width(interval),
      19L
    )
    
  }
)


testthat::test_that(
  "compute_interval rejects invalid features",
  {
    
    feature1 <- list(
      start = 100L,
      end = 200L,
      gene = "geneA"
    )
    
    feature2 <- list(
      start = NA_integer_,
      end = 300L,
      gene = "geneB"
    )
    
    testthat::expect_error(
      compute_interval(
        feature1,
        feature2
      )
    )
    
    testthat::expect_error(
      compute_interval(
        list(
          start = 100L,
          end = 200L
        ),
        feature1
      )
    )
    
    testthat::expect_error(
      compute_interval(
        "not a feature",
        feature1
      )
    )
    
  }
)


testthat::test_that(
  "extract_subsequence extracts non-wrapping intervals",
  {
    
    sequence <- "AACCGGTT"
    
    interval <- new_interval(
      start = 2L,
      end = 5L,
      wraps_origin = FALSE,
      left_gene = "left",
      right_gene = "right"
    )
    
    dna <- extract_subsequence(
      sequence,
      interval
    )
    
    testthat::expect_equal(
      dna,
      "ACCG"
    )
    
  }
)


testthat::test_that(
  "extract_subsequence extracts origin-wrapping intervals",
  {
    
    sequence <- "AACCGGTT"
    
    interval <- new_interval(
      start = 7L,
      end = 2L,
      wraps_origin = TRUE,
      left_gene = "left",
      right_gene = "right"
    )
    
    dna <- extract_subsequence(
      sequence,
      interval
    )
    
    testthat::expect_equal(
      dna,
      "TTAA"
    )
    
  }
)


testthat::test_that(
  "extract_subsequence rejects invalid inputs",
  {
    
    interval <- new_interval(
      start = 1L,
      end = 4L,
      wraps_origin = FALSE,
      left_gene = "left",
      right_gene = "right"
    )
    
    testthat::expect_error(
      extract_subsequence(
        "ATGX",
        interval
      )
    )
    
    testthat::expect_error(
      extract_subsequence(
        "ATGC",
        list()
      )
    )
    
  }
)


testthat::test_that(
  "extract_gene extracts a gene marker from a cpGenome object",
  {
    
    genome <- make_test_extract_genome()
    
    marker_info <- MARKER_DICTIONARY$RBCL
    
    barcode <- extract_gene(
      genome,
      marker_info,
      PIPELINE_OPTIONS
    )
    
    expected_sequence <- substr(
      get_sequence(genome),
      101L,
      120L
    )
    
    testthat::expect_true(
      inherits(
        barcode,
        "cpBarcode"
      )
    )
    
    testthat::expect_equal(
      get_marker(barcode),
      marker_info$marker
    )
    
    testthat::expect_equal(
      get_sequence(barcode),
      expected_sequence
    )
    
    testthat::expect_equal(
      get_coordinates(barcode),
      list(
        start = 101L,
        end = 120L
      )
    )
    
    testthat::expect_equal(
      get_strand(barcode),
      STRAND_FORWARD
    )
    
  }
)


testthat::test_that(
  "extract_gene reverse-complements reverse-strand genes when requested",
  {
    
    genome <- make_test_extract_genome()
    
    marker_info <- MARKER_DICTIONARY$MATK
    
    options <- PIPELINE_OPTIONS
    options$reverse_complement <- TRUE
    
    barcode <- extract_gene(
      genome,
      marker_info,
      options
    )
    
    expected_sequence <- reverse_complement(
      substr(
        get_sequence(genome),
        301L,
        330L
      )
    )
    
    testthat::expect_equal(
      get_sequence(barcode),
      expected_sequence
    )
    
    testthat::expect_equal(
      get_strand(barcode),
      STRAND_REVERSE
    )
    
  }
)


testthat::test_that(
  "extract_intergenic extracts an intergenic marker from a cpGenome object",
  {
    
    genome <- make_test_extract_genome()
    
    marker_info <- MARKER_DICTIONARY$PSBA_TRNH
    
    barcode <- extract_intergenic(
      genome,
      marker_info,
      PIPELINE_OPTIONS
    )
    
    expected_sequence <- substr(
      get_sequence(genome),
      521L,
      600L
    )
    
    testthat::expect_true(
      inherits(
        barcode,
        "cpBarcode"
      )
    )
    
    testthat::expect_equal(
      get_marker(barcode),
      marker_info$marker
    )
    
    testthat::expect_equal(
      get_feature_type(barcode),
      FEATURE_INTERGENIC
    )
    
    testthat::expect_equal(
      get_sequence(barcode),
      expected_sequence
    )
    
    testthat::expect_equal(
      get_coordinates(barcode),
      list(
        start = 521L,
        end = 600L
      )
    )
    
    testthat::expect_false(
      get_metadata(barcode)$wraps_origin
    )
    
  }
)


testthat::test_that(
  "extract_marker dispatches gene and intergenic marker extraction",
  {
    
    genome <- make_test_extract_genome()
    
    rbcL <- extract_marker(
      genome,
      "RBCL"
    )
    
    psbA_trnH <- extract_marker(
      genome,
      "PSBA_TRNH"
    )
    
    testthat::expect_true(
      inherits(
        rbcL,
        "cpBarcode"
      )
    )
    
    testthat::expect_true(
      inherits(
        psbA_trnH,
        "cpBarcode"
      )
    )
    
    testthat::expect_equal(
      get_feature_type(psbA_trnH),
      FEATURE_INTERGENIC
    )
    
  }
)


testthat::test_that(
  "extract_marker rejects unknown markers",
  {
    
    genome <- make_test_extract_genome()
    
    testthat::expect_error(
      extract_marker(
        genome,
        "UNKNOWN_MARKER"
      )
    )
    
  }
)


extract_test_file <- test_genbank_file(
  "NC_000932.gb"
)


testthat::test_that(
  "extract_marker extracts real GenBank gene markers",
  {
    
    genome <- read_cpGenome(
      extract_test_file
    )
    
    rbcL <- extract_marker(
      genome,
      "RBCL"
    )
    
    matK <- extract_marker(
      genome,
      "MATK"
    )
    
    testthat::expect_true(
      inherits(
        rbcL,
        "cpBarcode"
      )
    )
    
    testthat::expect_true(
      inherits(
        matK,
        "cpBarcode"
      )
    )
    
    testthat::expect_equal(
      get_marker(rbcL),
      MARKER_DICTIONARY$RBCL$marker
    )
    
    testthat::expect_equal(
      get_marker(matK),
      MARKER_DICTIONARY$MATK$marker
    )
    
    testthat::expect_gt(
      nchar(
        get_sequence(rbcL)
      ),
      1000L
    )
    
    testthat::expect_gt(
      nchar(
        get_sequence(matK)
      ),
      500L
    )
    
  }
)


testthat::test_that(
  "extract_marker extracts a real GenBank intergenic marker",
  {
    
    genome <- read_cpGenome(
      extract_test_file
    )
    
    psbA_trnH <- extract_marker(
      genome,
      "PSBA_TRNH"
    )
    
    testthat::expect_true(
      inherits(
        psbA_trnH,
        "cpBarcode"
      )
    )
    
    testthat::expect_equal(
      get_marker(psbA_trnH),
      MARKER_DICTIONARY$PSBA_TRNH$marker
    )
    
    testthat::expect_equal(
      get_feature_type(psbA_trnH),
      FEATURE_INTERGENIC
    )
    
    testthat::expect_gt(
      nchar(
        get_sequence(psbA_trnH)
      ),
      1L
    )
    
  }
)