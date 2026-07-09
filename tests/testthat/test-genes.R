make_test_feature_table <- function(){
  
  data.frame(
    type = c(
      "gene",
      "CDS",
      "gene",
      "gene",
      "gene",
      "gene"
    ),
    gene = c(
      "rbcL",
      "rbcL",
      "matK",
      "trnH-GUG",
      NA_character_,
      "psbA"
    ),
    start = c(
      100L,
      100L,
      500L,
      900L,
      1200L,
      1500L
    ),
    end = c(
      400L,
      400L,
      800L,
      1000L,
      1300L,
      1800L
    ),
    strand = c(
      STRAND_FORWARD,
      STRAND_FORWARD,
      STRAND_REVERSE,
      STRAND_FORWARD,
      STRAND_FORWARD,
      STRAND_REVERSE
    ),
    stringsAsFactors = FALSE
  )
  
}


testthat::test_that(
  "standardize_gene_name standardizes gene names",
  {
    
    testthat::expect_equal(
      standardize_gene_name(" rbcL "),
      "RBCL"
    )
    
    testthat::expect_equal(
      standardize_gene_name("trnH_gug"),
      "TRNH-GUG"
    )
    
    testthat::expect_error(
      standardize_gene_name(
        c("rbcL", "matK")
      )
    )
    
    testthat::expect_error(
      standardize_gene_name(
        NA_character_
      )
    )
    
  }
)


testthat::test_that(
  "find_gene locates matching gene annotations",
  {
    
    feature_table <- make_test_feature_table()
    
    matches <- find_gene(
      feature_table,
      "rbcL"
    )
    
    testthat::expect_true(
      is.data.frame(matches)
    )
    
    testthat::expect_equal(
      nrow(matches),
      2L
    )
    
    testthat::expect_true(
      all(
        matches$gene == "rbcL"
      )
    )
    
  }
)


testthat::test_that(
  "find_gene recognises gene aliases",
  {
    
    feature_table <- make_test_feature_table()
    
    matches <- find_gene(
      feature_table,
      "trnH"
    )
    
    testthat::expect_equal(
      nrow(matches),
      1L
    )
    
    testthat::expect_equal(
      matches$gene,
      "trnH-GUG"
    )
    
  }
)


testthat::test_that(
  "find_gene rejects invalid inputs",
  {
    
    feature_table <- make_test_feature_table()
    
    testthat::expect_error(
      find_gene(
        list(),
        "rbcL"
      )
    )
    
    testthat::expect_error(
      find_gene(
        feature_table[
          ,
          setdiff(
            names(feature_table),
            "gene"
          )
        ],
        "rbcL"
      )
    )
    
    testthat::expect_error(
      find_gene(
        feature_table,
        "unknownGene"
      )
    )
    
  }
)


testthat::test_that(
  "find_gene warns when no annotation is found",
  {
    
    feature_table <- make_test_feature_table()
    
    feature_table$gene[
      feature_table$gene == "matK"
    ] <- "other"
    
    testthat::expect_warning(
      matches <- find_gene(
        feature_table,
        "matK"
      )
    )
    
    testthat::expect_equal(
      nrow(matches),
      0L
    )
    
  }
)


testthat::test_that(
  "select_gene_feature selects the preferred feature type",
  {
    
    feature_table <- make_test_feature_table()
    
    matches <- find_gene(
      feature_table,
      "rbcL"
    )
    
    selected <- select_gene_feature(
      matches,
      "rbcL"
    )
    
    testthat::expect_true(
      is.data.frame(selected)
    )
    
    testthat::expect_equal(
      nrow(selected),
      1L
    )
    
    testthat::expect_equal(
      selected$type,
      "gene"
    )
    
  }
)


testthat::test_that(
  "select_gene_feature rejects unknown genes",
  {
    
    feature_table <- make_test_feature_table()
    
    testthat::expect_error(
      select_gene_feature(
        feature_table,
        "unknownGene"
      )
    )
    
  }
)


testthat::test_that(
  "validate_feature accepts a valid feature",
  {
    
    feature_table <- make_test_feature_table()
    
    feature <- feature_table[
      !is.na(feature_table$gene) &
        feature_table$gene == "matK",
      ,
      drop = FALSE
    ]
    
    validated <- validate_feature(
      feature,
      "MATK"
    )
    
    testthat::expect_true(
      is.data.frame(validated)
    )
    
    testthat::expect_equal(
      nrow(validated),
      1L
    )
    
  }
)


testthat::test_that(
  "validate_feature rejects invalid features",
  {
    
    feature_table <- make_test_feature_table()
    
    empty_feature <- feature_table[
      FALSE,
      ,
      drop = FALSE
    ]
    
    testthat::expect_error(
      validate_feature(
        empty_feature,
        "RBCL"
      )
    )
    
    invalid_coordinates <- feature_table[
      1,
      ,
      drop = FALSE
    ]
    
    invalid_coordinates$start <- NA_integer_
    
    testthat::expect_error(
      validate_feature(
        invalid_coordinates,
        "RBCL"
      )
    )
    
    invalid_strand <- feature_table[
      1,
      ,
      drop = FALSE
    ]
    
    invalid_strand$strand <- "?"
    
    testthat::expect_error(
      validate_feature(
        invalid_strand,
        "RBCL"
      )
    )
    
    missing_type <- feature_table[
      1,
      ,
      drop = FALSE
    ]
    
    missing_type$type <- NA_character_
    
    testthat::expect_error(
      validate_feature(
        missing_type,
        "RBCL"
      )
    )
    
  }
)


testthat::test_that(
  "validate_feature warns and keeps the first record when multiple annotations are supplied",
  {
    
    feature_table <- make_test_feature_table()
    
    duplicate <- rbind(
      feature_table[1, , drop = FALSE],
      feature_table[1, , drop = FALSE]
    )
    
    testthat::expect_warning(
      validated <- validate_feature(
        duplicate,
        "RBCL"
      )
    )
    
    testthat::expect_equal(
      nrow(validated),
      1L
    )
    
  }
)


testthat::test_that(
  "get_gene_feature retrieves a validated gene feature",
  {
    
    feature_table <- make_test_feature_table()
    
    feature <- get_gene_feature(
      feature_table,
      "psbA"
    )
    
    testthat::expect_true(
      is.data.frame(feature)
    )
    
    testthat::expect_equal(
      nrow(feature),
      1L
    )
    
    testthat::expect_equal(
      feature$gene,
      "psbA"
    )
    
    testthat::expect_equal(
      feature$strand,
      STRAND_REVERSE
    )
    
  }
)


genes_test_file <- test_genbank_file(
  "NC_000932.gb"
)


testthat::test_that(
  "get_gene_feature retrieves real GenBank gene annotations",
  {
    
    gb_lines <- read_genbank(
      genes_test_file
    )
    
    feature_table <- parse_features(
      gb_lines
    )
    
    rbcL <- get_gene_feature(
      feature_table,
      "rbcL"
    )
    
    matK <- get_gene_feature(
      feature_table,
      "matK"
    )
    
    psbA <- get_gene_feature(
      feature_table,
      "psbA"
    )
    
    testthat::expect_equal(
      nrow(rbcL),
      1L
    )
    
    testthat::expect_equal(
      nrow(matK),
      1L
    )
    
    testthat::expect_equal(
      nrow(psbA),
      1L
    )
    
    testthat::expect_equal(
      standardize_gene_name(rbcL$gene),
      "RBCL"
    )
    
    testthat::expect_equal(
      standardize_gene_name(matK$gene),
      "MATK"
    )
    
    testthat::expect_equal(
      standardize_gene_name(psbA$gene),
      "PSBA"
    )
    
  }
)