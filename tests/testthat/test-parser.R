test_file <- test_genbank_file(
  "NC_000932.gb"
)


testthat::test_that(
  "is_feature_start detects GenBank feature headers",
  {
    
    testthat::expect_true(
      is_feature_start(
        "     gene            1..100"
      )
    )
    
    testthat::expect_false(
      is_feature_start(
        "                     /gene=\"rbcL\""
      )
    )
    
  }
)


testthat::test_that(
  "new_feature creates an empty feature record",
  {
    
    feature <- new_feature()
    
    testthat::expect_true(
      is.list(feature)
    )
    
    testthat::expect_true(
      all(
        c(
          "type",
          "start",
          "end",
          "strand",
          "location_type",
          "location",
          "gene",
          "product",
          "locus_tag",
          "protein_id",
          "note"
        ) %in% names(feature)
      )
    )
    
  }
)


testthat::test_that(
  "parse_feature_header parses a feature header",
  {
    
    header <- parse_feature_header(
      "     gene            123..456"
    )
    
    testthat::expect_equal(
      header$type,
      "gene"
    )
    
    testthat::expect_equal(
      header$location,
      "123..456"
    )
    
  }
)


testthat::test_that(
  "parse_qualifier parses qualifiers",
  {
    
    qualifier <- parse_qualifier(
      '/gene="rbcL"'
    )
    
    testthat::expect_equal(
      qualifier$key,
      "gene"
    )
    
    testthat::expect_equal(
      qualifier$value,
      "rbcL"
    )
    
  }
)


testthat::test_that(
  "parse_location parses a simple interval",
  {
    
    location <- parse_location(
      "123..456"
    )
    
    testthat::expect_equal(
      location$start,
      123
    )
    
    testthat::expect_equal(
      location$end,
      456
    )
    
    testthat::expect_equal(
      location$strand,
      "+"
    )
    
  }
)


testthat::test_that(
  "parse_location parses a complement interval",
  {
    
    location <- parse_location(
      "complement(123..456)"
    )
    
    testthat::expect_equal(
      location$start,
      123
    )
    
    testthat::expect_equal(
      location$end,
      456
    )
    
    testthat::expect_equal(
      location$strand,
      "-"
    )
    
  }
)


testthat::test_that(
  "parse_location recognises joined locations",
  {
    
    location <- parse_location(
      "join(1..100,201..300)"
    )
    
    testthat::expect_true(
      location$location_type == "join"
    )
    
  }
)


testthat::test_that(
  "parse_location recognises complement joins",
  {
    
    location <- parse_location(
      "complement(join(1..100,201..300))"
    )
    
    testthat::expect_equal(
      location$strand,
      "-"
    )
    
  }
)


testthat::test_that(
  "split_feature_blocks returns feature blocks",
  {
    
    lines <- read_genbank(
      test_file
    )
    
    features <- extract_feature_section(
      lines
    )
    
    blocks <- split_feature_blocks(
      features
    )
    
    testthat::expect_true(
      is.list(blocks)
    )
    
    testthat::expect_gt(
      length(blocks),
      0
    )
    
  }
)


testthat::test_that(
  "parse_feature parses one feature block",
  {
    
    lines <- read_genbank(
      test_file
    )
    
    features <- extract_feature_section(
      lines
    )
    
    blocks <- split_feature_blocks(
      features
    )
    
    feature <- parse_feature(
      blocks[[1]]
    )
    
    testthat::expect_true(
      is.list(feature)
    )
    
    testthat::expect_true(
      "type" %in% names(feature)
    )
    
    testthat::expect_true(
      "start" %in% names(feature)
    )
    
    testthat::expect_true(
      "end" %in% names(feature)
    )
    
  }
)


testthat::test_that(
  "parse_features returns a feature table",
  {
    
    lines <- read_genbank(
      test_file
    )
    
    table <- parse_features(
      lines
    )
    
    testthat::expect_true(
      is.data.frame(table)
    )
    
    testthat::expect_gt(
      nrow(table),
      0
    )
    
    testthat::expect_true(
      all(
        c(
          "type",
          "start",
          "end",
          "strand",
          "location_type",
          "location",
          "gene"
        ) %in% names(table)
      )
    )
    
  }
)


testthat::test_that(
  "parse_location recognises origin-spanning complement joins",
  {
    
    location <- parse_location(
      "complement(join(70023..70028,1..69))"
    )
    
    testthat::expect_true(
      is.na(location$start)
    )
    
    testthat::expect_true(
      is.na(location$end)
    )
    
    testthat::expect_equal(
      location$strand,
      STRAND_REVERSE
    )
    
    testthat::expect_equal(
      location$location_type,
      LOCATION_JOIN
    )
    
  }
)


testthat::test_that(
  "parse_location recognises complex trans-spliced joins as joined locations",
  {
    
    location <- parse_location(
      "join(187159..188045,complement(167429..170407),complement(352420..356860))"
    )
    
    testthat::expect_true(
      is.na(location$start)
    )
    
    testthat::expect_true(
      is.na(location$end)
    )
    
    testthat::expect_equal(
      location$location_type,
      LOCATION_JOIN
    )
    
  }
)