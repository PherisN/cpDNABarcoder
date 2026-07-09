testthat::test_that(
  "MARKER_DICTIONARY validates successfully",
  {
    
    testthat::expect_invisible(
      validate_marker_dictionary(
        MARKER_DICTIONARY
      )
    )
    
  }
)


testthat::test_that(
  "required markers exist",
  {
    
    testthat::expect_true(
      "RBCL" %in% names(MARKER_DICTIONARY)
    )
    
    testthat::expect_true(
      "MATK" %in% names(MARKER_DICTIONARY)
    )
    
    testthat::expect_true(
      "PSBA_TRNH" %in% names(MARKER_DICTIONARY)
    )
    
  }
)


testthat::test_that(
  "marker classes are valid",
  {
    
    valid_classes <- c(
      MARKER_CLASS_GENE,
      MARKER_CLASS_INTERGENIC,
      MARKER_CLASS_NUCLEAR
    )
    
    for(marker in MARKER_DICTIONARY){
      
      testthat::expect_true(
        marker$class %in% valid_classes
      )
      
    }
    
  }
)


testthat::test_that(
  "marker names are single strings",
  {
    
    for(marker in MARKER_DICTIONARY){
      
      testthat::expect_true(
        is_single_string(
          marker$marker
        )
      )
      
    }
    
  }
)


testthat::test_that(
  "genome types are valid",
  {
    
    valid_types <- c(
      "plastid",
      "nuclear"
    )
    
    for(marker in MARKER_DICTIONARY){
      
      testthat::expect_true(
        marker$type %in% valid_types
      )
      
    }
    
  }
)


testthat::test_that(
  "gene markers define a target gene",
  {
    
    gene_markers <- names(
      Filter(
        function(x)
          x$class == MARKER_CLASS_GENE,
        MARKER_DICTIONARY
      )
    )
    
    for(name in gene_markers){
      
      testthat::expect_true(
        is_single_string(
          MARKER_DICTIONARY[[name]]$gene
        )
      )
      
    }
    
  }
)


testthat::test_that(
  "intergenic markers define flanking genes",
  {
    
    intergenic_markers <- names(
      Filter(
        function(x)
          x$class == MARKER_CLASS_INTERGENIC,
        MARKER_DICTIONARY
      )
    )
    
    for(name in intergenic_markers){
      
      testthat::expect_true(
        is_single_string(
          MARKER_DICTIONARY[[name]]$left_gene
        )
      )
      
      testthat::expect_true(
        is_single_string(
          MARKER_DICTIONARY[[name]]$right_gene
        )
      )
      
    }
    
  }
)


testthat::test_that(
  "nuclear markers define a locus",
  {
    
    nuclear_markers <- names(
      Filter(
        function(x)
          x$class == MARKER_CLASS_NUCLEAR,
        MARKER_DICTIONARY
      )
    )
    
    if(length(nuclear_markers) == 0){
      
      testthat::skip(
        "No nuclear markers are currently implemented."
      )
      
    }
    
    for(name in nuclear_markers){
      
      testthat::expect_true(
        is_single_string(
          MARKER_DICTIONARY[[name]]$locus
        )
      )
      
    }
    
  }
)