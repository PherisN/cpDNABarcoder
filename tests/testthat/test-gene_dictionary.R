testthat::test_that(
  "gene aliases are non-empty character vectors",
  {
    
    testthat::expect_true(
      is.character(
        GENE_DICTIONARY$RBCL$aliases
      )
    )
    
    testthat::expect_gte(
      length(
        GENE_DICTIONARY$RBCL$aliases
      ),
      1
    )
    
    testthat::expect_true(
      is.character(
        GENE_DICTIONARY$MATK$aliases
      )
    )
    
    testthat::expect_gte(
      length(
        GENE_DICTIONARY$MATK$aliases
      ),
      1
    )
    
    testthat::expect_true(
      is.character(
        GENE_DICTIONARY$PSBA$aliases
      )
    )
    
    testthat::expect_gte(
      length(
        GENE_DICTIONARY$PSBA$aliases
      ),
      1
    )
    
    testthat::expect_true(
      is.character(
        GENE_DICTIONARY$TRNH$aliases
      )
    )
    
    testthat::expect_gte(
      length(
        GENE_DICTIONARY$TRNH$aliases
      ),
      1
    )
    
  }
)