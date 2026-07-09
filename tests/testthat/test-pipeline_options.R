testthat::test_that(
  "default pipeline options validate successfully",
  {
    
    testthat::expect_identical(
      validate_pipeline_options(
        PIPELINE_OPTIONS
      ),
      PIPELINE_OPTIONS
    )
    
  }
)


testthat::test_that(
  "pipeline options are returned as a list",
  {
    
    testthat::expect_true(
      is.list(
        PIPELINE_OPTIONS
      )
    )
    
  }
)


testthat::test_that(
  "required option names exist",
  {
    
    required <- c(
      "reverse_complement",
      "include_flanks",
      "flank_length",
      "allow_wraparound",
      "export_coordinates",
      "verbose"
    )
    
    testthat::expect_true(
      all(
        required %in%
          names(PIPELINE_OPTIONS)
      )
    )
    
  }
)


testthat::test_that(
  "pipeline option types are correct",
  {
    
    logical_options <- c(
      "reverse_complement",
      "include_flanks",
      "allow_wraparound",
      "export_coordinates",
      "verbose"
    )
    
    for(name in logical_options){
      
      value <- PIPELINE_OPTIONS[[name]]
      
      testthat::expect_true(
        is.logical(value)
      )
      
      testthat::expect_length(
        value,
        1
      )
      
      testthat::expect_false(
        is.na(value)
      )
      
    }
    
    testthat::expect_true(
      is.numeric(
        PIPELINE_OPTIONS$flank_length
      )
    )
    
    testthat::expect_length(
      PIPELINE_OPTIONS$flank_length,
      1
    )
    
    testthat::expect_gte(
      PIPELINE_OPTIONS$flank_length,
      0
    )
    
  }
)


testthat::test_that(
  "validate_pipeline_options returns validated options",
  {
    
    options <- validate_pipeline_options(
      PIPELINE_OPTIONS
    )
    
    testthat::expect_identical(
      options,
      PIPELINE_OPTIONS
    )
    
  }
)