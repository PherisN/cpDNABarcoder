testthat::test_that(
  "is_single_string accepts valid single strings",
  {
    
    testthat::expect_true(
      is_single_string("rbcL")
    )
    
    testthat::expect_true(
      is_single_string("matK")
    )
    
    testthat::expect_true(
      is_single_string("")
    )
    
  }
)


testthat::test_that(
  "is_single_string rejects non-character objects",
  {
    
    testthat::expect_false(
      is_single_string(1)
    )
    
    testthat::expect_false(
      is_single_string(TRUE)
    )
    
    testthat::expect_false(
      is_single_string(list("rbcL"))
    )
    
    testthat::expect_false(
      is_single_string(factor("rbcL"))
    )
    
  }
)


testthat::test_that(
  "is_single_string rejects vectors",
  {
    
    testthat::expect_false(
      is_single_string(
        c("rbcL", "matK")
      )
    )
    
    testthat::expect_false(
      is_single_string(
        character(0)
      )
    )
    
  }
)


testthat::test_that(
  "is_single_string rejects missing values",
  {
    
    testthat::expect_false(
      is_single_string(NA_character_)
    )
    
    testthat::expect_false(
      is_single_string(NULL)
    )
    
  }
)