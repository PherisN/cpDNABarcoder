make_test_interval <- function(){
  
  new_interval(
    start = 10L,
    end = 20L,
    wraps_origin = FALSE,
    left_gene = "psbA",
    right_gene = "trnH"
  )
  
}


testthat::test_that(
  "new_interval constructs an Interval object",
  {
    
    interval <- make_test_interval()
    
    testthat::expect_true(
      inherits(
        interval,
        "Interval"
      )
    )
    
    testthat::expect_equal(
      interval$start,
      10L
    )
    
    testthat::expect_equal(
      interval$end,
      20L
    )
    
    testthat::expect_equal(
      interval$width,
      11L
    )
    
    testthat::expect_false(
      interval$wraps_origin
    )
    
    testthat::expect_equal(
      interval$left_gene,
      "psbA"
    )
    
    testthat::expect_equal(
      interval$right_gene,
      "trnH"
    )
    
  }
)


testthat::test_that(
  "new_interval constructs an origin-wrapping Interval",
  {
    
    interval <- new_interval(
      start = 150000L,
      end = 100L,
      wraps_origin = TRUE,
      left_gene = "trnH",
      right_gene = "psbA"
    )
    
    testthat::expect_true(
      inherits(
        interval,
        "Interval"
      )
    )
    
    testthat::expect_true(
      interval$wraps_origin
    )
    
    testthat::expect_true(
      is.na(
        interval$width
      )
    )
    
  }
)


testthat::test_that(
  "new_interval rejects invalid inputs",
  {
    
    testthat::expect_error(
      new_interval(
        start = 0L,
        end = 20L,
        wraps_origin = FALSE,
        left_gene = "psbA",
        right_gene = "trnH"
      )
    )
    
    testthat::expect_error(
      new_interval(
        start = 10L,
        end = NA_integer_,
        wraps_origin = FALSE,
        left_gene = "psbA",
        right_gene = "trnH"
      )
    )
    
    testthat::expect_error(
      new_interval(
        start = 10L,
        end = 20L,
        wraps_origin = c(TRUE, FALSE),
        left_gene = "psbA",
        right_gene = "trnH"
      )
    )
    
    testthat::expect_error(
      new_interval(
        start = 10L,
        end = 20L,
        wraps_origin = FALSE,
        left_gene = c("psbA", "matK"),
        right_gene = "trnH"
      )
    )
    
    testthat::expect_error(
      new_interval(
        start = 10L,
        end = 20L,
        wraps_origin = FALSE,
        left_gene = "psbA",
        right_gene = NA_character_
      )
    )
    
  }
)


testthat::test_that(
  "validate_interval rejects invalid Interval objects",
  {
    
    interval <- make_test_interval()
    
    missing_field <- interval
    missing_field$width <- NULL
    
    testthat::expect_error(
      validate_interval(
        missing_field
      )
    )
    
    invalid_coordinates <- interval
    invalid_coordinates$start <- 0L
    
    testthat::expect_error(
      validate_interval(
        invalid_coordinates
      )
    )
    
    reversed_interval <- interval
    reversed_interval$start <- 30L
    reversed_interval$end <- 20L
    reversed_interval$width <- -9L
    
    testthat::expect_error(
      validate_interval(
        reversed_interval
      )
    )
    
    inconsistent_width <- interval
    inconsistent_width$width <- 99L
    
    testthat::expect_error(
      validate_interval(
        inconsistent_width
      )
    )
    
  }
)


testthat::test_that(
  "assert_interval accepts Interval objects and rejects others",
  {
    
    interval <- make_test_interval()
    
    testthat::expect_identical(
      assert_interval(
        interval
      ),
      interval
    )
    
    testthat::expect_error(
      assert_interval(
        list()
      )
    )
    
  }
)


testthat::test_that(
  "Interval accessors return expected values",
  {
    
    interval <- make_test_interval()
    
    testthat::expect_equal(
      get_interval_start(interval),
      10L
    )
    
    testthat::expect_equal(
      get_interval_end(interval),
      20L
    )
    
    testthat::expect_equal(
      get_interval_width(interval),
      11L
    )
    
    testthat::expect_equal(
      get_left_gene(interval),
      "psbA"
    )
    
    testthat::expect_equal(
      get_right_gene(interval),
      "trnH"
    )
    
    testthat::expect_false(
      interval_wraps_origin(interval)
    )
    
  }
)


testthat::test_that(
  "print.Interval prints a concise interval summary",
  {
    
    interval <- make_test_interval()
    
    output <- capture.output(
      print(
        interval
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Interval Object",
          output
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "psbA",
          output
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "trnH",
          output
        )
      )
    )
    
  }
)