## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error = (.Platform$OS.type == "windows"),
  eval = requireNamespace("usethis", quietly = TRUE)
)

set.seed(20201218)

## ----setup--------------------------------------------------------------------
library(mockr)

## ----fun-def------------------------------------------------------------------
access_resource <- function() {
  message("Trying to access resource...")
  # For some reason we can't access the resource in our tests.
  stop("Can't access resource now.")
}

work_with_resource <- function() {
  resource <- access_resource()
  message("Fetched resource: ", resource)
  invisible(resource)
}

## ----example-error, error = TRUE----------------------------------------------
try({
work_with_resource()
})

## ----example-remedy-----------------------------------------------------------
access_resource_for_test <- function() {
  # We return a value that's good enough for testing
  # and can be computed quickly:
  42
}

local({
  # Here, we override the function that raises the error
  local_mock(access_resource = access_resource_for_test)

  work_with_resource()
})

## ----work-around-desc-bug-1, echo = FALSE-------------------------------------
# Fixed in https://github.com/r-lib/desc/commit/daece0e5816e17a461969489bfdda2d50b4f5fe5, requires desc > 1.4.0
desc_options <- options(cli.num_colors = 1)

## ----create-package-----------------------------------------------------------
pkg <- usethis::create_package(file.path(tempdir(), "mocktest"))

## ----work-around-desc-bug-2, echo = FALSE-------------------------------------
options(desc_options)

## ----set-focus, include = FALSE-----------------------------------------------
wd <- getwd()

knitr::knit_hooks$set(
  pkg = function(before, options, envir) {
    if (before) {
      wd <<- setwd(pkg)
    } else {
      setwd(wd)
    }

    invisible()
  }
)

knitr::opts_chunk$set(pkg = TRUE)

## ----pkg-location-------------------------------------------------------------
usethis::proj_set()

## ----dir-tree-----------------------------------------------------------------
fs::dir_tree()

## ----run-pkg, error = TRUE----------------------------------------------------
try({
pkgload::load_all()
work_with_resource_pkg()
})

## ----test---------------------------------------------------------------------
usethis::use_testthat()

## ----error = TRUE-------------------------------------------------------------
try({
testthat::test_local(reporter = "location")
})

## ----test-manually------------------------------------------------------------
test_that("Can work with resource", {
  mockr::local_mock(access_resource_pkg = function() {
    42
  })

  expect_message(
    expect_equal(work_with_resource_pkg(), 42)
  )
})

## -----------------------------------------------------------------------------
pkgload::load_all()

## ----test-runif---------------------------------------------------------------
test_that("d6() works correctly", {
  seq <- c(0.32, 5.4, 5, 2.99)
  my_runif_mock <- function(...) {
    on.exit(seq <<- seq[-1])
    seq[[1]]
  }

  mockr::local_mock(my_runif = my_runif_mock)

  expect_equal(d6(), 1)
  expect_equal(d6(), 6)
  expect_equal(d6(), 6)
  expect_equal(d6(), 3)
})

