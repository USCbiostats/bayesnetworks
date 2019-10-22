test_that("inconsistant input yields errors", {

  expect_error(create_network(1, "a"), "same type")
  expect_error(create_network(1:2, 1:3), "same length")
  expect_error(create_network(1, 1), "same")
  expect_error(create_network(node_type = "sink"), "cannot be specified if")
  expect_error(create_network(1, 2, "A"))

  # test that source and target is consistant with node_type
})

test_that("empty network creation", {
  network <- create_network(numeric(), numeric())

  expect_s3_class(network, "bayesnetworks_network")
  expect_equal(vapply(network, typeof, FUN.VALUE = character(1), USE.NAMES = FALSE),
               c("integer", "integer", "double", "character"))
  expect_equal(vapply(network, length, FUN.VALUE = numeric(1), USE.NAMES = FALSE),
               rep(0, 4))

  network <- create_network(character(), character())

  expect_s3_class(network, "bayesnetworks_network")
  expect_equal(vapply(network, typeof, FUN.VALUE = character(1), USE.NAMES = FALSE),
               c("integer", "integer", "character", "character"))
  expect_equal(vapply(network, length, FUN.VALUE = numeric(1), USE.NAMES = FALSE),
               rep(0, 4))
})

test_that("network with no egdes creation", {
  network <- create_network(node_labels = 1)

  expect_equal(vapply(network, typeof, FUN.VALUE = character(1), USE.NAMES = FALSE),
               c("integer", "integer", "double", "character"))
  expect_equal(vapply(network, length, FUN.VALUE = numeric(1), USE.NAMES = FALSE),
               c(0, 0, 1, 1))

  network <- create_network(node_labels = "1")

  expect_equal(vapply(network, typeof, FUN.VALUE = character(1), USE.NAMES = FALSE),
               c("integer", "integer", "character", "character"))
  expect_equal(vapply(network, length, FUN.VALUE = numeric(1), USE.NAMES = FALSE),
               c(0, 0, 1, 1))

  network <- create_network(node_labels = 1:100)

  expect_equal(vapply(network, typeof, FUN.VALUE = character(1), USE.NAMES = FALSE),
               c("integer", "integer", "integer", "character"))
  expect_equal(vapply(network, length, FUN.VALUE = numeric(1), USE.NAMES = FALSE),
               c(0, 0, 100, 100))
})

test_that("1 egde network creation", {
  network <- create_network(1, 2)

  expect_equal(vapply(network, typeof, FUN.VALUE = character(1), USE.NAMES = FALSE),
               c("integer", "integer", "double", "character"))
  expect_equal(vapply(network, length, FUN.VALUE = numeric(1), USE.NAMES = FALSE),
               c(1, 1, 2, 2))

  network <- create_network("A", "B")

  expect_equal(vapply(network, typeof, FUN.VALUE = character(1), USE.NAMES = FALSE),
               c("integer", "integer", "character", "character"))
  expect_equal(vapply(network, length, FUN.VALUE = numeric(1), USE.NAMES = FALSE),
               c(1, 1, 2, 2))
})


test_that("larger network creation", {
  network <- create_network(LETTERS[-1], rep("A",  25))

  expect_equal(vapply(network, typeof, FUN.VALUE = character(1), USE.NAMES = FALSE),
               c("integer", "integer", "character", "character"))
  expect_equal(vapply(network, length, FUN.VALUE = numeric(1), USE.NAMES = FALSE),
               c(25, 25, 26, 26))
})
