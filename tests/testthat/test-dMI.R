############
## set up ##
############
# valid, symmetric n x n distance setup (n = 4 samples)
set.seed(1)
valid_D.x_matrix = as.matrix(dist(matrix(rnorm(4 * 2), nrow = 4)))
valid_D.y_matrix = as.matrix(dist(matrix(rnorm(4 * 2), nrow = 4)))
valid_D.x_dist   = as.dist(valid_D.x_matrix)
valid_D.y_dist   = as.dist(valid_D.y_matrix)

# square but not symmetric matrix (for symmetry check)
asymmetric_matrix = matrix(c(0, 1, 2,
                             3, 0, 4,
                             5, 6, 0), nrow = 3, byrow = TRUE)

# non-square matrix (for square check)
nonsquare_matrix = matrix(1:6, nrow = 2, ncol = 3)

# square, symmetric distance matrix but with a different n (for the D.x/D.y dimension-match check)
mismatched_n_matrix = as.matrix(dist(matrix(rnorm(3 * 2), nrow = 3)))

# neither class 'dist' nor 'matrix' (for class check)
not_a_matrix = data.frame(x = 1:3, y = 4:6)


###############
## D.x tests ##
###############

test_that("D.x must be class 'dist' or 'matrix'", {
  expect_error(
    dMI(D.x = not_a_matrix, D.y = valid_D.y_matrix, kbar = 1:3),
    "D.x must be a distance matrix",
    fixed = TRUE
  )
})

test_that("D.x must be square", {
  expect_error(
    dMI(D.x = nonsquare_matrix, D.y = valid_D.y_matrix, kbar = 1:3),
    "D.x must be a square",
    fixed = TRUE
  )
})

test_that("D.x must be symmetric", {
  expect_error(
    dMI(D.x = asymmetric_matrix, D.y = valid_D.y_matrix, kbar = 1:3),
    "D.x must be a symmetric",
    fixed = TRUE
  )
})


###############
## D.y tests ##
###############

test_that("D.y must be class 'dist' or 'matrix'", {
  expect_error(
    dMI(D.x = valid_D.x_matrix, D.y = not_a_matrix, kbar = 1:3),
    "D.y must be a distance matrix",
    fixed = TRUE
  )
})

test_that("D.y must be square", {
  expect_error(
    dMI(D.x = valid_D.x_matrix, D.y = nonsquare_matrix, kbar = 1:3),
    "D.y must be a square",
    fixed = TRUE
  )
})

test_that("D.y must be symmetric", {
  expect_error(
    dMI(D.x = valid_D.x_matrix, D.y = asymmetric_matrix, kbar = 1:3),
    "D.y must be a symmetric",
    fixed = TRUE
  )
})


######################################
## D.x/D.y matching dimension tests ##
######################################
# mismatched_n_matrix is valid (square, symmetric, correct class)
# only has a different n than valid_D.x_matrix/valid_D.y_matrix.

test_that("D.x and D.y must have the same dimension", {
  expect_error(
    dMI(D.x = valid_D.x_matrix, D.y = mismatched_n_matrix, kbar = 1:3),
    "same dimension",
    fixed = TRUE
  )
})

test_that("D.x and D.y dimension check is symmetric (order doesn't matter)", {
  expect_error(
    dMI(D.x = mismatched_n_matrix, D.y = valid_D.y_matrix, kbar = 1:3),
    "same dimension",
    fixed = TRUE
  )
})


################
## kbar tests ##
################

test_that("kbar must have length > 1", {
  expect_error(
    dMI(D.x = valid_D.x_matrix, D.y = valid_D.y_matrix, kbar = 3),
    "kbar must be a vector",
    fixed = TRUE
  )
})

test_that("kbar of length 1 passed as a vector still errors", {
  # length() == 1 regardless of how it's constructed
  expect_error(
    dMI(D.x = valid_D.x_matrix, D.y = valid_D.y_matrix, kbar = c(3)),
    "kbar must be a vector",
    fixed = TRUE
  )
})

test_that("kbar cannot contain values larger than n(n-1)/2", {
  n <- nrow(valid_D.x_matrix)
  too_large_k <- (n * (n - 1) / 2) + 1
  expect_error(
    dMI(D.x = valid_D.x_matrix, D.y = valid_D.y_matrix, kbar = 1:too_large_k),
    "kbar cannot contain values larger than n(n-1)/2",
    fixed = TRUE
  )
})

test_that("kbar at exactly n(n-1)/2 is allowed (boundary is inclusive)", {
  n <- nrow(valid_D.x_matrix)
  max_k <- n * (n - 1) / 2
  expect_failure(
    expect_error(
      dMI(D.x = valid_D.x_matrix, D.y = valid_D.y_matrix, kbar = 1:max_k),
      "kbar cannot contain values larger than n(n-1)/2",
      fixed = TRUE
    )
  )
})

#######################
## valid input tests ##
#######################

test_that("valid matrix inputs do not trigger any of the known validation errors", {
  expect_no_error(
    dMI(D.x = valid_D.x_matrix, D.y = valid_D.y_matrix, kbar = 1:3)
  )
})

test_that("valid 'dist' class inputs pass the class/square/symmetric checks", {
  expect_no_error(
    dMI(D.x = valid_D.x_dist, D.y = valid_D.y_dist, kbar = 1:3)
  )
})
