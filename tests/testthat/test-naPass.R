# The formula methods of the k-sample tests default to na.action = na.pass
# and leave the missing values to the default method (design rules, NA
# policy, scope: tests). That must not change a single number compared with
# removing them beforehand.

.naPassData <- data.frame(
  val = c(2.9, NA,  2.5, 2.6, 3.2,
          3.8, 2.7, 4.0, 2.4,
          2.8, 3.4, 3.7, 2.2, 2.0, 3.1),
  grp = factor(c(rep(c("X", "Y", "Z"), c(5, 4, 5)), NA))
)

.naPassTests <- list(
  conoverTest            = conoverTest,
  dscfTest               = dscfTest,
  dunnTest               = dunnTest,
  dunnettTest            = dunnettTest,
  jonckheereTerpstraTest = jonckheereTerpstraTest,
  nemenyiTest            = nemenyiTest,
  steelTest              = steelTest,
  vanDerWaerdenTest      = vanDerWaerdenTest
)


test_that("na.pass is the default of the k-sample formula methods", {

  for (nm in names(.naPassTests))
    expect_identical(
      formals(getS3method(nm, "formula"))$na.action, quote(na.pass),
      info = nm)
})


test_that("the default method drops incomplete cases: na.pass equals na.omit", {

  for (nm in names(.naPassTests)) {

    f <- .naPassTests[[nm]]

    viaPass <- suppressWarnings(f(val ~ grp, data = .naPassData))
    viaOmit <- suppressWarnings(f(val ~ grp, data = .naPassData,
                                  na.action = na.omit))

    expect_equal(viaPass, viaOmit, info = nm)
  }
})
