# Naming audit, design rules section 3. The rules are checked by
# bedrock::auditNames(); what is listed here is accepted, each entry with
# its reason. An entry starting with "open:" is a name still to be
# decided - it documents a debt and is removed when the name changes.

.auditNamesExceptions <- c(
  "hosmerLemeshowTest.default(type = \"C\")"   = "the names of the two statistics in Hosmer and Lemeshow",
  "hosmerLemeshowTest.default(type = \"H\")"   = "the names of the two statistics in Hosmer and Lemeshow",
  "hosmerLemeshowTest.glm(type = \"C\")"       = "the names of the two statistics in Hosmer and Lemeshow",
  "hosmerLemeshowTest.glm(type = \"H\")"       = "the names of the two statistics in Hosmer and Lemeshow"
)


test_that("exported names and arguments follow the naming rules", {

  res <- bedrock::auditNames("lumen",
                             exceptions = names(.auditNamesExceptions))

  expect_identical(
    nrow(res), 0L,
    info = paste(sprintf("%s [%s] %s", res$key, res$rule, res$detail),
                 collapse = "\n"))
})


test_that("every accepted exception still matches a finding", {

  res <- bedrock::auditNames("lumen",
                             exceptions = names(.auditNamesExceptions))

  expect_identical(attr(res, "unused"), character(0))
})
