# densify without cols specified (should show a warning)

    Code
      result <- densify(WALS[1:10, 1:10])
    Condition
      Warning in `densify()`:
      ! in `densify()`: no `cols` argument specified, using all columns as variables
      i use `cols = <tidy column spec>` to silence this warning

# densify without cols specified (should show a warning and guess Glottocode)

    Code
      result <- densify(WALS[1:10, 1:10], taxonomy = glottolog_languoids)
    Condition
      Warning in `densify()`:
      ! in `densify()`: using column `Glottocode` as taxon id
      i specify `taxon_id = <column name>` to silence this warning
      Warning in `densify()`:
      ! in `densify()`: no `cols` argument specified, using all columns as variables
      i use `cols = <tidy column spec>` to silence this warning

# densify with empty data (should give a warning)

    Code
      result <- densify(WALS[0, , drop = FALSE])
    Condition
      Warning in `densify()`:
      ! in `densify()`: empty dataset supplied
      Warning in `densify()`:
      ! in `densify()`: no `cols` argument specified, using all columns as variables
      i use `cols = <tidy column spec>` to silence this warning

# densify with empty column selection (should give a warning)

    Code
      result <- densify(WALS, tidyselect::starts_with("no_such_column"))
    Condition
      Warning in `densify()`:
      ! in `densify()`: no variables selected, returning empty data frame
      i please adjust `cols = <tidy column spec>`

# densify with invalid columns (should result in an error)

    Code
      result <- densify(WALS, c(no_such_column1, no_such_column2))
    Condition
      Error in `densify()`:
      ! Can't subset columns that don't exist.
      x Column `no_such_column1` doesn't exist.

