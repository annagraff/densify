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
      Warning:
      ! in `densify()`: using column `Glottocode` as taxon id
      i specify `taxon_id = <column name>` to silence this warning
      Warning in `densify()`:
      ! in `densify()`: no `cols` argument specified, using all columns as variables
      i use `cols = <tidy column spec>` to silence this warning

