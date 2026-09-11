#!/usr/bin/env Rscript
# "99 Bottles of Beer" -- traditional counting song, public domain.
# Not part of mytools; a standalone warm-up script.

# Count as a noun phrase: 0 -> "no more bottles", 1 -> "1 bottle", n -> "n bottles".
bottles <- function(n) {
  if (n == 0) "no more bottles" else sprintf("%d %s", n, if (n == 1) "bottle" else "bottles")
}

# Same, but capitalised for the start of a line.
Bottles <- function(n) {
  s <- bottles(n)
  paste0(toupper(substring(s, 1, 1)), substring(s, 2))
}

start <- 99
out <- character(0)

for (n in start:1) {
  out <- c(out,
    sprintf("%s of beer on the wall, %s of beer.", Bottles(n), bottles(n)),
    sprintf("Take one down and pass it around, %s of beer on the wall.", bottles(n - 1)),
    "")
}

out <- c(out,
  sprintf("%s of beer on the wall, %s of beer.", Bottles(0), bottles(0)),
  sprintf("Go to the store and buy some more, %s of beer on the wall.", bottles(start)))

writeLines(out)
