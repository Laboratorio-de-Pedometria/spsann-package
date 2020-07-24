# Check the proportion of accepted jitters in the first chain
#
# COMMAND:
# eval(.check_first_chain())
#
# SUMMARY
# 1. Check if we are at the end of the first chain;
# 2. Compute the proportion of accepted jitters;
# 3. Break if the proportion of accepted jitters is lowwer than
#    'initial.acceptance';
# 4. Continue and return the number of accepted jitters otherwise.
#
# NOTES
# 1.
#
if (i == 1) {
  x <- round(n_accept / c(n_pts * schedule$chain.length), 2)
  if (x < schedule$initial.acceptance[1]) {
    cat(
      "\nthe initial temperature is too low: ", round(x * 100, 2),
      "% of acceptance in the 1st Markov chain\n", sep = "")
    break
  } else if (x > schedule$initial.acceptance[2]) {
    cat(
      "\nthe initial temperature is too high: more than ",
      round(schedule$initial.acceptance[2] * 100, 2),
      "% of acceptance in the 1st Markov chain\n", sep = "")
    break
  } else {
    cat(
      "\n", round(x * 100, 2),
      "% of acceptance in the 1st Markov chain\n", sep = "")
  }
}
