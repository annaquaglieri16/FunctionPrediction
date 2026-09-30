# Using Claude subagents alongside pytest on this project

Two complementary layers of "testing," matched to the learning goals for
this repo:

1. **Deterministic biology logic -> pytest.** Anything with a
   known-correct answer (enzyme cleavage rules, missed-cleavage counts,
   fragment coordinates, FASTA parsing) gets a pytest test in `tests/`.
   This is ground truth -- an agent's opinion doesn't override a failing
   assertion here.

2. **The engineering loop -> subagents.** For everything less
   deterministic (is this the right embedding-pooling strategy, does this
   evaluation script actually answer the research question, is this
   Docker setup sane), use a small set of roles:
   - a *research/explore* agent to survey the codebase or literature
     before changing something,
   - a *build* agent that writes the change,
   - a *verify/review* agent, run separately from the build step, that
     checks the diff against the actual requirement (not just "does it
     run") -- this is the one worth keeping independent, since a reviewer
     that shares context with the builder tends to rubber-stamp its own
     work.

For a solo learning project the simplest version of this is: ask for a
change, then in a *separate* follow-up ask explicitly for a critical
review of what was just built before treating it as done.
