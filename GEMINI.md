# Instructions for AI coding agents

Read `AGENT.md` before making changes. Section 6 ("Shared Infrastructure: CI, Tests & Commits") is mandatory for every commit:

- Start from the current `origin/main` (`git fetch origin && git rebase origin/main`). Never copy files in from an older checkout or worktree; that reverts other people's changes.
- Before committing, run `git diff --stat origin/main` and restore any file you did not mean to change (`git checkout origin/main -- <file>`).
- Do not edit `.github/workflows/` or the `deep` test labels in `CMakeLists.txt` unless the task is about CI. CI is two-tier on purpose: a fast tier on every push and PR (`ctest -LE deep`) and a deep tier nightly or on PRs labelled `full-ci`.
- New tests that take a minute or more get the `deep` label.
- Test data must be inside this repo (`tests/data/`), never `../hyphy/...`.
