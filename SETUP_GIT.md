# Git setup for this repo

Git is now installed as a **portable MinGit** at
`%LOCALAPPDATA%\Tools\MinGit\cmd\git.exe` (added to your **User PATH**).
**Open a brand-new terminal/VS Code window** for the PATH change to take
effect — then plain `git` commands work without a full path.

This repo has already been initialized and committed locally (`git init`,
`git add .`, first commit). What's left is pushing it somewhere.

## Push to GitHub

This repo does not yet have a remote configured — no GitHub URL was given
for it (only the TCAMS app repo had an existing URL). To publish it:

```powershell
# 1. Create a new, empty repository on GitHub first (e.g. via github.com/new
#    or `gh repo create`), WITHOUT initializing it with a README/.gitignore.

# 2. From this folder:
git remote add origin https://github.com/<your-username>/<repo-name>.git
git push -u origin main
```

## Notes specific to this session's environment

- A stray double-quote had corrupted the `HOME` **User** environment
  variable (`H:"` instead of `H:\`), which made every `git` invocation fail
  with `fatal: unable to access '...': Invalid argument`. This was fixed at
  the registry level. If you ever see that error again in a *new* terminal,
  check `[Environment]::GetEnvironmentVariable("HOME","User")` in PowerShell.
- `G:\` and other mapped/network drives can trigger git's "detected dubious
  ownership" safety check. This repo's path was already added as an
  exception via `git config --global --add safe.directory
  G:/Thesis/3rdChapter/PEP_QC_HMSC_Analysis`. If you move the repo to a
  different path, re-run that command with the new path.
- You may see `LF will be replaced by CRLF` warnings on commit — harmless,
  just Git normalizing line endings on Windows (`core.autocrlf=true` is set
  globally).
