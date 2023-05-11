
# A Snakekemake repository template

## Usage

Copy the content of the repository to an empty directory.
Then link your upstream repo and add the template:

```bash
git remote set-url origin https://git-repo/new-repository.git
git branch -M main
git add *
git add .gitignore .github/* .tests/*
git commit -am "Created repo"
git push -u origin main
```

Note that gh-actions will fail until you create some proper tests and setup 
gh-pages for the documentation.
