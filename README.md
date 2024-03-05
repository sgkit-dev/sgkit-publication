# vcf-zarr-publication

This repo contains the manuscript for the publication describing the vcf-zarr 
specification and its compression and query performance on several datasets.
All code required to generate figures and example analyses.

## Contribution

Writing occurs either via this repo and the standard GitHub workflows,
or using the linked
[Overleaf project](https://www.overleaf.com/project/659c19afc265913df38df1ab).
Both are useful for different things, and so will be used in parallel.

Please contact Jerome Kelleher if you want to get edit access on
Overleaf.

## Contributing via Overleaf

GitHub syncing on Overleaf can be somewhat fragile, and the process needs
to be manually curated (Jerome Kelleher will manage).
From previous experience, the following protocols work reasonably well:

- **DO NOT USE** Overleaf's "Visual Editor". This can change the source LaTeX
  in arbitrary ways, and really muck things up for those of us using
  line-oriented editors. Please stick to the default "Code Editor" view.
- Only edit the files ``paper.tex`` and ``paper.bib`` on Overleaf. For anything
  else, please use GitHub.


## General contributing guidelines

- Please insert line-breaks before 90-100 characters, ideally at semantic
  break points. (Less necessary on Overleaf, but it makes keeping track of changes
  on git easier, and also helps people using old-school editors).

- When adding to the bibliography, please use Google Scholar's key format
  (e.g. `garud2015recent`), or preferably use Google Scholar's bibtex export
  (click the "cite" link for a particular paper, and the Bibtex option is on the
  bottom).


## Use of Issues/Discussions/Overleaf comments

All of these routes are useful for different things. Here are some rough
guidelines on when they might be used:

- GitHub Discussions should be used for high-level topics, where there isn't
  any particular action to be taken and the threaded nature can be useful
  (e.g., https://github.com/sgkit-dev/sgkit-publication/discussions/51)
- GitHub issues should be designed to be closed by some action---issues
  keep track of specific problems that need to be resolved. Ultimately
  we should have 0 issues open when the paper is submitted/published.
- Overleaf comments are useful for short comments on the text. These should
  also be designed to be closed, as the value of text comments drops quickly
  when there are lots lying around unresolved.
  Any longer discussions arising from the text should move to GitHub.

