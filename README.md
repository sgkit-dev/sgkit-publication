# msprime-1.0-paper

This repo contains the manuscript for the publication describing the 
msprime 1.0 paper. 

## Abstract

Coalescent simulation is a key tool in population genetics and
and has been integral to coalescent theory since its earliest days.
Because of the ease at which the basic model can be simulated,
a large number of different simulators have been developed. However,
the coalescent with recombination is far more challenging to simulate,
and, until recently, it was not possible to simulate efficiently.
The msprime software has revolutionised
coalescent simulation, making it possible to simulate millions
of whole genomes for the first time. We summarise the many features
of msprime 1.0, which is built around a core model of efficiently
implementing recombination via the succinct tree sequence data
structure. We advocate a community oriented open-source development
process as a way to reduce duplication of effort and increase
the quality of simulation software.

## Writing process

Like the development of msprime itself, the process of writing this 
article is managed through an open source development model via
GitHub. Those responsible for particular features will be 
assigned sections of writing via issues. Those wishing to make 
edits should also follow these guidelines:

- For significant edits, please use the GitHub workflow described 
  [here](https://stdpopsim.readthedocs.io/en/latest/development.html#github-workflow)
  to create a fork, feature branches and make pull requests.

- For smaller edits, using the 
  [edit file](https://help.github.com/en/github/managing-files-in-a-repository/editing-files-in-your-repository) 
  button on the paper.tex file in GitHub may be simpler than checking out 
  the repo locally, etc.

- Please keep edits localised to a single section to make the 
  review process simpler.

- Please insert line-breaks before 90-100 characters, ideally at semantic
  break points.
  
- For high-level/structural issues, please open an issue to discuss before making changes.


