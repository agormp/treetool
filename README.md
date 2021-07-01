# TreeTool: python script for

Script used for data analysis in paper "XXX".

Converts between different tree file formats. Performs various manipulations and analyses on trees

## Installation

Requires installation of [treelib.py](https://github.com/agormp/treelib) module

Place `treetool.py` in directory that is pointed to by the PATH environment variable, and make script executable. For instance (on MacOS):
```
mv treetool.py /usr/local/bin
chmod 755 /usr/local/bin/treetool.py
```

## Usage

The script should be run on the command line. It requires one or more tree files.

The following command prints help for command line usage:
```
treetool.py -h
```

For example, the command below reads a Newick formatted tree file, and outputs a NEXUS formatted tree file to stdout (which is then redirected to a file):

```
treetool.py -I newick -O nexus mytreefile.newick > mytreefile.nexus
```

Below is a list of available options (output from -h option):

```
Usage: treetool.py [options] TREEFILE

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -I FORMAT             Input format: nexus or newick  [default: nexus]
  -O FORMAT             Output format: nexus or newick
  -m, --midpoint        Root tree at midpoint
  --root=FILE           Root tree using taxa in FILE as outgroup (one name per
                        line)
  --clustcut=ROOTDIST   Split leaves into clusters by cutting across tree
                        ROOTDIST from root. Results are written to file
                        clusterinfo.txt in dir specified using option
                        --clustdir
  --clustn=N            Split leaves into N clusters by cutting across tree at
                        suitable distance from root. Results are written to
                        file clusterinfo.txt in dir specified using option
                        --clustdir
  --clustrees           (Requires --clustn or clustcut): output treefiles for
                        each cluster, and one treefile showing cluster
                        placement
  --clustdir=DIRNAME    Name of directory in which cluster result files are
                        placed. [default: clusterdir]
  --cladeinfo=FILE      For leaves listed in FILE, reports: Is this a
                        monophyletic group? What is distance from root to
                        enveloping clade's MRCA and MRCA's parent?
  --rename=FILE         Rename leaves using information in FILE (format: old
                        new)
  --shortnames=B[,E]    Abbreviate names: Keep B first (and E last) characters
                        in each name
  --nameprune=SEP[,KEEP]
                        Prune clades containing multiple taxa with identical
                        name starts (up to first SEP). KEEP: do not prune
                        names containing this pattern
  --prune=FILE          Remove leaves listed in FILE (one leaf name per line)
  --debug               Print longer error messages
  --overwrite           Overwrite identically named previous resultfiles with
                        no warning (e.g. for clustering options)
```