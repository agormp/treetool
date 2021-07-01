#!/usr/bin/env python3
# By Anders Gorm Pedersen, agpe@dtu.dk, Technical University of Denmark
# Converts between different tree file formats. Performs various manipulations and analyses on trees

import sys, os.path, math, shutil, treelib
from optparse import OptionParser

################################################################################################

def build_parser():
    parser = OptionParser(usage="usage: %prog [options] TREEFILE",
                          version="0.5")

    parser.add_option("-I", type="choice", dest="informat",
                    choices=["NEXUS", "nexus", "NEWICK", "newick"], metavar="FORMAT",
                    help="Input format: nexus or newick  [default: nexus]")

    parser.add_option("-O",type="choice", dest="outformat",
                    choices=["NEXUS", "nexus", "NEWICK", "newick"], metavar="FORMAT",
                    help="Output format: nexus or newick")

    parser.add_option("-m", "--midpoint", action="store_true", dest="midroot",
                          help="Root tree at midpoint")

    parser.add_option("--root", action="store", type="string", dest="rootfile", metavar="FILE",
                      help="Root tree using taxa in FILE as outgroup (one name per line)")

    parser.add_option("--clustcut", action="store", type="float", dest="clustcut", metavar="ROOTDIST",
                      help="Split leaves into clusters by cutting across tree ROOTDIST from root. Results are written to file clusterinfo.txt in dir specified using option --clustdir")

    parser.add_option("--clustn", action="store", type="int", dest="clustn", metavar="N",
                      help="Split leaves into N clusters by cutting across tree at suitable distance from root. Results are written to file clusterinfo.txt in dir specified using option --clustdir")

    parser.add_option("--clustrees", action="store_true", dest="clustrees",
                      help="(Requires --clustn or clustcut): output treefiles for each cluster, and one treefile showing cluster placement")

    parser.add_option("--clustdir", action="store", type="string", dest="clustdir", metavar="DIRNAME",
                      help="Name of directory in which cluster result files are placed. [default: clusterdir]")

    parser.add_option("--cladeinfo", action="store", type="string", dest="cladeinfo", metavar="FILE",
                      help="For leaves listed in FILE, reports: Is this a monophyletic group? What is distance from root to enveloping clade's MRCA and MRCA's parent?")

    parser.add_option("--rename", action="store", type="string", dest="rename", metavar="FILE",
                    help="Rename leaves using information in FILE (format: old new)")

    parser.add_option("--shortnames", action="store", type="string", dest="shortnames", metavar="B[,E]",
                    help="Abbreviate names: Keep B first (and E last) characters in each name")

    parser.add_option("--nameprune", action="store", type="string", dest="nameprune", metavar="SEP[,KEEP]",
                    help="Prune clades containing multiple taxa with identical name starts (up to first SEP). KEEP: do not prune names containing this pattern")

    parser.add_option("--prune", action="store", type="string", dest="prune", metavar="FILE",
                    help="Remove leaves listed in FILE (one leaf name per line)")

    parser.add_option("--debug", action="store_true", dest="debug",
                      help="Print longer error messages")

    parser.add_option("--overwrite", action="store_true", dest="overwrite",
                      help="Overwrite identically named previous resultfiles with no warning (e.g. for clustering options)")

    parser.set_defaults(informat="newick", outformat="newick", midroot=False, rootfile=None, cladeinfo=None, shortnames=False, rename=None,
                        prune=None, nameprune=False, clustn=None, clustcut=None, clustrees=False, clustdir="clusterdir",
                        debug=False, overwrite=False)

    return parser

################################################################################################

def check_commandline(options, args, parser):

    # Set printing of tree to True. May be changed by other options
    options.printtree = True

    # Parse and check filenames
    if len(args) < 1:
        filename = ["-"]
    elif len(args) > 1:
        parser.error("More than one filename listed: {}".format(args))
    else:
        filename = args[0]
        if not os.path.isfile(filename):
            parser.error("File %s not found." % filename)

    # Check that outgroup file exists if outgroup rooting requested:
    if options.rootfile:
        if not os.path.isfile(options.rootfile):
            parser.error("Outgroup-file for rooting not found: {}".format(options.rootfile))

    # Check that rename file exists if renaming is requested:
    if options.rename:
        if not os.path.isfile(options.rename):
            parser.error("File for renaming (%s) not found." % options.rename)

    # Check that cladeinfo file exists if requested. Also turn off printing of tree:
    if options.cladeinfo:
        options.printtree = False
        if not os.path.isfile(options.cladeinfo):
            parser.error("File for cladeinfo analysis ({}) not found.".format(options.cladeinfo))

    # If options.clustrees: check that clustn or clustcut specified. Also turn off printing of tree:
    if options.clustrees and not (options.clustn or options.clustcut):
        parser.error("Option --clustrees requires either --clustn or --clustcut")

    # If cluster analysis requested: Create clustrees folder if it does not already exist, or if overwrite specified.
    # Also turn off printing of tree
    if options.clustn or options.clustcut:
        options.printtree = False
        if not os.path.exists(options.clustdir):
            os.makedirs(options.clustdir)
        elif not options.overwrite:
            parser.error("Folder for cluster subtrees ({}) already exists.".format(options.clustdir))
        else:
            shutil.rmtree(options.clustdir)
            os.makedirs(options.clustdir)

    # Check that prune file exists if pruning is requested:
    if options.prune:
        if not os.path.isfile(options.prune):
            parser.error("File for renaming (%s) not found." % options.prune)

    return (options, args, filename)

################################################################################################

def read_treefile(options, filename):
    """Takes filename as input, returns Tree object"""

    if options.informat.lower() == "nexus":
        treefile = treelib.Nexustreefile(filename)
    else:
        treefile = treelib.Newicktreefile(filename)

    tree = next(treefile)
    return tree

################################################################################################

def read_names(filename):
    """File with name "filename" assumed to contain one leafname per line. Read and return set of names"""

    names = set()
    with open(filename, "r") as infile:
        for line in infile:
            leaf = line.strip()
            if leaf:
                names.add(leaf)
    return names

################################################################################################

def analyze_tree(options, tree):

    # Midpoint root tree if requested
    if options.midroot:
        tree.rootmid()

    # Root tree on outgroup listed in rootfile if requested
    elif options.rootfile:
        outgroup = read_names(options.rootfile)
        tree.rootout(outgroup)

    # Abrreviate all leaf names
    if options.shortnames:
        numlist = options.shortnames.split(",")
        b = int(numlist[0])
        if len(numlist) == 2:
            e = int(numlist[1])
        else:
            e = 0
        leafnames = tree.leaflist()
        for oldname in leafnames:
            newname = oldname[:b] + "_" + oldname[-e:]    # First b + last e characters
            tree.rename_leaf(oldname, newname, fixdups=True)

    # Prune leafs based on name redundancy: Find subtrees where all leafs have same start of name
    if options.nameprune:
        patlist = options.nameprune.split(",")
        sep = patlist[0]
        if len(patlist) == 2:
            keep = patlist[1]
        else:
            keep=None
        tree.nameprune(sep, keep)

    # Rename leafs based on "old new" pairs in FILE
    if options.rename:
        tree.transname(options.rename)

    # Remove leaves if requested
    if options.prune:
        remleaves = read_names(options.prune)
        tree.remove_leaves(remleaves)

    # If requested: check if leaves listed in options.cladeinfo form a monophyletic group
    # Regardless: find smallest containing clade and report distance from root to basenode and basenode's parent
    # Root tree on outgroup listed in rootfile if requested
    if options.cladeinfo:
        leaves = read_names(options.cladeinfo)
        mrca = tree.findMRCA(leaves)
        if tree.remote_children(mrca) == leaves:
            print ("Leaves in file '{}' form a monophyletic group".format(options.cladeinfo))
        else:
            print ("Warning: Leaves in file '{}' do NOT form a monophyletic group".format(options.cladeinfo))
        print("Distance from root to mrca:          {}".format(tree.nodedist(mrca, tree.root)))
        print("Distance from root to mrca's parent: {}".format(tree.nodedist(tree.parent(mrca), tree.root)))

    # If requested: Split tree into clusters by cutting across tree at suitable distance from root
    # Print clusters to stdout, one leaf per line: clusterID leafname
    if options.clustn or options.clustcut:
        if options.clustn:
            clusterlist, cluster_basenodes, unclassified = tree.cluster_n(nclust=options.clustn)
        elif options.clustcut:
            clusterlist, cluster_basenodes, unclassified = tree.cluster_cut(cutoff=options.clustcut)

        # Print cluster info to tab separated file: first column is leaf naem, second is cluster ID or "UNCLASSIFIED"
        fname = os.path.join(options.clustdir, "clusterinfo.txt")
        with open(fname, "w") as outfile:
            for leaf in unclassified:
                outfile.write("{}\tUNCLASSIFIED\n".format(leaf))
            clusterno = 0
            for cluster in clusterlist:
                clusterno += 1
                for leaf in cluster:
                    outfile.write("{}\t{}\n".format(leaf, clusterno))

        # If requested: also output files for each subtree, and one for the stub of the tree showin cluster placements
        # Note: unclassified leaves are also shown on stub
        # Note 2: cluster IDs are implicit in list indices for clusterlist AND cluster_basenodes (= index + 1)
        if options.clustrees:
            if options.outformat == "nexus":
                suffix = "nexus"
            else:
                suffix = "newick"

            # Write treefiles for each cluster
            ndigits = int( math.log10( len(cluster_basenodes) ) ) + 1       # Number of digits required to write cluster IDs
            clusterno = 0
            for basenode in cluster_basenodes:
                clusterno += 1
                clusID = str(clusterno).zfill(ndigits)
                fname = os.path.join(options.clustdir, "cluster_{}.{}".format(clusID, suffix) )
                with open(fname, "w") as outfile:
                    cluster_tree = tree.subtree(basenode)
                    if options.outformat == "nexus":
                        outfile.write(cluster_tree.nexus())
                    else:
                        outfile.write(cluster_tree.newick())

            # Create treefile for remaining stub (what is left after cutting across branches)
            # Location of clusters (now cut off) is indicated by special leaves having clusterID as names
            clusterno = 0
            for basenode in cluster_basenodes:
                clusterno += 1
                clusID = "Cluster_" + str(clusterno).zfill(ndigits)
                leaflist = tree.remote_children(basenode)
                tree.collapse_clade(leaflist, clusID)
            fname = os.path.join(options.clustdir, "stub.{}".format(suffix) )
            with open(fname, "w") as outfile:
                if options.outformat == "nexus":
                    outfile.write(tree.nexus())
                else:
                    outfile.write(tree.newick())

################################################################################################

def print_tree(options, tree):
    """Accepts either Tree or Tree_set as input. Prints all trees on stdout"""

    # Print tree on standard out
    if options.outformat.lower() == "nexus":
        print(tree.nexus())
    else:
        print(tree.newick())

################################################################################################

def main():
    parser = build_parser()
    (options, args) = parser.parse_args()

    try:
        (options, args, treefilename) = check_commandline(options, args, parser)
        tree = read_treefile(options, treefilename)

        analyze_tree(options, tree)

        if options.printtree:
            print_tree(options, tree)

    except treelib.TreeError as exc:
        if options.debug:
            import traceback
            traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write("Error: %s\n" % exc.errormessage)


################################################################################################

if __name__ == "__main__":
    main()
