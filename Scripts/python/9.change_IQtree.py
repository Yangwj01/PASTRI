import os
import argparse
import time
from ete3 import Tree

def Parsers():
    parser = argparse.ArgumentParser()
    #parser.add_argument("-o", "--output", type=str, help="The output iqtree file", required=True)
    parser.add_argument("-n", "--prefix", type=str, help="Give the prefix", required=True)
    parser.add_argument("-d","--directory", type=str, help="Set up an directory to run mix", required=True)
    parser.add_argument("-i","--input", type=str, help="The file of iqtree input", required=True)
    args = parser.parse_args()
    return args

def changeTree(inTree, outTree):
    ##Nexus format
    try:
        rawTree = Tree(inTree, format=1)
    except ete3.parser.newick.NewickError:
        print("Error parsing newick tree file. Please check the file format.")
        sys.exit(1)
        
    for node in rawTree.iter_descendants():
        node.dist = 1
    ## Convert Nexus format tree to Newick format
    rawTree.write(outfile=outTree, format=0)


if __name__ == "__main__":
    options = Parsers()
    start = time.time()
    ## change distance to 1
    rawTree = open(options.input,'r').read() #.treefile
    newTree = os.path.join(options.directory, options.prefix + ".nwk")
    changeTree(rawTree, newTree)
    end = time.time()
    print("{}\t{}\n".format(options.prefix,end - start))    
  


    
