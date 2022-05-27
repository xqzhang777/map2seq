import os
import tempfile
import sys
from graph import parse_file
from fms import fms
import argparse

#defaults
#dir = "../Testing"
#mapFile = "spalsy1.mrc"
#pdbFile = "spalsy1f.pdb"
#reverse_pdb = "spalsy1b.pdb"
tophits = "--tophits 30000"
#db = "human.fa.gz"
hmmer_out = "hmmer_output.txt"

def main(map, pdb, db = "tempDir/human.fa.gz", dir = "tempDir", graphname1 = "fms"):

    map = os.path.abspath(map)
    pdb = os.path.abspath(pdb)
    db = os.path.abspath(db)
    dir = os.path.abspath(dir)

    if dir[-1] != "/":
        dir += "/"

    print('start')
    #create temp directory for hmmer_output.txt
    with tempfile.TemporaryDirectory() as tmpdir:
    #with tempfile.mkdtmp() as tmpdir:
        orig_stdout = sys.stdout
        f = open(os.path.join(dir, 'log.txt'), 'w')
        #sys.stdout = f
        tmpdir = os.path.abspath(tmpdir)
        tmpdir += '/'
        #run fms on forward pdb
        print(map)
        print(pdb)
        print(db)
        print(tmpdir)
        fms(map, pdb, db, tmpdir)
        #print('Making figure')
        #graph fms output
        #try:
        num = parse_file("{}{}.png".format(dir, graphname1), "{}{}".format(tmpdir, hmmer_out))
        #print('Done')
        sys.stdout = f
        if num == -1:
            print("No Matches")
        elif num == 1:
            print("Success")
        #except:
        #    print("Error")
        
        sys.stdout = orig_stdout
        f.close()
        #reverse pdb
        # revpdb = reverse(dir, pdb)
        #run fms on backward pdb
        # fms(map, revpdb, db, tmpdir)
        #graph fms output
        # parse_file("{}{}.png".format(dir, graphname2), "{}{}".format(tmpdir, hmmer_out))
        # os.system("open {}{}.png".format(dir, graphname1))
        # os.system("open {}{}.png".format(dir, graphname2))

if __name__ == "__main__":
    print('running..')
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-map", "--map", help = "Map file path", required = True)
    requiredNamed.add_argument("-pdb", "--pdb", help = "Pdb file path", required = True)
    args = parser.parse_args()
    main(args.map, args.pdb)
