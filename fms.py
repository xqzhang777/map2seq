import os
import sys

#testing
#dir = "../Map_Files/"
#mapFile = "spalsy1.mrc"
#pdbFile = "spalsy1f.pdb"
tophits = "--tophits 30000"
#db = "human.fa.gz"

def fms(mapPath, pdbPath, databasePath, tmpdirPath):
    cmd_template = "findmysequence --mapin {} --modelin {} --db {} --tmpdir {} {}".format(mapPath, pdbPath, databasePath, tmpdirPath, tophits)
    os.system(cmd_template)

def main():
    fms()

if __name__ == "__main__":
    main()