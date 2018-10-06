import os
import sys
import argparse
def main():
    # Read the two BLAST output and find the queries where the top hit has changed
    parser = argparse.ArgumentParser(description="Finds queries where bitscore changes")
    parser.add_argument("-first","--first_file", help="BLAST output file 1",required=True)
    parser.add_argument("-second","--sec_file", help="BLAST output file 2", required=True)
    parser.add_argument("-p","--prefix", help="output file prefix", required = True)
    args = parser.parse_args()

    top_map = {}
    topsubject_map = {}
    blast_line = {}
    seen = {}
    with open(args.first_file) as f:
        for line in f:
            val = line.strip().split('\t')
            if val[0] not in seen:
                seen[val[0]] = 1
                top_map[val[0]] = val[-1]
                topsubject_map[val[0]] = val[1]
                blast_line[val[0]] = line.strip()

   
    fw = open(''.join([args.prefix, "_changes"]), 'w')
    seen2 = {}
    with open(args.sec_file) as f:
        for line in f:
            val = line.strip().split('\t')
            if val[0] not in seen2:
                seen2[val[0]] = 1
                
                if float(val[-1]) != float(top_map[val[0]]):
                    fw.write(''.join([line.strip(),'\n',blast_line[val[0]],'\n\n']))

if __name__ == '__main__':
    main()