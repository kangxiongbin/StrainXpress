import sys
import re

in_file = sys.argv[1]
out_file = sys.argv[2]

i = 0
with open(in_file,'r') as fr:
    with open(out_file,'w') as fw:
        for line in fr:
            if i % 4 == 0 and re.search("/1", line):
                line.strip()
                fw.write(line[1:-3]+'\n')
            i += 1
