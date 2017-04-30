import sys

infile=sys.argv[1]
outfile=open(sys.argv[2], 'w')

for line in open(infile):
        parts = line.split('\t')
        try:
        	parts[4] == ''
        except IndexError:
        	continue
        if parts[0] == "CHROM":
        	outfile.write(line)
        elif parts[7][0].isalpha():
        	outfile.write(parts[0]+'\t'+parts[1]+'\t'+parts[2]+'\t'+parts[3]+'\t'+parts[4]+'\t'+parts[5]+'\t'+parts[6]+'\t'+parts[7]+'\t'+'NA'+'\t'+parts[8]+'\t'+'\n')
        elif parts[7][0].isdigit():
        	outfile.write(parts[0]+'\t'+parts[1]+'\t'+parts[2]+'\t'+parts[3]+'\t'+parts[4]+'\t'+parts[5]+'\t'+parts[6]+'\t'+'NA'+'\t'+parts[7]+'\t'+parts[8]+'\t'+'\n')
        else:
        	outfile.write(line)
