#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, getopt, gc, gzip

#func: How to use this script
def usage():
    print("")
    print("Usage is as follows:")
    print("")
    print("  eg:\tpython phredCheck.py -i sample.fq -n 10 &")
    print("  -i:\tSpecify the input FASTQ file(gzip format is OK).")
    print("  -n:\tSpecify the number of sequences to extract(default 10000).")
    print("  -h:\tPrint this page.")
    print("\nDetails are as follow:\n")
    for v in range(1,6):
        print("  " + scType[v])
    print("")
    sys.exit()

scType = ["Unknown data score format!",
          "Sanger | Phred+33 | Qual[33,73] | Val(0,40)",
          "Solexa | Solexa+64 | Qual[59,104] | Val(-5,40)",
          "Illumina1.3 | Phred+64 | Qual[64,104] | Val(0,40)",
          "Illumina1.5 | Phred+64 | Qual[66,104] | Val(3,40)",
          "Illumina1.8 | Phred+33 | Qual[33,74] | Val(0,41)"]

def readFile(fq_read, list, num):
    flag = 1
    for line in fq_read:
        if flag > num * 4:
            break
        elif flag % 4 == 0:
            list.append(line.strip())
        else:
            pass
        flag += 1

def judgeFile(input_fq, num, list):
    basename = str(input_fq).split('.')[-1]
    if basename == "fq" or basename == "fastq":
        with open(input_fq, 'r') as fq_read:
            readFile(fq_read, list, num)
    elif basename == "gz":
        with gzip.open(input_fq, 'rb') as fq_read:
            readFile(fq_read, list, num)
    else:
        print("Please input the FASTQ file!")
        usage()
        sys.exit()

def checkFq(input_fq, list):
    score = []
    for line in list:
        for ABC in line:
            score.extend([str(ord(ABC))])
    MIN = int(min(score))
    MAX = int(max(score))

    if MIN < 33 or MAX > 104:
        print ("%s [%s,%s] ==> %s"%(input_fq, MIN, MAX, scType[0]))
    elif MIN >= 33 and MAX <= 73:
        print ("%s [%s,%s] ==> %s"%(input_fq, MIN, MAX, scType[1]))
    elif MIN >= 59 and MAX <= 104:
        print ("%s [%s,%s] ==> %s"%(input_fq, MIN, MAX, scType[2]))
    elif MIN >= 64 and MAX <= 104:
        print ("%s [%s,%s] ==> %s"%(input_fq, MIN, MAX, scType[3]))
    elif MIN >= 66 and MAX <= 104:
        print ("%s [%s,%s] ==> %s"%(input_fq, MIN, MAX, scType[4]))
    elif MIN >= 33 and MAX <= 74:
        print ("%s [%s,%s] ==> %s"%(input_fq, MIN, MAX, scType[5]))
    else:
        print ("%s [%s,%s] ==> %s"%(input_fq, MIN, MAX, scType[0]))

# main program is here
if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], "hi:n:o:")
    input_fq = ""
    num = 10000
    list = []
    if opts:
        for op, value in opts:
            if op == "-i":
                input_fq = value
            elif op == "-n":
                num = int(value)
            elif op == "-o":
                output = value
            elif op == "-h":
                usage()
                sys.exit()
    else:
        usage()

    judgeFile(input_fq, num, list)
    checkFq(input_fq, list)
    gc.collect()
