#!/usr/bin/env python
# encoding: utf-8
"""
calcFC.py

Created by Yong Wang on 2015-02-06.
Copyright (c) 2015 Yong Wang All rights reserved.
"""

import sys
import getopt
import re
import os

help_message = '''
    USAGE: calcFC -p <hotspot> -b <bam>
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def processHot(Hotfile, Bamfile):
    """
        Load a hotspot file and calculate the fold change for each region.
    """
    Out=open(Hotfile+'.fc','w')
    f = open(Hotfile,'r')
    Num=0
    for line in f:
        Num=Num+1
        peakLine = line.split('\t')
        chrom = peakLine[0]
        start = peakLine[1]
        end = peakLine[2]
        zscore = peakLine[3].strip('\n')
        Region=chrom+":"+start+"-"+end 
        #Command="samtools view "+Bamfile+ " "+Region +"| wc -l"
        Command="samtools view "+Bamfile+ " "+Region +" -c"
        #print Command
        #print Command
        count=os.popen(Command).read();
        a=float(count)

        Length=int(end)-int(start)
        NewStart=int(start)-25000
        if (NewStart<0):
            NewStart=1
        Region5K=chrom+":"+str(NewStart)+"-"+str(int(end)+25000)
        #Command5K="samtools view "+Bamfile+ " " + Region5K +"| wc -l"
        Command5K="samtools view "+Bamfile+ " " + Region5K +" -c"
        count=os.popen(Command5K).read()
        #count=os.popen(Command5K).read()
        b=float(count)
        if b<0.1:
            b=1

        NewStart=int(start)-500000
        if (NewStart<0):
                NewStart=1	
        Region1000K=chrom+":"+str(NewStart)+"-"+str(int(end)+500000)
        #Command1000K="samtools view "+Bamfile+ " " + Region1000K +"| wc -l"        
        Command1000K="samtools view "+Bamfile+ " " + Region1000K +" -c"
        count=os.popen(Command1000K).read()
        #count=os.popen(Command1000K).read()
        c=float(count)
        if c<0.1:
            c=1

        FC1=a*(50000+Length)/Length/b
        FC2=a*(1000000+Length)/Length/c
        Out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chrom, start, end, zscore, Length, a, b, FC1, c, FC2))
        if (Num%1000==0):
		print Num
	        print chrom, start, end, zscore, Length, a, b, FC1, c, FC2 

    f.close()
    Out.close()
    print 'In total %s lines processed'%Num


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hp:b:", ["help", "peak", "bam"])
        except getopt.error, msg:
            raise Usage(help_message)

        #print args
        #print opts
        # option processing
        bamFile = ""
        hotFile = ""
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-b", "--bam"):
                bamFile = value
            if option in ("-p", "--peak"):
                hotFile = value

        try:
           f = open(bamFile, 'r')
        except IOError, msg:
            raise Usage(help_message)
        try:
           f = open(hotFile, 'r')
        except IOError, msg:
            raise Usage(help_message)


    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        return 2

    # input a hotspot and BAM file.
    processHot(hotFile,bamFile)

if __name__ == "__main__":
    sys.exit(main())

