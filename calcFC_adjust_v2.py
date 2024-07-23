#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 19:16:18 2018

@author: Zhanying Feng
"""

import sys
import getopt
import re
import os
import multiprocessing
import time
help_message = '''
    USAGE: calcFC -p <hotspot> -b <bam>
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def Partition(n,x,y):
    s = 1; e = int(n/30)
    Par = []
    for i in range(29):
        Par.append((s+i*e,e+i*e,x,y))
    Par.append((s+29*e,n,x,y))
    return Par

def HideFile(x):
	return x[0:x.rfind('/')+1]+'.'+x[x.rfind('/')+1:len(x)]

def ProcessPartition(arg):
    s, e, Hotfile, Bamfile = arg
    Tmpfile = HideFile(Hotfile) + '_' + str(s) + '_' + str(e) + '.tmp'
    Out = open(Tmpfile,'w')
    for i in range(s,e+1):
        CommandRead = 'sed -n ' + str(i) + 'p ' + Hotfile
        line = os.popen(CommandRead).read()
        peakLine = (line.strip('\n')).split('\t')
        chrom = peakLine[0]
        start = peakLine[1]
        end = peakLine[2]
        #zscore = peakLine[3].strip('\n')
        Region=chrom+":"+start+"-"+end 
        Command="samtools view "+Bamfile+ " "+Region +" -c"
        count=os.popen(Command).read();
        a=float(count)

        Length=int(end)-int(start)
        NewStart=int(start)-25000
        if (NewStart<0):
            NewStart=1
        Region5K=chrom+":"+str(NewStart)+"-"+str(int(end)+25000)
        Command5K="samtools view "+Bamfile+ " " + Region5K +" -c"
        count=os.popen(Command5K).read()
        b=float(count)
        if b<0.1:
            b=1

        NewStart=int(start)-500000
        if (NewStart<0):
                NewStart=1	
        Region1000K=chrom+":"+str(NewStart)+"-"+str(int(end)+500000)
        Command1000K="samtools view "+Bamfile+ " " + Region1000K +" -c"
        count=os.popen(Command1000K).read()
        c=float(count)
        if c<0.1:
            c=1

        FC1=a*(50000+Length)/Length/b
        FC2=a*(1000000+Length)/Length/c
        Out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chrom, start, end, Length, a, b, FC1, c, FC2))
    Out.close()
 

def processHot(Hotfile, Bamfile):
    """
        Load a hotspot file and calculate the fold change for each region.
    """
    PeakNum = int((os.popen('wc -l ' + Hotfile).read()).split(' ')[0])
    if PeakNum > 30:
        Par = Partition(PeakNum, Hotfile, Bamfile)
        pool = multiprocessing.Pool(30)
        pool.map(ProcessPartition,Par)
        pool.close();pool.join()
        CommandMerge = 'cat '; CommandDel = 'rm -f '
        for i in range(30):
            CommandMerge += HideFile(Hotfile) + '_' + str(Par[i][0]) + '_' + str(Par[i][1]) + '.tmp '
            CommandDel += HideFile(Hotfile) + '_' + str(Par[i][0]) + '_' + str(Par[i][1]) + '.tmp '
        CommandMerge += ('> ' + Hotfile+'.fc')
        os.popen(CommandMerge);time.sleep(180);os.popen(CommandDel);
    else:
        ProcessPartition((1, PeakNum, Hotfile, Bamfile))
	CommandMv = 'mv ' + HideFile(Hotfile) + '_' + str(1) + '_' + str(PeakNum) + '.tmp ' + Hotfile+'.fc'
	os.popen(CommandMv)

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
