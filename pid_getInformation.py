################################################
# Skyler Kuhn
# Dr. Amanda Dickinson's Lab
################################################

from __future__ import division
import os


def parse(fastafile):
    sequence = ""
    header = ""
    fastaDict = {}
    with open(fastafile) as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                if sequence != "":
                    fastaDict[header] = sequence
                    sequence = ""
                header = line  # >seq1, >seq2, >seq3
            else:
                sequence += line
        else:  # we need to grab the last sequence and its annotation
            fastaDict[header] = sequence
        return fastaDict


def getlength(seq):
    return len(seq)


def translate(sequence):
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}
    aminoacid = ""
    for pos in range(0, len(sequence) - 2, 3):
        codon = sequence[pos:pos + 3]
        aminoacid += codontable[codon]
    return aminoacid


def findminseq(seq1,seq2):
    # returns a tuple, the zeroth index of the tuple is the smaller seq, the first index is the larger seq
    if len(seq1) == len(seq2):
        return (seq1,seq2)  # in this case it does not matter, they are the same length, we return seq1 as the smallest
    elif len(seq1) < len(seq2):
        return (seq1,seq2)
    elif len(seq1) > len(seq2):
        return (seq2,seq1)
    else:
        print("The sequences you are trying to compare are not strings. Please input valid data!")

def gaps_formatter(indexlist,minseq):
    gapstring = ""
    for i in range(len(minseq)):
        if i in indexlist:
            gapstring += "|"
        else:
            gapstring += " "
    return gapstring

def percentidentitiy(seq1,seq2):
    """This method will find the sequence identity of two given sequences."""
    minseq, maxseq = findminseq(seq1, seq2)  # created
    matchcount = 0
    indexlist = []
    for pos, base in enumerate(minseq):
        if base == maxseq[pos]:
            indexlist.append(pos)
            matchcount += 1
    gapstr = gaps_formatter(indexlist, minseq)
    print("-------------------------\n{}\n{}\n{}".format(minseq,gapstr,maxseq))
    return round((matchcount/len(minseq)) * 100, 3)

def getinformation():
    curdir = os.getcwd()
    queriesDict = {}  # contains the header & sequence of every file that is not dspL
    dspL_list = []
    for subdir, dirs, files in os.walk(curdir):
        for file in files:
            filepath = os.path.join(subdir, file)
            if "dspL" in filepath and "vs" not in filepath:  # this is what we are comparing everything to, we will create a seperate list
                dspL_FH = open(filepath)
                dspL_list = dspL_FH.readlines()
                headerdspL = dspL_list[0].split(" ")[0][1:]
                #print headerdspL
                dspL_list = [headerdspL , dspL_list[1].upper()]
                #print dspL_list
            elif ".fa" in filepath:
                FH = open(filepath)
                query_list = FH.readlines()
                headers = query_list[0].split(" ")[0][1:]
                queriesDict[headers] = query_list[1].upper()
    return dspL_list, queriesDict

def main():

    a = "AAGGCTT"
    b = "AAGGC"
    c = "AAGGCAT"

    try:
        print("\n\nCalculating Percent Identity of test sequences...")
        pid_AB = percentidentitiy(a,b)
        print "Percent Identity:",pid_AB

        pid_BC = percentidentitiy(b,c)
        print "Percent Identity:",pid_BC

        pid_AC = percentidentitiy(a,c)
        print "Percent Identity:",pid_AC
        print("-------------------------")
        print("Done: Process Complete without any errors!")
    except:
        print("An error occured, please review input data and proper usage!")

    # This print out hours billable
    dspL_list, queriesDict = getinformation() #returns a list containing target info and dictionary contain all the query data
    print("\n\nCalculating Percent Identity of the sequences...")
    for header, seq in queriesDict.items():
        try:
            pid_AB = percentidentitiy(seq,dspL_list[1])
            print("Comparing {} and {}".format(header,dspL_list[0]))
            print("Percent Identity: {}".format(pid_AB))
            print("-------------------------")
        except:
            print("An error occured, please review input data and proper usage!")

    print("Done: Process Complete without any errors!")

    # This print out hours billable
    hoursBillable()

def hoursBillable():
    # Hours Billable
    print("\n\n\n_____________________________________")
    print("         Hours Billable Report")
    print("_____________________________________")
    hoursworked = {'1/4/2016': 4,'1/5/2016': 4, '1/6/2016': 4, '1/9/2016': 4,
                   '1/10/2016': 4, '1/11/2016': 5, '1/12/2016': 1}   # Dictionary to keep track of hours worked!
    totalhours = 0
    for date, hour in hoursworked.items():
        print("Date: {}\t\tHours Worked: {}\n".format(date,hour))
        totalhours += hour
    print("Total Hours Worked: {}".format(totalhours))
    print("_____________________________________")


if __name__ == "__main__":

    # dspL translation and writing to fasta file
    fastadict = parse("dspL_mRNA.fa")
    print(fastadict)
    dspL_aa = translate(fastadict.values()[0].upper())
    print("{}\n!Length{}".format(dspL_aa, len(dspL_aa)))
    dspl_aaFH = open("dspL_aa.fa","w")
    dspl_aaFH.write("{}\n{}".format(fastadict.keys()[0], dspL_aa))
    dspl_aaFH.close()
    print

    # dspS translation and writing to fasta file
    fastadict1 = parse("dspS_mRNA.fa")
    print(fastadict1)
    dspS_aa = translate(fastadict1.values()[0].upper())
    print("{}\n!Length{}".format(dspS_aa, len(dspS_aa)))
    dspS_aaFH = open("dspS_aa.fa", "w")
    dspS_aaFH.write("{}\n{}".format(fastadict1.keys()[0], dspS_aa))
    dspS_aaFH.close()
    print

    # dsp_tropicalis translation and writing to fasta file
    fastadict2 = parse("dsp_tropicalis_mRNA.fa")
    print(fastadict2)
    dspT_aa = translate(fastadict2.values()[0].upper())
    print("{}\n!Length{}".format(dspT_aa, len(dspT_aa)))
    dspT_aaFH = open("dsp_tropicalis_aa.fa", "w")
    dspT_aaFH.write("{}\n{}".format(fastadict2.keys()[0], dspT_aa))
    dspT_aaFH.close()
    print

    # Parsing through the domain spreadsheet to create a nested dictionary containing the domain intverals
    domainDictionary = {}
    domainFh = open('amandaDomainfo.txt')
    headerlist = domainFh.readline().strip("\n").split("\t")[1:]  # contains column labels
    for line in domainFh:
        #print(line)
        data = line.strip("\n").split("\t")
        header2 = data[0]  # represents the row labels
        domainDictionary[header2] = {}
        domainInts = data[1:]
        for i,intv in enumerate(domainInts):
            ints = intv.split("-")
            start = int(ints[0]) - 1
            stop = int(ints[1]) - 1
            domainDictionary[header2][headerlist[i]] = (start,stop)
    print domainDictionary
    domainFh.close()

    curdir = os.getcwd()
    queriesDict = {}  # contains the header & sequence of every file that is not dspL
    dspL_list = []
    for key1 in domainDictionary:   # key1 is organism key
        #print(key1)
        for subdir, dirs, files in os.walk(curdir):
            for file in files:
                filepath = os.path.join(subdir, file)
                print(filepath)
                querypath = key1 +"_aa.fa"
                if querypath in filepath:
                    fastadictionary = parse(filepath)
                    sequence = fastadictionary.values()[0]  # there is only one sequence in this file
                    print("{}\n{}\nLength: {}".format(key1,sequence,len(sequence)))
                    for key2, interv in domainDictionary[key1].items():  # key2 is domain key
                        fastafilehandle = "./alignments/domains/"+key2+"/"+key1+"_"+key2+"_aaDomain.fa"
                        print fastafilehandle
                        fastaFH = open(fastafilehandle, "w")
                        #print "#############",key1,key2,interv
                        domainSeq = sequence[interv[0]:interv[1]]
                        #print "@@@@@@@",domainSeq
                        fastaFH.write(">{}\tdomain={}\tint={}:{}\n{}".format(key1,key2,interv[0],interv[1],domainSeq))
                        fastaFH.close()
    hoursBillable()