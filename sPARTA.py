## sPARTA: small RNA-PARE Target Analyzer public version 
## Updated: version-1.25 05/13/2017
## Author: kakrana@udel.edu
## Author: rkweku@udel.edu

## Copyright (c): 2016, by University of Delaware
##              Contributor : Atul Kakrana
##              Affilation  : Meyers Lab (Donald Danforth Plant Science Center, St. Louis, MO)
##              License copy: Available at https://opensource.org/licenses/Artistic-2.0


#### PYTHON FUNCTIONS ##############################
import sys,os,re,time,glob,shutil,operator,datetime,argparse,importlib,getpass
import subprocess, multiprocessing
from multiprocessing import Process, Queue, Pool
from operator import itemgetter
from os.path import expanduser


#### USER SETTINGS ################################
parser = argparse.ArgumentParser()
parser.add_argument('-annoType',  default='', help='GFF if the annotation file'\
    'file is a GFF file. GTF if the annotation file is a GTF file')
parser.add_argument('-annoFile',  default='', help='GFF/GTF file for the '\
    'species being analyzed corresponding to the genome assembly being used')
parser.add_argument('-genomeFile', default='', help='Genome file in FASTA '\
    'format')
parser.add_argument('-featureFile', default='', help='Feature file in FASTA '\
    'format')
parser.add_argument('-genomeFeature', required=True, help='0 if prediction is '\
    'to be done in genic region. 1 if prediction is to be done in intergenic '\
    'region')
parser.add_argument('-miRNAFile', default='', help='FASTA format of miRNA '\
    'sequences')
parser.add_argument('-tarPred', nargs='?', const='H', help='Mode of target '\
    'prediction. H for heuristic. E for exhaustive. H is default if no mode '\
    'is specified')
parser.add_argument('-tarScore', nargs='?', const='S', help='Scoring mode '\
    'for target prediction. S for seedless. N for normal. S is default if '\
    'no mode is specified')
parser.add_argument('-libs', nargs='*', default=[], help='List of PARE '\
    'library files in tag count format. Data can be converted into tag '\
    'count format using')
parser.add_argument('-minTagLen', default=20, type=int, help='Minimum length '\
    'of PARE tag. Tags shorter than minTagLen will be discarded. 20 is '\
    'default')
parser.add_argument('-maxTagLen', default=30, type=int, help='Maximum length '\
    'of PARE tag. Tags longer than maxTagLen will be chopped to the specified '\
    'length. 30 is default')
parser.add_argument('--tag2FASTA', action='store_true', default=False, help=
    'Convert tag count file for PARE libraries to FASTA files for mapping')
parser.add_argument('--map2DD', action='store_true', default=False, help=
    'Map the PARE reads to feature set')
parser.add_argument('--validate', action='store_true', default=False, help=
    'Flag to perform the validation of the potential cleave sites from '\
    'miRferno')
parser.add_argument('--repeats', action='store_false', default=True, help=
    'Flag to include PARE reads from repetitive regions')
parser.add_argument('--noiseFilter', action='store_false', default=True,
    help='Flag to include all PARE validations with p-value of <=.5, '\
    'irrespective of the noise to signal ratio at cleave site and category '\
    'of PARE read')
parser.add_argument('-accel', default='Y', help='Y to use '\
    'balanced multiple process scheme or else specify the number of '\
    'processors to be used. Y is default')
parser.add_argument('--standardCleave', action='store_true', default=False,
    help='Flag to use standard cleave locations (10th, 11th and 12th '\
    'positions), or rules more specific to miRNA size')

### Developer Options ###
parser.add_argument('--generateFasta', action='store_false', default=False,
    help=argparse.SUPPRESS)
parser.add_argument('--fileFrag', action='store_true', default=False,
    help=argparse.SUPPRESS)
parser.add_argument('--indexStep', action='store_true', default=False,
    help=argparse.SUPPRESS)
parser.add_argument('-splitCutoff', default=20, help=argparse.SUPPRESS)
parser.add_argument('-maxHits', default=30, help=argparse.SUPPRESS)
parser.add_argument('--cat4Show', action='store_false', default=True,
    help=argparse.SUPPRESS)

args = parser.parse_args()

### Various checks for dependencies within command line arguments

# If either annotation or genome file is given without the other and
# featureFile is not given, exit.
if(((args.annoFile and not args.genomeFile) or (args.genomeFile and not
        args.annoFile)) and (not args.featureFile)):
    print("annoFile and genomeFile both must be provided to extract seqeunces")
    exit()

# If annoType is provided and not GFF or GTF, report the error and exit
if(args.annoType and args.annoType != 'GFF' and args.annoType != 'GTF'):
    print("annoType must be either GFF3 or GTF")
    exit()

# If either the annotation file or annotation type is given without the other,
# exit.
if((args.annoType and not args.annoFile) or (args.annoFile and not
        args.annoType)):
    print("annoType and annoFile must both be give to parse either the GFF "\
    "or GTF file.")
    exit()

# If the user input both a genome and feature file, exit as both cannot be
# supplied for proper execution
if(args.genomeFile and args.featureFile):
    print("genomeFile and featureFile cannot both be supplied for execution")
    exit()

# If annoFile and genomeFile are given turn on extraction, frag and index steps
# must be set on
if(args.annoFile and args.genomeFile):
    args.generateFasta = True
    args.fileFrag = True
    args.indexStep = True

# If featureFile is given, frag and index steps must be set on
if(args.featureFile):
    # If featureFile is given and annoFile is given, give a warning letting
    # user know the annoFile will be ignored and the input fasta file may
    # have been intended as a genomeFile
    if(args.annoFile):
        print("Warning: You have input a annoFile but input a FASTA file as "\
        "the featureFile. If you intended for this to be used in conjunction "\
        "with the annotation file to create a feature file, please press "\
        "'ctrl+c' to cancel the execution and rerun with the FASTA file "\
        "under the argument 'genomeFile'. If this is in fact the feature "\
        "file, allow sPARTA to continue its execution.")
        time.sleep(10)
    args.fileFrag = True
    args.indexStep = True

# If indexStep is on and tarPred is off, turn tarPred and tarScore on
if(args.indexStep):
    if(not args.tarPred):
        args.tarPred = 'H'
    if(not args.tarScore):
        args.tarScore = 'S'

# If tarPred is on, then tarScore will default to S
if(args.tarPred and not args.tarScore):
    args.tarScore = 'S'

# If tarPred is on, then miRNAFile must be provided
if(args.tarPred and not args.miRNAFile):
    print("miRNA file must be given to perform target prediction")
    exit()

# If tag2FASTA is on, turn map2DD on
if(args.tag2FASTA and not args.map2DD):
    args.map2DD = True

# If tag2FASTA is on, then libraries must be defined
if(args.tag2FASTA and not args.libs):
    print("libs must be assigned to perform tag2FASTA")
    exit()

# If validate is on, then libraries must be input
if(args.validate and not args.libs):
    print("At least one library must be given to perfor the validate")
    exit()


# genomeFeature must be an integer
args.genomeFeature = int(args.genomeFeature)

####################################################################
#### sPARTA FUNCTIONS ##############################################

def checkUser():
    '''
    Checks if user is authorized to use script
    '''
    print ("\n#### Fn: Checking user #################################")
    auser       = getpass.getuser()
    localtime   = time.asctime( time.localtime(time.time()) )
    # print ("Local current time :", localtime)
    print("Hello '%s' - %s" % (auser,localtime))
    print("You can seek help or report issues at: https://github.com/atulkakrana/sPARTA/issues" )
    # if auser in allowedUser:
    #     print("Hello '%s' - Issues need to be reproted: https://github.com/atulkakrana/phasTER/issues \n" % (auser))
    # else:
    #     print("YOU ARE NOT AUTHORIZED TO USE DEVELOPMENT VERSION OF 'sPARTA'")
    #     print("Contact 'Atul Kakrana' at kakrana@gmail.com for permission\n")
    #     sys.exit()
    
    return None

def checkLibs():
    '''Checks for required components on user system'''

    print("\n#### Fn: Check libraries and components ################")
    goSignal    = 1 
    # isRpy2      = importlib.find_loader('rpy2')
    # if isRpy2 is None:
    #     print("--rpy2   : missing")
    #     goSignal    = 0
    # else:
    #     print("--rpy2   : found")
    #     pass

    isNumpy     = importlib.find_loader('numpy')
    if isNumpy is None:
        print("--numpy              : missing")
        goSignal    = 0
    else:
        print("--numpy              : found")
        pass

    isScipy     = importlib.find_loader('scipy')
    if isScipy is None:
        print("--scipy              : missing")
        goSignal    = 0
    else:
        print("--scipy              : found")
        pass

    # isPyfasta   = importlib.find_loader('pyfasta')
    # if isPyfasta is None:
    #     print("--pyfasta: missing")
    #     goSignal    = 0
    # else:
    #     print("--pyfasta: found")
    #     pass

    if goSignal == 0:
        print("\nPlease install missing libraries before running the analyses")
        print("See README for how to install these")
        print("sPARTA has unmet dependendies and will exit for now\n")
        print("You can seek help or report issues at: https://github.com/atulkakrana/sPARTA/issues")
        sys.exit()

    # isTally = shutil.which("tally")
    # if isTally:
    #     print("Found:Tally")
    #     pass
    # else:
    #     print("Please install 'Tally' before using the tool")
    #     print("See README for how to INSTALL")
    #     sys.exit()

    return None

def checkfiles():
    '''
    pre-analysis check for files in working directory
    '''

    print("\n#### Fn: Check input files #############################")
    goSignal    = 1 

    ## 
    if args.featureFile:
        ## Check for feature file
        if not os.path.isfile(args.featureFile):
            # print("'%s' genome file not found at:%s" % (genomeFile,os.getcwd()))
            print("--featureFile        : missing")
            goSignal = 0
        else:
            print("--featureFile        : found")
            pass 

    if args.miRNAFile:
        ## check for miRNA file
        if not os.path.isfile(args.miRNAFile):
            # print("'%s' file not found at:%s" % (args.miRNAFile,os.getcwd()))
            # print("Please check if miRNA fasta file exists in your directory\n")
            print("--miRNAFile          : missing")
            goSignal = 0
        else:
            print("--miRNAFile          : found")
            pass  

    if args.annoFile:
        ## Check for annotation file
        if not os.path.isfile(args.annoFile):
            # print("'%s' file not found at:%s" % (args.annoFile,os.getcwd()))
            # print("Please check if GTF/GFF3 file exists in your directory\n")
            print("--annoFile           : missing")
            goSignal = 0
        else:
            print("--annoFile           : found")
            pass
    
    if args.genomeFile:
        ## Check for genome file
        if not os.path.isfile(args.genomeFile):
            # print("'%s' genome file not found at:%s" % (genomeFile,os.getcwd()))
            print("--genomeFile         : missing")
            goSignal = -1
        else:
            print("--genomeFile         : found")
            pass  

    if args.libs:
        ## Check for PARE libs
        for alib in args.libs:
            if not os.path.isfile(alib):
                # print("'%s' file not found at:%s" % (alib,os.getcwd()))
                # print("Please check if GTF/GFF3 file exists in your directory\n")
                print("--PARE lib(s)        : missing")
                goSignal = -2
            else:
                print("--PARE lib(s)        : found")
                pass


    if goSignal == 0:
        print("\n** Some input file(s) couldn't be located !!")
        print("** sPARTA expects all input file to be in same directory with sPARTA program")
        print("** Please copy the missing input files in '%s' directory and rerun sPARTA" % (os.getcwd()))
        print("** sPARTA will exit for now \n")
        print("You can seek help or report issues at: https://github.com/atulkakrana/sPARTA/issues")
        sys.exit()
    
    elif goSignal == -1:
        print("\n** Genome FASTA file couldn't be located !!")
        print("** sPARTA expects all input file to be in same directory with sPARTA program")
        print("** Please copy the genome FASTA file in '%s' directory and rerun sPARTA" % (os.getcwd()))
        print("** Make sure that your genome FASTA has no digits/numbers in filename, otherwise it will be deleted in cleanup operation")
        print("** sPARTA will exit for now \n")
        print("You can seek help or report issues at: https://github.com/atulkakrana/sPARTA/issues")
        sys.exit()
    
    elif goSignal == -2:
        print("\n** All PARE/Degradome libraries couldn't be located !!")
        print("** sPARTA expects all libraries in the same directory with sPARTA program")
        print("** At least one of the input library is missing from '%s' directory" % (os.getcwd()))
        print("** Please check and and rerun sPARTA")
        print("** sPARTA will exit for now \n")
        print("You can seek help or report issues at: https://github.com/atulkakrana/sPARTA/issues")
        sys.exit()

    else:
        ## Good to go !!
        pass

    return None

def genomeReader(genomeFile):
    '''
    Reads Genome FASTA file
    '''

    print("\n#### Fn: genomeReader ##################################")
    
    #### Check for genome file
    if not os.path.isfile(genomeFile):
        print("'%s' genome file not found at:%s" % (genomeFile,os.getcwd()))
        print("\nPlease check the genomeFile - Is it in sPARTA directory? Did you input wrong name?")
        print("If it exists and you input correct name then it could have been deleted while last cleanup operation - Sorry!!")
        print("If deleted then please rename genomeFile (without any integers in name) and re-copy in sPARTA folder\n")
        print("Script will exit for now\n")
        sys.exit()
    else:
        fh_in = open(genomeFile, 'r')
        pass

    print("Caching genome fasta")
    genomeFile = fh_in.read()
    genomeList = genomeFile.split('>')[1:] 
    chromoDict = {}
    chromoLenDict = {}
    for i in genomeList:
        chromoInfo  = i.partition('\n') 
        chrid       = chromoInfo[0].split()[0] 
        chrSeq      = chromoInfo[2].replace("\n", "").strip()
        # print ("head",chrSeq[:20].strip())
        # print ("excise",chrSeq[9:50].strip())
        # print ("tail",chrSeq[-20:].strip())
        # sys.exit()
        chromoDict[chrid]       = chrSeq ## Bug 05-12-17a - The value for chromosome dictionary was being stored as list here,  and list convetred to string in EXTRACTFASTA1, therefore adding [' to 5'-end and '] to 3'-end - Fixed
        chromoLenDict[chrid]    = int(len(chrSeq)) ## To convert coords to python format subtract 1nt - reverted in v1.25 as no such offset required for length

    print("Genome dict prepared for %s chromosome/scaffolds" % (len(chromoDict)))
    
    return chromoDict,chromoLenDict

def gtfParser(gtfFile):

    '''This function parses Trinity and Rocket GTF file
    to give a trascript entry and entry for all exons - Basically this is the parser for gffread geenrated files'''
    print("\n#### Fn: gtfParser #####################################")

    #### Check if file exists
    if not os.path.isfile(gtfFile):
        print("'%s' file not found at:%s" % (gtfFile,os.getcwd()))
        print("Please check if GTF file exists in your directory\n")
        sys.exit()
    else:
        pass
    
    ### Read GTF ####
    mode = 3 ## This is the mode from real module that corresponts to TopHat/Cufflink GTF files
    with open(gtfFile) as fh_in:
        lines = (line.rstrip() for line in fh_in) 
        gtfRead = list(line for line in lines if line) # Non-blank lines in a list
    fh_in.close()

    ### PARSE GTF file #####
    gtfList     = [] ## List to hold parsed GTF entries
    tempName    = [] ## Stores current trans name
    tempCoords  = [] ## Temp coords
    
    for i in gtfRead:
        # print(i)
        ent = i.split("\t")
        # print("\nEnt",ent)
        gScore  = ent[5] ## Score provided for mapping accuracy in Trinity GTF from GMAP 0 to 100
        gtype   = ent[2]
        if gtype == "exon": #### <exon/CDS> Use CDS to get similar number of transcripts as in GFF3 file - but CDS doesn't include UTRs so exon is choosen to represent protein-coding genes
            # print("\nExon:",ent)
            gchr    = ent[0]
            gstart  = int(ent[3])
            gend    = int(ent[4])
            gstrand = ent[6].translate(str.maketrans("+-","wc"))
            info    = ent[8].strip("\n").split(";")
            # print(info)
            gid     = info[1].split()[1].replace('"','') ## Gene ID
            # print('-',gid,gchr,gstart,gend,gstrand)
            
            ## Parse the info and add the exon entries  #################
            #############################################################

            # if mode == 2 or (mode == 4 and len(info)==4): ## This is  trinity GTF info
            #     ## Protein coding gene with a version number
            #     tid         = info[0].split()[1].replace('"','').split(".")[0] ## Transcript ID
            #     gid         = info[2].split()[1].replace('"','').rpartition("_")[0] ## Gene ID
            #     # aflag   = 'T' ## Trinity
            #     print('-',gid,gchr,gstart,gend,gstrand,geneType)
            #     gtfList.append((gid,tid,gchr,gstart,gend,gstrand,geneType))

            # elif mode == 3 or (mode ==4 and len(info) >= 7): ## This is rocket GTF info, with and w/o p_id "P1" in the end
            #     tid         = info[0].split()[1].replace('"','') ## Transcript ID
            #     gid         = info[1].split()[1].replace('"','') ## Gene ID
            #     # aflag   = 'R' ## Rocket
            #     geneType    = 'gene'
            #     # print('-',gid,gchr,gstart,gend,gstrand,geneType)
            #     gtfList.append((gchr,gstrand,gid,gstart,gend,geneType)) ###@@@@@@@@@@@@@@@@@@@ CHECK @@@@@@@@@@@@@@@@@@

            # else:
            #     print("-This entry has more than expected info")
            #     print("-Info block",info)
            #     print("-Info block len:%s" % len(info))
            #     print("-Debug!!")
            #     sys.exit()

            ## Check trascript change #########
            ###################################

            ## Check if transcript entry needs to be added i.e. all the exons from previous one has been captured
            if not tempName or tempName[0] == gid:
                tempName.append(gid)
                tempCoords.append(gstart)
                tempCoords.append(gend)
                tempStrand          = gstrand   ## Strand of last transcript
                tempgid             = gid       ## Gene name of last transcript
                tempChr             = gchr      ## Chr of last transcript
            
            elif tempName[0] != gid: ## Check if the new transcript is being read. if yes then add transcript entry for old exons using tempName and tempCoords
                # print("-New transcript read - summarizing transcript entry for earlier one")
                gstartProcess       = min(tempCoords)
                gendProcess         = max(tempCoords)
                gidProcess          = tempName[0]
                geneType            = 'gene'
                # print('-',tempChr,tempStrand,gidProcess,gstartProcess,gendProcess,geneType)
                gtfList.append((tempChr,tempStrand,gidProcess,int(gstartProcess),int(gendProcess),geneType)) ### @@@@@@@@@@@@@@@@@ CHECK @@@@@@@@@@@@@@@@@@@@@@ gene start and gene end
                
                ## Empty lists and fill with current exon from new transcript
                tempName            = []            ## Empty gene names
                tempCoords          = []            ## Empty exon coords from earlier gene
                tempName.append(gid)        ## Add new gene name 
                tempCoords.append(gstart)   ## Add exons coords for new gene
                tempCoords.append(gend)     ## Add exons coords for new gene
                tempStrand          = gstrand       ## Strand of last transcript
                tempgid             = gid           ## Gene name of last transcript
                tempChr             = gchr          ## Chr of last transcript
                # sys.exit()
            else:
                print("-Unforseen scenario encountered")
                print("There seems to be a problem with GTF file - Please check integrity of your GTF file")
                sys.exit()

        else:
            # print("We don't need this info for current script") ## CDS or GENE
            pass

    if len(gtfList) == 0:
        print("Check if the feature used for extraction i.e. gene, transcript, exon is correct")
        print("Debug!!")
        sys.exit()
    else:
        # print("Total entries fetched from GTF file:%s" % (str(len(gtfList))))
        # print("First 10 entries of gtfList list: %s" % (gtfList[:10]))
        pass

    ### Sort entries to prepare for extraction
    genome_info = list(gtfList)
    # genome_info_sorted = sorted(gtfList, key=operator.itemgetter(0,1,3)) ## CHR_ID, strand and start 
    genome_info.sort(key=lambda k:(k[0],k[1],k[3])) #
    genome_info_inter = genome_info
    print("Entries in genome_info:%s" % (len(genome_info)))
    print("First 10 entries of genome_info list: %s" % (genome_info[:10]))
    # sys.exit()
    
    return genome_info,genome_info_inter

def gffParser(gffFile):
    '''
    Extract coordinates from GFF file. tested with GFF3 file from phytozome
    '''
    print("\n####Fn: gffParser ######################################")
    
    #### Check for GFF file
    if not os.path.isfile(gffFile):
        print("'%s' file not found at:%s" % (gffFile,os.getcwd()))
        print("Please check if GFF file exists in your directory\n")
        sys.exit()
    else:
        pass

    #### Parse
    fh_in       = open(gffFile,'r')
    # fh_in.readline() ## GFF version - Fixes bug 05-12-17b - First entry was being missed, any line with # is skipped below
    gffRead     = fh_in.readlines()
    genome_info = []
    chromo_dict = {}
    
    for i in gffRead:
        if i.rstrip(): ## If line not empty
            ent = i.strip('\n').split('\t')
            # print (ent)
            if i.startswith("#"): ## Fixes Bug: 010115
                pass
            else:
                if ent[2] == 'gene':
                    chrID       = ent[0]
                    strand      = ent[6].translate(str.maketrans("+-","wc"))
                    geneName    = ent[8].split(';')[0].split('=')[1]
                    geneType    = 'gene'
                    # print(chrID,strand,geneName,ent[3],ent[4],geneType)
                    genome_info.append((chrID,strand,geneName,int(ent[3]),int(ent[4]),geneType))
        
    #genome_info_inter = genome_info ##
    # genome_info_sorted = sorted(genome_info, key=operator.itemgetter(0,1,3)) ## CHR_ID, strand and start 
    genome_info.sort(key=lambda k:(k[0],k[1],k[3])) #
    genome_info_inter = genome_info
    print("Entries in genome_info:%s" % (len(genome_info)))
    # print("First 10 entries of genome_info list: %s" % (genome_info[:10]))
    # sys.exit()

    fh_in.close()
    return genome_info,genome_info_inter

def extractFeatures(genomeFile,chromoDict,chromoLenDict,genome_info,genome_info_inter):
    '''
    extract coordinates of genes and intergenic regions 
    '''
    print("\n#### Fn: extractFeatures ###############################")
    alist = []##
    for i in range(0, int(len(genome_info))+1): #
        #print (i)
        gene1       = (genome_info[i])
        gene2       = (genome_info[i+1])
        gene_type   = 'inter' #
        achr        = gene1[0]
        bchr        = gene2[0]
        # print(gene1,gene2)
        # sys.exit()

        
        if gene1[3] == gene2[3] and gene1[4] == gene2[4]:
            ## Gene is same/overlapping consider next gene
            pass
        else:
            if tuple(gene1[0:2]) not in alist:#
                print ('--Caching gene coords for chromosome: %s and strand: %s' % (gene1[0], gene1[1]))
                alist.append((gene1[0:2]))
                inter_start1    = 1
                inter_end1      = int(gene1[3])-1#
                
                ## If both the genes belong to same chr and strand i.e. chromosome has atleast two genes
                if gene1[0]         == gene2[0] and gene1[1] == gene2[1]: 
                    inter_start2    = int(gene1[4])+1 ## From end of first gene of chromosome
                    inter_end2      = int(gene2[3])-1 ## Till start of second gene
                    
                    if gene1[1]     == 'w': ##The gene is on positive strand so upstream
                        inter_name1 = ('%s_up'      % (gene1[2]))
                        inter_name2 = ('%s_up'      % (gene2[2]))
                    else: ##Its on negative strand
                        inter_name1 = ('%s_down'    % (gene1[2]))
                        inter_name2 = ('%s_up'      % (gene1[2]))

                    genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))
                    genome_info_inter.append((gene2[0],gene2[1],inter_name2,inter_start2,inter_end2,gene_type))

                ## If the first two genes are not from same strand i.e. there is only one gene on this chromosme and strand
                else: ## Intergenic from end of chromosome/scaffold
                    inter_start2    = int(gene1[4])+1 ## From end of first gene of chromosome
                    # inter_end2      = '-' ### Till end of chromosome

                    try:
                        inter_end2      = int(chromoLenDict[achr]) ## Till end of chromosome
                    except KeyError:
                        print("** '%s' is either missing from FASTA file or" % (achr))
                        print("** the name does not matches with the GFF/GTF file")
                        print("**Check if the '%s' exists in your genome FASTA file" % (achr))
                        print("**If '%s' is present in your genome FASTA file then make sure the names matches\n" % (achr))
                        sys.exit()

                    if gene1[1] == 'w': ##The gene is on positive strand so upstream
                        inter_name1 = ('%s_up'      % (gene1[2]))
                        inter_name2 = ('%s_down'    % (gene1[2]))
                    else: ##Its on negative strand
                        inter_name1 = ('%s_down'    % (gene1[2]))
                        inter_name2 = ('%s_up'      % (gene1[2]))                    
                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))
                    genome_info_inter.append((gene1[0],gene1[1],inter_name2,inter_start2,inter_end2,gene_type))

            else:
                if gene1[0] == gene2[0] and gene1[1] == gene2[1]: ### If chr_id and strands are equal than find intergenic.
                    inter_start = int(gene1[4])+1#
                    inter_end   = int(gene2[3])-1 #
                    if gene2[1] == 'w': #
                        inter_name = ('%s_up' % (gene2[2]))
                    else:
                        ## reverse strand
                        inter_name = ('%s_up' % (gene1[2]))
                    genome_info_inter.append((gene2[0],gene2[1],inter_name,inter_start,inter_end,gene_type))
                
                else: #
                    inter_start = int(gene1[4])+1#
                    # inter_end   = '-' ## Different from MPPP as no table ro query for end of chromosome in public version
                    try:
                        inter_end2  = int(chromoLenDict[achr]) ## Till end of chromosome
                    except KeyError:
                        print("** '%s' is either missing from FASTA file or" % (achr))
                        print("** the name does not matches with the GFF/GTF file")
                        print("**Check if the '%s' exists in your genome FASTA file" % (achr))
                        print("**If '%s' is present in your genome FASTA file then make sure the names matches\n" % (achr))
                        sys.exit()
                    
                    if gene1[1] == 'w':
                        inter_name = ('%s_down' % (gene1[2]))
                    else: 
                        inter_name = ('%s_up' % (gene1[2]))                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name,inter_start,inter_end,gene_type)) ##Chr_id, strand
    

    ## Additional check for scaffolded genomes, if there are no genes in a scaffold it's whole seqeunce will be fetched as intergenic
    if args.genomeFeature == 1:
        for i in chromoDict.keys():
            alen = len(chromoDict[i])
            # print("Chr:%s | Length:%s" % (i,alen))
            if tuple((i,'c')) in alist:
                # print("Found")
                # sys.exit()
                pass
            else:
                # print("Get the tuple")
                inter_name = ('chr%s_c_all' % (i))
                genome_info_inter.append((i,'c',inter_name,1,alen,'inter')) ## Chr_id, strand, name, start, stop, length

            if tuple((i,'w')) in alist:
                # print("Found")
                # sys.exit()
                pass
            else:
                # print("Get the tuple")
                inter_name = ('chr%s_w_all' % (i))
                genome_info_inter.append((i,'w',inter_name,1,alen,'inter')) ## Chr_id, strand, name, start, stop, length

    
    ## Sort the list after adding intergenic regions on on basis of chr_id and strand that is essential while caching chromosme during slicing sequences
    genome_info_inter_sort = sorted(genome_info_inter, key=operator.itemgetter(0,1))
        
    gene_coords_file    = './coords'
    coords_out          = open(gene_coords_file, 'w')
    coords              = [] ## List to hold extracted coords        
    
    if args.genomeFeature == 2: ## Both gene and inter
        for ent in genome_info_inter_sort:
            print(ent)
            if ent[4] == '-': ## End of chromosome
                coords.append(ent[0:])
                coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
                
            elif int(ent[4])-int(ent[3]) > 25:
                coords.append(ent[0:])
                coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
    
    else:
        if args.genomeFeature == 0:
            genomeFilter = 'gene'
        elif args.genomeFeature == 1:
            genomeFilter = 'inter'
        else:
            pass

        for ent in genome_info_inter_sort:
            if (ent[5] == genomeFilter):
                #print(ent)
                if ent[4] == '-': ## End of chromosome
                    coords.append(ent[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
                    
                elif int(ent[4])-int(ent[3]) > 25:
                    coords.append(ent[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))

                else:
                    ## Too short to use for analysis
                    pass
    
    print ("Number of coords in 'coords' list: %s" % (len(coords)))
    coords_out.close()

    return coords

def getFASTA1(genomeFile,coords,chromoDict):

    '''
    Extracts Genes or intergenic regions based on coords list from genome seqeunce - New proposed name get features
    '''
    print("\n#### Fn: getFASTA ######################################")
    fastaOut    = './genomic_seq.fa'
    fh_out      = open(fastaOut, 'w')

    print("Fetching genes or intergenic regions")
    fastaList   = [] ## Stores name and seq for fastFile
    chromo_mem  = []
    for i in coords: ## Coords is list from annotation parser
        #print (i)
        chr_id  = i[0]
        strand  = i[1]
        gene    = i[2]
        start   = i[3] - 1 
        end     = i[4]
        # print('\nGene:%s | Start:%s End:%s Chr:%s Strand:%s' % (gene,start,end,chr_id,strand))
        
        if tuple(i[0:2]) not in chromo_mem: 
            chromo_mem.append(tuple(i[0:2]))   ## First entry of chromosome
            print("--Reading chromosome: %s and strand: %s" % (i[0],i[1]) )
            # print('--Gene:%s | Start:%s End:%s Chr:%s Strand:%s' % (gene,start,end,chr_id,strand))
            # print("-- Fetching gene:%s" % (gene))
            chrKey = i[0]
            
            try:
                chromo = str(chromoDict[chrKey]) ## Bug 05-12-17a - The value for chromosome dictionary was a list, being convetred to string, therfore adding [' to 5'-end and '] to 3'-end - Fixed
                # print ("head",chromo[:20].strip())
                # print ("excise",chromo[9:50].strip())
                # print ("tail",chromo[-20:].strip())
            except KeyError:
                print("** '%s' is either missing from FASTA file or" % (achr))
                print("** the name does not matches with the GFF/GTF file")
                print("**Check if the '%s' exists in your genome FASTA file" % (achr))
                print("**If '%s' is present in your genome FASTA file then make sure the names matches\n" % (achr))
                sys.exit()
            
            if end == '-':
                gene_seq = chromo[start:].translate(str.maketrans("autgcn","AUTGCN")) ## Usually first entry (recorded in this loop) will be from chr-start to gene-start, but if this coord/region is few nts (<20nt) 
                                          ## it is skipped, and you will encouter directly the gene-end to chr-end region, and you need this loop - Very rare case
            else:
                gene_seq = chromo[start:end].translate(str.maketrans("autgcn","AUTGCN"))

            ncount = gene_seq.count('N') ### Check of 'Ns' that cause bowtie to hang for a very long time. This is observed in chickpea genome.
            if ncount < len(gene_seq):
                if strand == 'c':
                    gene_seq_rev = gene_seq[::-1].translate(str.maketrans("TAGC","ATCG"))
                    fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                    fastaList.append((gene,gene_seq_rev))
                else:
                    fh_out.write('>%s\n%s\n' % (gene,gene_seq))
                    fastaList.append((gene,gene_seq))
    
        elif end == '-': 
            # print('--Gene:%s | Start:%s End:%s Chr:%s Strand:%s' % (gene,start,end,chr_id,strand))
            # print("--Fetching gene:%s" % (gene))
            gene_seq = chromo[start:].translate(str.maketrans("autgcn","AUTGCN")) ##
            
            ncount = gene_seq.count('N') ### Check of 'Ns' that cause bowtie to hang for a very long time. This is observed in chickpea genome.
            if ncount < len(gene_seq):
                if strand == 'c':
                    gene_seq_rev = gene_seq[::-1].translate(str.maketrans("TAGC","ATCG"))
                    fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                    fastaList.append((gene,gene_seq_rev))
                    
                else:
                    fh_out.write('>%s\n%s\n' % (gene,gene_seq))
                    fastaList.append((gene,gene_seq))
            
        else:
            # print('--Gene:%s | Start:%s End:%s Chr:%s Strand:%s' % (gene,start,end,chr_id,strand))
            # print("--Fetching gene:%s" % (gene))
            gene_seq = chromo[start:end].translate(str.maketrans("autgcn","AUTGCN")) ##
            
            ncount = gene_seq.count('N') ### Check of 'Ns' that cause bowtie to hang for a very long time. This is observed in chickpea genome.
            if ncount < len(gene_seq):
                if strand == 'c':
                    gene_seq_rev = gene_seq[::-1].translate(str.maketrans("TAGC","ATCG"))
                    fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                    fastaList.append((gene,gene_seq_rev))
                    
                else:
                    fh_out.write('>%s\n%s\n' % (gene,gene_seq))
                    fastaList.append((gene,gene_seq))
    
        # ## Test extracted coords - IDs correspond to ChickPea genome
        # if gene == "Ca_99998":
        #     print('\nGene:%s | Start:%s End:%s Chr:%s Strand:%s' % (gene,start,end,chr_id,strand))
        #     print(gene_seq_rev)
        # elif gene == "Ca_99999":
        #     print('\nGene:%s | Start:%s End:%s Chr:%s Strand:%s' % (gene,start,end,chr_id,strand))
        #     print(gene,gene_seq)
        #     sys.exit()
        # else:
        #     pass


    time.sleep(3)
    fh_out.close()
    # sys.exit()
    
    return fastaOut,fastaList

def fastaReader(fastaFile):
    
    '''Cleans FASTA file - multi-line fasta to single line, header clean, empty lines removal'''

    print("\nFn - fastaReader")
    ## Read seqeunce file
    print ('+Reading "%s" FASTA file' % (fastaFile))
    fh_in       = open(fastaFile, 'r')
    fasta       = fh_in.read()
    fasta_splt  = fasta.split('>')
    acount      = 0 ## count the number of entries
    empty_count = 0

    fastaList = [] ## Stores name and seq for fastFile
    # fastaDict = {}

    acount +=1
    for i in fasta_splt[1:]:
        acount  +=1
        ent     = i.split('\n')
        name    = ent[0].split()[0].strip()
        seq     = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
        # alen    = len(seq)
        if seq:
            fastaList.append((name,seq))
            # fastaDict[name] = seq

    print("--Total entries in phased fastaFile:%s" % (str(acount)))
    print("--fastaList generated with %s entries\n" % (str(len(fastaList)))) ## Does not counts an empty enry from split

    return fastaList

def fragmentor(FASTA,fastaList,nseq,nfrags):
    '''
    fragments and writes frags
    '''
    aslice      = int(nseq/nfrags) ## not round to float, float added to remainder
    remainder   = nseq%nfrags
    print("--Number of seqeunces:%s | required fragments:%s | computed splice:%s | remainder seqeunces:%s" % (nseq,nfrags,aslice,remainder))

    sliceStart      = 0 ## Keep track of slice end
    sliceEnd        = aslice
    acount          = 0
    for x in range (int(nfrags+1)): ## Extra one for remainder
        if x < 10:
            afile   = "%s.frag.0%s.fa" % (FASTA.rpartition('.')[0],x)
            fh_out  = open(afile,'w')
        else:
            afile   = "%s.frag.%s.fa" % (FASTA.rpartition('.')[0],x)
            fh_out  = open(afile,'w')

        if x < nfrags: ## '<' because the 'x' starts from zero and goes till nfrags-1. if nfrags = 30, x will go from 0 to 29. The last bin (nfrag+1) will have value 30 and is captured below
            ## Slice list
            print("--Fragment:%s | sliceStart:%s | sliceEnd:%s" % (acount,sliceStart,sliceEnd))
            alist       = fastaList[sliceStart:sliceEnd]
            #write fasta
            for ent in alist:
                name    = ent[0].split()[0].strip()
                seq     = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
                fh_out.write('>%s\n%s\n' % (name,seq))

            ## Slice update
            acount      += 1
            sliceStart  += aslice
            sliceEnd    += aslice

        else: ## Write remaineder seqeunces in last file
            print("--Fragment:%s | sliceStart:%s | sliceEnd: Till end" % (acount,sliceStart))
            alist       = fastaList[sliceStart:]
            #write fasta
            for ent in alist:
                name    = ent[0].split()[0].strip()
                seq     = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
                fh_out.write('>%s\n%s\n' % (name,seq))
                acount      += 1

        fh_out.close()

    return None

def fragFASTA(FASTA,fastaList):

    print('Fn - fragFASTA')
    ## Purge fragmented files from earlier run
    shutil.rmtree('./genome', ignore_errors=True)
    os.mkdir('./genome')
    
    pattern = ".*\.frag\.[0-9].*\.fa" #
    print ("+Purging older files")
    for file in os.listdir():
        if re.search(pattern,file):
            print ('--%s is being deleted' % (file))
            os.remove(file)
    
    ## Compute nfrags and seq. per frag user FASTA file ######
    statInfo = os.stat(FASTA)
    filesize = round(statInfo.st_size/1048576,2)
    print('\n+Input FASTA size: %sMB**' % (filesize))##
    
    if filesize <= args.splitCutoff: ## 
        fls = []
        fls.append(FASTA)
        print ('--No fragmentation performed for file %s' % (fls))
        
    else:
        # fh_in       = open(FASTA, 'r')
        # seq_count   = fh_in.read().count('>')
        nseq        = len(fastaList)
        print('--Number of headers in file: %s'% (nseq))
        #if genome == 'N':
        if nseq >= 30: #
            if filesize     <= 3072:
                nfrags    = int(args.maxHits)
            elif filesize   > 3072 and filesize <= 5120:
                nfrags    = int(round(args.maxHits*1.25))
            else:
                nfrags    = int(round(args.maxHits*1.5))        
            
            print ("--Size based fragmentation in process for '%s' file" % (FASTA))

            # retcode = subprocess.call(["pyfasta","split", "-n", nfrags, FASTA]) - Deprecated

            fragmentor(FASTA,fastaList,nseq,nfrags)
            fls = [file for file in os.listdir() if re.search(r'.*\.frag\.[0-9].*\.fa', file)] ## file list using regex
            #fls = glob.glob(r'%s.[0-9]{1-3}.fa' % (FASTA.split('.')[0])) ## fragmented file list ##
            #print ('The fragments: %s' % (fls))
               
        
        else:
            nfrags = str(nseq) ## 
            if fragFasta == 'Y':
                print ("--Header based fragmentation in process for '%s' file" % (FASTA))
                # retcode = subprocess.call(["pyfasta","split", "-n", nfrags, FASTA]) - Deprecated
                
                fragmentor(FASTA,fastaList,nseq,frags)
                fls = [file for file in os.listdir() if re.search(r'.*\.frag\.[0-9].*\.fa', file)]

            
            #print ('The fragments: %s' % (fls))
    
    memFile = open('frag.mem','w')
    memFile.write("fasta=%s\n" % (FASTA))
    memFile.write("genomeFeature=%s\n" % (str(args.genomeFeature)))
    memFile.write("size=%sMB\n" % (filesize))
    memFile.write("frags=%s" % (','.join(fls)))
    memFile.close()
             
    return fls

def miRinput():
    '''
    Cleans miRNA fasta file and processes for target prediction
    '''

    print("\nFn: miRNAProcessor#########################################")

    #### Check for file
    if not os.path.isfile(args.miRNAFile):
        print("'%s' file not found at:%s" % (args.miRNAFile,os.getcwd()))
        print("Please check if miRNA fasta file exists in your directory\n")
        sys.exit()
    else:
        pass    

    mirL                = [] ## List to store miRNAs
    outfile             = ('%s_clean_revcomp.fa' % (args.miRNAFile))
    miRNA_file_clean    = cleanHeader(args.miRNAFile,outfile)#
    fh_miRNA            = open(miRNA_file_clean, 'r')
    fh_out2             = open('miRinput_RevComp.fa', 'w')
    mir_base            = fh_miRNA.read()
    mir_blocks          = mir_base.split('>')
    for i in mir_blocks[1:]:
        #print (i)
        block       = i.strip('\n')
        ent         = block.split('\n')
        miRname     = ent[0]
        miRseq_rc   = ent[1].translate(str.maketrans("AUTGC","TAACG"))[::-1]

        mirL.append((miRname,'None','None',miRseq_rc))## 
        fh_out2.write('>%s\n%s\n' % (miRname,miRseq_rc))
    
    mirTable = 'None'#
    print('Total number of miRNAs in given file: %s\n' % (len(mirL)))
    
    fh_miRNA.close()
    fh_out2.close()
        
    #for i in miRs:
    #    print (i)
        
    return mirL ## miRs holds the list of miRNA name and query where as miRtable holds flag -table name or local

## New version added - Apr1/15
def tarFind4(frag):
    
    file_out = './predicted/%s.targ' % (frag.rpartition('.')[0]) ## Result File
    
    ### Make or select index
    index = "./index/%s_index" % (frag)
    if args.indexStep:
        print('**Creating index of cDNA/genomic sequences:%s\n**' % (index))
        retcode = subprocess.call(["bowtie2-build", frag, index])

    else: #
        if os.path.isfile('%s.1.bt2' % index): #
            retcode = 0
            print('**Found index of cDNA/genomic sequences:%s\n**' % (index))
        else:
            print('**Could not find index of cDNA/genomic sequences:%s\n**' % (index))
            sys.exit()
            
    if retcode == 0: ### Index creation sucessful or index already exists
        print ('Predicting targets for frag:%s using index:%s' % (frag,index))
        nspread2 = str(nspread)
        if args.tarPred == 'H': ## Heurustic
            print ("You chose 'Heuristic mode' for target identification")
            intervalFunc    = str("L,4,0.1")
            minScoreFunc    = str("L,-24,-0.5") ###~34.5 - 35 i.e 34 - Stable v1.08 
            refGap          = str("10,8") ## Bulge in miR  + 3MM, 2 seprate or consequite bulges filtered out later - - updated-Mar-23-15 - Can be (10,12) to gain speed - it will effect 1 gap and MM but make 2 gaps impossible
            readGap         = str("10,6") ## Bulge in tar + 3MM, 2bulge in tar + 2MM and no bulge 6MM updated-Mar-23-15,
            misPen          = str("6,2") ## Mismatch penalty, w/o phred score mx is used i.e. first one - can be changed to (5,2) to improve sensistivity - It will only effect 1 gap scenario and 2 gaps cases will be as is
            matScore        = str("0") ## Match score
            retcode2        = subprocess.call(["bowtie2","-a","--end-to-end","-D 3","-R 2","-N 1","-L 8","-i",intervalFunc,"--rdg",readGap,"--rfg",refGap,"--mp",misPen,"--ma",matScore,"--min-score",minScoreFunc,"--norc","--no-hd","--no-unal","-p",nspread2, "-x", index, "-f" ,"miRinput_RevComp.fa","-S", file_out])
        
        elif args.tarPred == 'E': ##Exhaustive
            print ("You chose 'Exhaustive mode' for target identification - Please be patient")
            intervalFunc    = str("L,2,0.1")
            minScoreFunc    = str("L,-24,-0.5") ### ~34.5 - 35 i.e 34 - Stable v1.08                      
            refGap          = str("10,8") ## Bulge in miR  + 3MM, 2 seprate or consequite bulges filtered out later - updated-Mar-23-15
            readGap         = str("10,4") ## Bulge in tar + 4MM, 2bulge in tar + 3MM and no bulge 6MM - updated-Mar-23-15,
            misPen          = str("5,2") ## Mismatch penalty, w/o phred score mx is used i.e. first one
            matScore        = str("0") ## Match score
            retcode2        = subprocess.call(["bowtie2","-a","--end-to-end","-D 4","-R 2","-N 1","-L 6","-i",intervalFunc,"--rdg",readGap,"--rfg",refGap,"--mp",misPen,"--ma",matScore,"--min-score",minScoreFunc,"--norc","--no-hd","--no-unal","-p",nspread2, "-x", index,"-f", "miRinput_RevComp.fa","-S", file_out])
        
        else:
            print ('''\nPlease choose correct target prediction mode - Heuristic (H) or Exhaustive (E)\n
                   miRferno will exit now''')
            sys.exit()
    else:
        print("There seems to be a problem with index generation of locating them - miRferno will exit now")
        print("Try reruning the analysis with all steps to generate fresh files")
        sys.exit()
    
    ### Check for proper completion of Target prediction
    if retcode2 == 0:## The bowtie mapping exit with status 0, all is well
                    print('\n miRNAs mapped to Fragment: %s' % (frag))
    else:
        print ("There is some problem with miRNA mapping '%s' to cDNA/genomic seq index" % (frag))
        print ("Script exiting.......")
        sys.exit()

## Deprecated - Apr-1 [Retained for backward compatibility]
def tarParse3(targComb):
    
    print ('\n**Target prediction results are being generated**')
    #
    print ("File for parsing: '%s' in predicted folder\n" % (targComb))
    fh_in = open(targComb,'r')
    TarPred =  './predicted/%s.parsed.csv' % (targComb.rpartition('/')[-1]) ### Similar to parsed target finder format
    fh_out = open(TarPred,'w')
    fh_out.write('miRname,Target,BindSite,miRseq,tarSeq,Score,Mismatch,CIGAR\n')
    
    
    acount = 0 #
    parseCount = 0 #
    for i in fh_in:
        acount += 1
        ent = i.strip('\n').split('\t')
        #print('\n%s\n' % ent)
        miRrevcomp = ent[9] #
        miRrev = miRrevcomp.translate(str.maketrans("TACG","AUGC")) # 
        tarHash = list(miRrevcomp) #
        
        gapinfo = ent[5]
        gappos = re.split("[A-Z]",gapinfo) ## 
        gapNuc = re.findall("[A-Z]",gapinfo)
        posCount = 0
        for x,y in zip(gappos[:-1],gapNuc):## I
            #print(x,y)
            if y == 'I':
                #tarHash.insert(posCount,'-') ## 
                tarHash[posCount] = '-' ###OK
                posCount += int(x)
            else:
                posCount += int(x)       

        misinfo = ent[-2].split(':')[-1] #
        mispos = re.split("[A,T,G,C,N]",misinfo) #
        misposCorrect = [int(x)+1 for x in mispos] # 
        misNuc = re.findall("[A,T,G,C,N]",misinfo) #
        posCount = 0
        for x,y in zip(misposCorrect,misNuc):
            #print(x,y)
            posCount += x #
            gaps = tarHash[:posCount-1].count('-') #
            #print ('Position of mismatch:%s' % (posCount))
            tarHash[posCount-1+gaps] = y

        tar = ''.join(tarHash).replace("T","U") ##
        bindsite = '%s-%d' % (ent[3],int(ent[3])+(len(miRrev)-1))

        
        ##
        gap = [] #
        mis = [] #
        wobble = [] #
        nt_cnt = 1 #

        for x,y in zip(miRrevcomp[::-1].replace("T","U"),tar[::-1]):#
            if x == '-' or y == '-':
                #print('gap')
                gap.append(nt_cnt)
                if y == '-':
                    nt_cnt+=1
                
            elif x == 'A' and y == 'G': #
                #print ('wobble')
                wobble.append(nt_cnt)
                nt_cnt+=1
            elif x == 'C' and y == 'U': #
                #print ('wobble')
                wobble.append(nt_cnt)
                nt_cnt+=1
            elif x == y:
                #print('match')
                nt_cnt+=1
            else:
                #print('mismatch')
                mis.append(nt_cnt)
                nt_cnt+=1
                

        score = 0   ## Initialize
        #print (mis)
        
        if args.tarScore == 'S': #
            mis2 = list(mis)
            #if set([10,11]).issubset(mis): 
            if 10 in mis and 11 in mis: 
                score += 2.5
                #print('Removing 10')
                mis2.remove(10)
                #print ('Removing 11')
                mis2.remove(11) #
                
            for i in mis2:
                    score += 1
            for i in gap:
                score += 1.5
            for i in wobble:
                if (i+1 in mis) or (i-1 in mis): #
                    score += 1.5
                elif (i+1) in mis and (i-1 in mis):
                    score += 2
                else:
                    score += 0.5
        else:
            ##Heuristic and Exhaustive
            for i in mis:
                if i>= 2 and i<=13:
                    score += 2
                else:
                    score += 1
            for i in gap:
                if i>= 2 and i<=13:
                    score += 2
                else:
                    score += 1
            for i in wobble:
                if i>= 2 and i<=13:
                    score += 1
                    #print ('Wobble pos:%s' % (i))
                else:
                    score += 0.5
                    #print ('Wobble pos:%s' % (i))
        ###################
            
        #print(ent[0],ent[2],bindsite,miRrev,tar,score,misinfo,gapinfo)## MiRname, Tarname, mirSeq,Taerseq,binding site
        fh_out.write('>%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0],ent[2],bindsite,miRrev,tar,score,misinfo,gapinfo))
        parseCount  += 1
    
    
    print("Total number of interactions from 'miRferno':%s AND total interactions scored: %s" % (acount,parseCount))
    fh_in.close()
    fh_out.close()

    return TarPred

#### New Version added - Apr1/15
def tarParse4(targComb):

    ''' Modifying this function is worst nightmare of life - Needs cleaning
    cutoffs w/o bulge or gap - 5MM + 1 wobble 
    bulge in miRNA - 1 bulge + 3MM
    bulge in reference - 1bulge+4mm or 2 consequite bulges +3MM'''
    
    print ('\n**Target prediction results are being generated**')
    ## Input / Output file ######
    print ("File for parsing: '%s' in predicted folder\n" % (targComb))
    fh_in   = open(targComb,'r')
    TarPred =  './predicted/%s.parsed.csv' % (targComb.rpartition('/')[-1]) ### Similar to parsed target finder format
    fh_out  = open(TarPred,'w')
    fh_out.write('miRname,Target,BindSite,miRseq,tarSeq,Score,Mismatch,CIGAR\n')
    
    #acount = 0
    #for i in fh_in:
    #    print ('Check lines:%s' % (i))
    #    acount +=1
    #print ('Total lines read:',acount)
    #sys.exit()
    
    #### Regenerate Target sequence with all features #####
    acount      = 0 ##Total number of interactions from predictions
    parseCount  = 0 ## Total number of interactions scores and written to result file
    for i in fh_in:
        # print("\nEntry:",i.strip("\n"))
        acount      += 1
        ent         = i.strip('\n').split('\t')
        #print('\n%s\n' % ent)
        miRrevcomp  = ent[9]                 ## miRNA complemented and reversed to map genome using bowtie. That is target sequence if mismatches and gaps are added
        tarHash     = list(miRrevcomp)          ## Strings are immutable convert to list - To rebuilt a traget seq
        # print("\nTarHash",tarHash)
        
        miRrev      = miRrevcomp.translate(str.maketrans("TACG","AUGC")) ## Re-translated to get miR but still in reverse orientation - OK
        mirHash     = list(miRrev)
        
        #print('Original read mapped i.e miRNA revcomp',miRrevcomp)
        
        ## Gap and Bulges (with reference to miRNA) - Identify gaps and bulges and modify miRNA read used for mapping to regenerate target as well as miRNA
        ## Add gap to target sequence  first to make miR length comparable to target
        gapinfo     = ent[5]
        gappos      = re.split("[A-Z]",gapinfo)   ## In python format - gap in target seq and bulge in miRNAseq
        gapNuc      = re.findall("[A-Z]",gapinfo)
        # print("gappos:",gappos,"| gapNuc:",gapNuc)
        
        ###########################################################################################
        ## SECTION - A - FIND GAPS AND BULGES AND ADD INDICATORS TO MIRNA OR TARGET SEQUENCES
        ###########################################################################################

        posCount = 0
        ## At this point both mirHash and tarHash has same length and perfect complementrity as tarHash is essentially reverse complemented miRNA used for matching
        for x,y in zip(gappos[:-1],gapNuc): ## In gap pos list which is made of alphabet splitting there is always am empty value at end because string has alphabet at last
            # print(x,y)
            if y == 'I':                    ## There was an insertion in miRNA and gap in reference and bulge in miRNA
                for i in range(int(x)):     ## For as many as bulges in miRNA - like 11M 2I 11M              
                    tarHash[posCount] = '-' ## Replace existing nucleotide (from revcomp miRNA) with a gap
                posCount += int(x)          ## This only effects consequitve bulges in miR, if there are multiple insertions like 3I, then "3" needs to be added once and not in every iterneration
            
            elif y == 'D':                  ## There was a deletion in miRNA i.e. gap in miRNA and bulge in reference - In this case length of both miRNA and target will increase
                for i in range(int(x)):             ## For as many as gaps in miRNA                
                    mirHash.insert(posCount,'-')    ## When counted in python insertion will occur after posCount value -  Tested OK
                    tarHash.insert(posCount,'^')    ## Add bulge markers in target sequence as well - Tested OK
                posCount += int(x)                  ## This only effects consequitve gaps in miR, if there are multiple insertions like 3I, then "3" needs to be added once and not in every iterneration
            
            else:
                posCount += int(x)
        
        # print('Target %s seq after manipulation: %s' % (ent[2],''.join(tarHash))) ## Has '-' at gap and '^' at extra nucleotide position (i.e. gap in MiR) -OK
        # print('miRNA %s after after maipulations: %s' % (ent[0],''.join(mirHash))) ## Has '-' at gap pos -OK
        
        #########################################################################################
        ## SECTION -B - GET CORRECT POSITIONS FOR MISMATCHES,GAPS AND BULGES 
        ## AND REGENERATE TARGET BY INSERTING CORRECT NUCLETIDES AT EDITS AND BULGES
        #########################################################################################
        
        ## Mismatches - Identify mismatches and modify miRNA read used for mapping to regenerate target
        misinfoBlock = ent[-2].split(':')[-1] ## Reverse index because XS:i is optional column ## MD:Z:16C3 - these positions are from references - so if there is an insertion/bulge in miRNA i.e. gap that it should be added to these positions before editing miRNA to tar
        # print ('This is the mismatch info block:%s' % (misinfoBlock))
        
        ## Deletion (gap) in miRNAS i.e y= D - which has been added but replace the bulge '^' in target with actual sequence
        ## Should work if there is dletion in miRNA and deletion in traget i.e. two bulges one in miRNA and one in target - NO MM possible
        if '^' in misinfoBlock: 
            misinfo = misinfoBlock.replace('^','')          ## 11^A13 - Here miRNA was 24nt but misinfo shows 25nt as 11+A+13 - Replace the '^' inserted in target with 'A'
            mispos = re.split("[A,T,G,C,N]",misinfo)        ## Found N in one case so included, N confimed in sequence too, will be counted as mismatch
            # print('Mismatch info:%s | Mismatch pos:%s'%(misinfo,mispos))       

            ## Add one to every position to get position where mismatch occured instead of position after which mismatch occured - This is an index and not position
            misposCorrect = []                              ## List hold corrected positions, because edit is a nucleotide next to integers in misPos
            # posIndex = -1                                   ## Index of position for misposCorrect, -1 because after first addition to list it will be incremented to 0
            for x in mispos:
                # print(x)
                ## Assumption: There will be no empty entry in mispos at the begining, because there will be a position to indicate comsequtive edits, others are handled here
                if x:
                    misposCorrect.append(int(x)+1)          ## Add one to every position to get position where mismatch occured instead of position after which mismatch occured - This is an index and not position
                    # posIndex +=1
                else: ## If 'x' is empty like in case of three consequentive gaps in miRNA - Mismatch info:11GTA8 | Mismatch pos:['11', '', '', '8']
                    # y = misposCorrect[posIndex]             ## Get the last corrected position, add one to get position for empty entry
                    # misposCorrect.append(int(y)+1)          ## In consequitve edits this is position just next to last one
                    # posIndex +=1
                    misposCorrect.append(int(0)+1)

            misNuc = re.findall("[A,T,G,C,N]",misinfo)      ## Found N in one case so included, N confirmed in sequence too, will be counted as mismatch
            # print('Misafter:',mispos,'Mispos', misposCorrect,'MisNuc',misNuc)

            ## Replace the nucleotides at bulges(^) and mismatches in tarHash to give actual targets
            ## And also add nucleotide to target (replace ^) if gap in miRNA
            #print('Unedited target:%s-%s' % (''.join(tarHash),len(''.join(tarHash))))
            posCount = 0
            for x,y in zip(misposCorrect,misNuc):
                posCount += x                               ## Convert bowtie positions to python format
                # print(x,y,posCount)
                ## Account for gap before replacing the nucleotide with that in target
                gaps = tarHash[:posCount-1-1].count('-')      ## -1 to convert to python, -1 because - count at positions in target before the current mismatch/bulge position
                tarHash[posCount-1+gaps] = y                  ## Replaced the bulge or mimatch with nucleotide in target - OK
        
        else:   ## Normal i.e y = I - With insertion(bulge) in miRNA and gap in target - OK - What id there is a bulge in miRNA???
                ## In this case miRNA already had inserted nucleltides and target has been added '-' in section-A. Just replace mimatches at correct postions of target
            
            misinfo = misinfoBlock
            mispos = re.split("[A,T,G,C,N]",misinfo)        ## Found N in one case so included, N confimed in sequence too, will be counted as mismatch   
            misposCorrect = [int(x)+1 for x in mispos]      ## Add one to every position to get position where mismatch occured instead of position after which mismatch occured - This is an index and not position
            misNuc = re.findall("[A,T,G,C,N]",misinfo)      ## Found N in one case so included, N confirmed in sequence too, will be counted as mismatch
            # print('Misafter:',mispos,'Mispos', misposCorrect,'MisNuc',misNuc)
            
            ## Count for gaps, since they are added by us and replace MM nucleotides to give actual targets
            posCount = 0
            for x,y in zip(misposCorrect,misNuc):
                posCount += x                                 ## Keep adding the positions, as these are cumulative - MD:Z:2T12C4 - Misafter: ['2', '12', '4'] Mispos [3, 13, 5] MisNuc ['T', 'C'] - Tested OK
                # print(x,y,posCount)
                ## Account for gaps before replacing the nucleotide with that in target
                ## Can give problem if more than one gap? But more than one gap not allowed V07 modification?
                gaps = tarHash[:posCount-1-1].count('-')      ## -1 to convert to python,-1 because - count at positions in target before the current mismatch/bulge position
                tarHash[posCount-1+gaps] = y                  ## TESTED - OK

        tar = ''.join(tarHash).replace("T","U")             ## Target converted to RNA format, will help in catching wobbles ahead
        mir = ''.join(mirHash)
        bindsite = '%s-%d' % (ent[3],int(ent[3])+(len(miRrev)-1))
        
        # print ("Target:%s-%s | miRNA:%s-%s" % (tar,len(tar),mir,len(mir)))
        
        ###################################################################################
        ## SECTION-C - GET POSITIONAL INFORMATION ON GAPS, MISMATCHES, WOBBLES AND BULGES 
        ## AND CHOOSE VALID INTERACTIONS 
        ###################################################################################

        mirGap = [] ## List to hold gaps in miRNA
        tarGap = [] ## List  to hold gap in targets
        mis = []    ## List to hold mismatch position
        wobble = [] ## List to hold Wobble pos
        # print('miRNA: %s\n%s' % (miRrevcomp[::-1],miRrevcomp[::-1].replace("T","U") ))

        ## Read from mapping file -> miRrevComp -> uncomplement -> miRrev -> morhash-> miR 
        valid = 1 ## Validity flag [0] - Invalid and [1] - valid, if more then 1 bulges in miR or 2bulges in tar or bulges in 9,10,11 = invalid
        nt_cnt = 1 ## Keep track of actual position,
        for x,y in zip(mir.translate(str.maketrans("AUGC","UACG",))[::-1],tar[::-1]): ## Orientation changed to read from 5' miRNA - OK
            
            #print(miRrev[::-1][nt_cnt-1],x,y)## Print miRNA, rev complemmnetry miRNA used for matching, target
            # if x == '-' or y == '-' or x == '^' or y == '^':
            #     #print('gap')
            #     gap.append(nt_cnt)
            #     if y == '-':
            #         nt_cnt+=1

            if x == '-' or x == '^': ## Don't think '^' would be here, since it has been replaced with nucleotide
                #print('miRNA gap')
                mirGap.append(nt_cnt)
                if nt_cnt == 9 or nt_cnt == 10 or nt_cnt == 11:
                    # print("@miRNA has gap/bulge in 9th and 10th position")
                    valid = 0 ## This is invalid as it has bulge in miRNA at 9th and 10th pos, this interaction will be skipped

            elif y == '-' or y == '^':
                #print('target gap')
                tarGap.append(nt_cnt)
                if y == '-': ## Don't think '^' would be here, since it has been replaced with nucleotide
                    nt_cnt+=1
                    if nt_cnt == 9 or nt_cnt == 10 or nt_cnt == 11:
                        # print("@target has gap/bulge in 9th and 10th position")
                        valid = 0 ## This is invalid as it has bulge in miRNA at 9th and 10th pos, this interaction will be skipped

            elif x == 'A' and y == 'G': ### If in reference its 'G' than miRNA should have 'U' i.e. T but this is revcomplememnt of miRNA so complement of 'U' is A - Tested OK - v08 modifcation
                #print ('wobble')
                wobble.append(nt_cnt)
                nt_cnt+=1
            elif x == 'C' and y == 'U': ### If in reference its 'U' than miRNA should have 'G' but this is rev complememnt of miRNA so complement of 'G' is C - Tested OK - v08 modification
                #print ('wobble')
                wobble.append(nt_cnt)
                nt_cnt+=1
            elif x == y:
                #print('match')
                nt_cnt+=1
            else:
                #print('mismatch')
                mis.append(nt_cnt)
                nt_cnt+=1
        # print('MismatchList:%s | miRGapList = %s | tarGapList = %s | WobbleList = %s' % (mis, mirGap, tarGap, wobble)) ## Poistion of mismatch gap and wobble

        ## Check if there are more then mpermitted gaps in miRNA and target
        if len(mirGap) > 2:
            # print("@miRNA %s has more then two gaps" % (ent[0]))
            valid = 0
        elif len(tarGap) > 1:
            # print("@target has more then one gap")
            valid = 0
        elif len(mis) > 5:
            # print("More then 5 mismatches not allowed")
            valid = 0
        elif len(mis) + len(wobble) > 6: 
            # print("Six edits are not allowed - It's too much")
            valid = 0
        else:
            pass

        ## Decide to report this interaction if its valid
        ## Validity flag, if more then 1 bulges in miR or 2bulges in tar or bulges in 9,10,11 = invalid
        if valid == 0: 
            # print("Skipping this entry - It's biologically invalid\n")
            continue ## Go to main for loop
        else:
            ## Go for scoring
            pass

        gap = mirGap+tarGap ## Combine gaps or scoring
        
        #####################################################################################
        ## SECTION-D - SCORE THE INTERACTIONS 
        #####################################################################################

        score = 0   ## Initialize
        #print (mis)
        if args.tarScore == 'S': ## Allowed 3 MM, 2 Wob, 1 Gap
            mis2 = list(mis)
            #if set([10,11]).issubset(mis): ## Works well but took 1 sec more than below in Rice timed test
            if 10 in mis and 11 in mis: ## Check for sunsequent mismatch at 10 and 11 if yes than strict penalty ## if set(['a','b']).issubset( ['b','a','foo','bar'] )
                score += 2.5
                #print('Removing 10')
                mis2.remove(10)
                #print ('Removing 11')
                mis2.remove(11) ## So that they are not counted again
                
            for i in mis2:
                    score += 1
            for i in gap:
                score += 1.5
            for i in wobble:
                if (i+1 in mis) or (i-1 in mis): ## Mismatches around wobble - Strong penalty
                    score += 1.5
                elif (i+1) in mis and (i-1 in mis): ## Mismatches on both sides - Stronger penalty
                    score += 2
                else:
                    score += 0.5
        else:
            ##Heuristic and Exhaustive
            for i in mis:
                if i>= 2 and i<=13:
                    score += 2
                else:
                    score += 1
            for i in gap:
                if i>= 2 and i<=13:
                    score += 2
                else:
                    score += 1
            for i in wobble:
                if i>= 2 and i<=13:
                    score += 1
                    #print ('Wobble pos:%s' % (i))
                else:
                    score += 0.5
                    #print ('Wobble pos:%s' % (i))
        
        ## Correctly output mismatches - If there is no mismatch then misinfo is just the length of match like - 15G2A1T1 (if there is mismatch) and 21(if no mismatch)
        ## Since the second ouput which originally is part of mapping file confussing in mismatch column - take mismatch info from 'mis' list
        
        if mis or gap: ## If list of mismathes has positions of mismatches, then output the block from mapping file
            mismatches = misinfo
        else:
            mismatches = 0

        #print(ent[0],ent[2],bindsite,miRrev,tar,score,misinfo,gapinfo)## MiRname, Tarname, mirSeq,Taerseq,binding site
        fh_out.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0],ent[2],bindsite,mir.replace("U","T"),tar.replace("U","T"),score,mismatches,gapinfo))
        parseCount  += 1
    
    
    print("Total number of interactions from 'miRferno':%s AND total interactions scored: %s" % (acount,parseCount))
    fh_in.close()
    fh_out.close()

    return TarPred

def tag2FASTA2(lib):
    print("'%s' tag count file being converted to FASTA format" % (lib))
    fh_in = open(lib,'r') #
    fh_out = open('./PARE/%s_PARE_tags.fa' % (lib), 'w')#
    tag_num = 1 ##
    for tag in fh_in:#
        #print(tag.strip('\n').split('\t'))
        ent = tag.strip('\n').split('\t')
        tag = ent[0]
        if len(tag) >= args.minTagLen: ##
            fh_out.write('>%s\n%s\n' % (tag_num, tag[:args.maxTagLen]))
            tag_num += 1
        else:
            #print ('Length is not 20nt')
            pass
    fh_out.close()

def mapdd2trans(anIndex):# 
    mismatch = str(0) #
    nspread2 = str(nspread)
    index = anIndex.rsplit('.', 2)[0]
    indexLoc = './index/%s' % index
    #for lib in libs:
    dd_file = ('./PARE/%s_PARE_tags.fa' % (templib))
    map_out = ('./dd_map/%s_%s_map' % (templib,index))
    print ('\n**The library %s is being mapped to transcriptome index file: %s**\n' % (dd_file,indexLoc))
    
    retcode2 = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nspread2, "-x", indexLoc, "-f", dd_file,"-S",map_out]) #
    #retcode2 = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "-t","-p",nspread2, "-f", indexLoc, dd_file,"-S",map_out]) #
    
    if retcode2 == 0:##
        print('\nDegradome from PARE lib: %s mapped to cDNA/Trascript file' % (templib)) 
    else:
        print ("There is some problem with mapping of PARE lib: %s to cDNA/genomic seq index" % (templib))
        print ("Script exiting.......")
        sys.exit()

def FileCombine():

    print('\n****************************************')
    targ_fls = [file for file in os.listdir('./predicted') if file.endswith ('.targ')]
    print ('Target files:',targ_fls)
    print ('\nCombining all the target prediction files for parsing and scoring\n')
    
    targComb = './predicted/All.targs'
    targ_out = open(targComb ,'w')
    
    for x in targ_fls:
        print (x)
        targfile = open('./predicted/%s' % (x), 'r')
        #targfile.readline()
        data = targfile.read()
        targfile.close()
        targ_out.write(data)
    
    targ_out.close()
        
    return targComb

def cleanHeader(filename,outfile):

    fh_in   = open(filename, 'r')
    fh_out  = open(outfile, 'w')
    
    print ('\nProcessing "%s" file to clean FASTA headers\n' % (filename))
    
    acount = 0 #
    for i in fh_in:
        if re.match('>', i):
            header = i.split()#
            new_head = header[0].split('|')[0]#
            fh_out.write('%s\n' % new_head)
            acount+=1

        else:
            fh_out.write('%s' % i)
        
    fh_in.close()
    fh_out.close()
    return outfile

    print('The fasta file with reduced header: "%s" with total entries %s has been prepared\n' % (out_file, acount))

def PP(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()
    nprocPP = round((args.accel/int(nspread))+1) #
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)
    
def PPmultiple(module,alist1,alist2):
    start = time.time()
    npool = Pool(int(args.accel))
    npool.map(lambda args: module(*args), alist2)

def PPResults(module,alist):
    npool = Pool(int(args.accel))    
    res = npool.map_async(module, alist)
    results = (res.get())
    npool.close()
    return results
    
def feed(queue, parlist):
    print ('Feeder function started')
    for par in parlist:
        print ('Echo from Feeder: %s' % (par))
        queue.put(par)
    print ('**Feeder finished queing**')

def calc(queueIn, queueOut):
    print ('Worker function started')
    while True:
        try:
            par = queueIn.get(block = False)
            print ('Echo from Worker \n Dealing with:', par)
            res = function(par)
            queueOut.put((par,res))
        except:
            break
    print ('**Worker finished **')

def write(queue, fname):
    print ('Writer function started')
    fhandle = open(fname, "w")
    while True:
        
        try:
            par, res = queue.get(block = False)
            print >>fhandle, par, res
        except:
            break
    fhandle.close()

def readFile(filename):
    """Read a file into memory

    Args:
        filename: Name of the file to be read
    Returns:
        Variable wholeFile containing the entire file

    """
    # home            = expanduser("~")
    # cwd             = os.getcwd()
    # print ("$HOME:%s | CWD:%s" % (home,cwd))
    f = open(filename, 'r')
    wholeFile = f.readlines()
    f.close()


    for i in range(len(wholeFile)) :
        wholeFile[i] = wholeFile[i].replace('\r', '')
        wholeFile[i] = wholeFile[i].replace('\n', '')
    
    return(wholeFile)

def createTargetFinderDataStructure(targetFinderFile):
    """Create data structure for targetFinder input

    Args:
        targetFinderFile: wholeFile as output from readFile function
    Returns:
        2D Array of the entire targetFinder input file

    """

    targetFinderList = []

    for i in range(len(targetFinderFile)):
        targetFinderList.append(targetFinderFile[i].split(','))

    #
    targetFinderListSorted = sorted(targetFinderList, 
        key=operator.itemgetter(0,5))

    return targetFinderListSorted

def pValueCalculator(target, targetFinderList, proportion):
    """ Function to calculate p-value

    Args:
        proportion: Proportion for given category

    Returns:
        p-value

    """

    # r = robjects.r ## Use scipy/numpy for this http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.binomial.html
    miRNAName   = target[0]
    score       = target[5]

    #
    if(float(score) % 1):
        n = sum(x.count(miRNAName) for x in targetFinderList if x[5].split('.')[0] == score.split('.')[0]) ## RH Change to allow X.0 and X.5 as equal.
    else:
        n = sum(x.count(miRNAName) for x in targetFinderList if float(x[5]) ==
            float(score))
    # pval = 1-(r.pbinom(0,n,proportion))[0]
    pval = 1-(scipy.stats.binom.pmf([0],n,proportion)[0])
    return pval

def validatedTargetsFinder(PAGeDict):
    """Perform the mapping. Take all entries from targetFinderList and
       identify if a target location matches to the 10th or 11th position
       of the miRNA.

    Args:
        PAGeDict: dictionary of PAGe file
        targetFinderList: list of target finder file
        categoryList: list of category proportions
    Returns:
        validated targets in form of targetFinder file with appended
        cleavage site and p-value.

    """ 

    cleaveLocations = {}
    cleaveLocations[21] = [9, 10, 11]
    cleaveLocations[22] = [9, 10, 11]
    cleaveLocations[23] = [9, 10, 11, 12]
    cleaveLocations[24] = [9, 10, 11, 12]
    cleaveStandard = [9, 10, 11]

    validatedTargets = []
    for target in targetFinderList:
        gene = target[1]
        # Get the start location of the target
        end = int(target[2].split('-')[1])
        # Get the length of the mIRNA
        length = int(len(target[3]))

        if(gene in PAGeDict.keys()):
            currDict = PAGeDict[gene]
            cleavageSite = []
            location = 0
            targetAbundances = []
            targetCategories = []
            targetLocations = []
            #
            #
            #

            # If the length of the miRNA is in our dictionary and the cleave
            # location exists in the dictionary, add it to the target lists
            if(not args.standardCleave and length in cleaveLocations.keys()):
                for cleaveLocation in cleaveLocations[length]:
                    if(str(end-cleaveLocation) in currDict):
                        targetAbundances.append(currDict[str(end -
                            cleaveLocation)][0])
                        targetCategories.append(currDict[str(end -
                            cleaveLocation)][1])
                        targetLocations.append(end - cleaveLocation)

            # If the length of the miRNA is not in our dictionary, we will 
            # just investigate the 10th 11th and 12th positions 
            else:
                for cleaveLocation in cleaveStandard:
                    if(str(end-cleaveLocation) in currDict):
                        targetAbundances.append(currDict[str(end -
                            cleaveLocation)][0])
                        targetCategories.append(currDict[str(end -
                            cleaveLocation)][1])
                        targetLocations.append(end - cleaveLocation)
                

            # If there is a PARE cleavage at any of the above positions,
            # find the best candidate for retainer.
            if(targetCategories):
                ## Debugging statement retained for quick analysis. Use with
                ## Debugging code below
                #print(targetAbundances, targetCategories, targetLocations)

                # If there is only one minimum category, use this target as
                # most probably cleavage location
                if(targetCategories.count(min(targetCategories)) == 1):
                    cleavageIndex = targetCategories.index(min(
                        targetCategories))
                    location = targetLocations[cleavageIndex]
                    cleavageSite = currDict[str(location)]

                # If there is more than one minimum category, we must filter
                # further to base our target of interest on greatest read
                elif(targetCategories.count(min(targetCategories)) > 1):
                    # Get all indices of the minimum category
                    cleavageIndices = [i for i, x in enumerate(
                        targetCategories) if x == min(targetCategories)]

                    # Get list of abundances with minimum categories
                    abundances = [targetAbundances[index] for index in 
                        cleavageIndices]
                    
                    # Cleavage index will be the read with the greatest
                    # abundance. If there is a still a tie, use the 
                    # index of the category with greatest read abundance.
                    # If there is a tie, the lowest index will always be used.
                    cleavageIndex = targetAbundances.index(max(abundances))
                    location = targetLocations[cleavageIndex]
                    cleavageSite = currDict[str(location)]
            
            if(location):
                ## Debugging statement retained in conjunction with above.
                ## Shows if cleavage is 10th, 11th or 12 position. (not
                ## coordinated with locations that output prior.)
                #print(cleavageSite, -(location-end) + 1)
                windowSum = 0
                toAppend = list(target)
                # The category score is the 2nd position.
                categoryScore = cleavageSite[1]
                # Calulate the p-value.
                pValue = pValueCalculator(target, targetFinderList, 
                    categoryList[int(categoryScore)])

                for i in range(location-5, location+6):
                    if(str(i) in currDict):
                        windowSum += int(currDict[str(i)][0])

                # Add one to the location because we need to account for the
                # 0th position of the index     
                toAppend.append(str(location)) ### Editing for biological understanding - original (str(location+1))
                # Add PARE abundance
                toAppend.append(cleavageSite[0])
                # Add sum of reads within 5 bp of cleavage site in each
                # direction.
                toAppend.append(str(windowSum))
                # Add ratio of abundance of cleavage site to sum within 5 bp 
                # of the cleavage site in each direction
                toAppend.append(str("%f" % float(int(cleavageSite[0])/windowSum)))
                # Add category at cleavage site
                toAppend.append(str(categoryScore))
                # Append the p-value to the toAppend list
                toAppend.append(str("%f" % pValue))
                validatedTargets.append(toAppend)
    
    return(validatedTargets)

def createPAGeIndex(bowtieFilename):
    """Create data structure for targetFinder input. Then, performs the
       mapping. Take all entries from tagCountFile and identify all tags
       that map to 
    Args:
         bowtieFilename: name of map file
    Returns:
        Dictionary with the proper information for dd_PAGe output. Key
        is gene name and the value is a list with tuples as values. The first
        element in the tuple is the location of the hit, and the second
        element is the abundance for that sequence.
        A list of all abundance values are also returned
    """

    #
    bowtieDict = {}
    
    #
    bowtieFilename = 'dd_map/' + bowtieFilename
    bowtieFile = readFile(bowtieFilename)

    if args.repeats:
        #
        for entry in bowtieFile:
            if(entry.split('\t')[4] != '255'):
                pass
            else:
                #
                gene = entry.split('\t')[2]
                location = entry.split('\t')[3]
                sequence = entry.split('\t')[9]

                #
                try:
                    bowtieDict[sequence]
                #
                except:
                    bowtieDict[sequence] = {}

                # 
                try:
                    bowtieDict[sequence][gene]

                # 
                except:
                    bowtieDict[sequence][gene] = []

                #
                bowtieDict[sequence][gene].append(location)        

    else:
        for entry in bowtieFile:
            #
            gene = entry.split('\t')[2]
            location = entry.split('\t')[3]
            sequence = entry.split('\t')[9]

            #
            try:
                bowtieDict[sequence]
            #
            except:
                bowtieDict[sequence] = {}

            #
            try:
                bowtieDict[sequence][gene]

            #
            except:
                bowtieDict[sequence][gene] = []

            #
            bowtieDict[sequence][gene].append(location)

    PAGeDict = {}
    allHits = []

    #
    for entry in tagCountFile:
        #
        sequence = entry.split('\t')[0][:args.maxTagLen]
        hits = entry.split('\t')[1]

        # 
        # 
        try:
            #
            for key in bowtieDict[sequence].keys():
                # 
                # 
                try:
                    PAGeDict[key]
                except:
                    PAGeDict[key] = {}

                # 
                # 
                for location in bowtieDict[sequence][key]:
                    # 
                    PAGeDict[key][location] = int(hits)
                   
                    # Append the hits to the hits if hits > 2
                    if(int(hits) > 2):
                        allHits.append(int(hits))

        except:
            pass

    return(PAGeDict, allHits)

def unambiguousBaseCounter(transcriptomeFilename, minTagLen):
    """Get the counts of ambiguous bases in the transcriptome file as well
       as counts of ambiguous bases that are within the ends of the
       transcriptome - the minTagLen.

    Args:
        transcriptomeFilename: Name of the trnascritome fasta file
        minTagLen: Number of bases that an N is allowed to be away from the ends
            of the gene in order to be counted

    Returns:
        Total number of ambiguous bases and ambiguous bases minTagLen-bp away from
        the ends of the gene.

    """
    baseCounts, baseCountsOffTagLen = 0, 0

    transcriptomeFile = readFile(transcriptomeFilename)
    for i in range(1,len(transcriptomeFile),2):
        currentLine = transcriptomeFile[i]
        baseCounts += len(currentLine) - currentLine.count('N')
        baseCountsOffTagLen += (len(currentLine) - 2 * minTagLen) - currentLine[
            minTagLen:len(currentLine)-minTagLen].count('N')

    f_output = open('baseCounts.mem', 'w')
    f_output.write(str(baseCounts) + '\n' + str(baseCountsOffTagLen))
    f_output.close()

def writePAGeFile(PAGeDict, mode, allHits, baseCounts, baseCountsOffTagLen,
    outputFile, transcriptomeFilename, library):
    """Write validated targets to an output file

    Args:
        PAGeDict: Dictionary of genes with target locations and hits for
            each tag mapped to those locations.
        mode: 0 (genic) or 1 (intergenic)
        allHits: List of all abundance values
        baseCounts: Total number of unambiguous bases
        baseCountsOffTagLen: Total number of unambiguous bases minTagLen-bp away from
            the ends of the gene.
        outputFile: file to output PAGe inforation
        transcriptomeFilename: Name of transcriptome file
        library: Name of library being analyzed

    """
    # 
    categoryCounts = [0,0,0,0,0]
    categoryList = []
    f = open(outputFile,'w')
    #g = open('output/%s_categoryCounts.csv' % library,'w')
    # Get the count of the total genes
    numGenes = len(PAGeDict.keys())

    if(mode == 1):
        globalMedian = numpy.median(allHits)
        seventyFivePercentile = stats.scoreatpercentile(allHits, 75)
        ninetyPercentile = stats.scoreatpercentile(allHits, 90)
        print('median = %s\nseventyFivePercentile = %s\nninetyPercentile = %s'
            % (globalMedian, seventyFivePercentile, ninetyPercentile))
        #g.write('gene,cat0,cat1,cat2,cat3\n')

    # Sort genes so that genes are in alphabetical order
    for gene in sorted(PAGeDict.keys()):
        # Dictionary for each gene which records all locations
        # with the same abundance and categories for output
        hitsDict = {}
        catTwoAbun = []
        catThreeAbun = []
        # Genic tracks category 0 and 1 as one
        if(mode == 0):
            catZeroOneAbun = []
        # Intergenic separates categories 0 and 1
        else:
            catZeroAbun = []
            catOneAbun = []
       
        f.write('>%s\n' % str(gene))
        #g.write(str(gene) + ',')
        if(mode == 0):
            geneHits = []
            multHitsFlag = 0

            # Store all of the abundance values for each location
            for hits in PAGeDict[gene].values():
                geneHits.append(hits)

            # Calculate median and max on gene
            median = numpy.median(geneHits)
            maxHit = max(geneHits)
            if(len([i for i, x in enumerate(geneHits) if x == maxHit]) > 1):
                multHitsFlag = 1

            # Sort all locations in which a tag maps in that gene so they can
            # be listed in order. Must use key as int because in dictionary,
            # locations are stored as strings
            for location in sorted(PAGeDict[gene].keys(),key=int):
                hits = PAGeDict[gene][location]

                # Calculate category
                if(hits == 1):
                    category = '4'
                    categoryCounts[4] += 1
                    PAGeDict[gene][location] = (hits, 4)
                elif(hits <= median):
                    category = '3'
                    categoryCounts[3] += 1
                    PAGeDict[gene][location] = (hits, 3)
                    catThreeAbun.append(hits)
                elif(hits > median and hits != maxHit):
                    category = '2'
                    categoryCounts[2] += 1
                    PAGeDict[gene][location] = (hits, 2)
                    catTwoAbun.append(hits)
                elif(hits > median and multHitsFlag):
                    category = '1'
                    categoryCounts[1] += 1
                    PAGeDict[gene][location] = (hits, 1)
                    catZeroOneAbun.append(hits)
                else:
                    category = '0'
                    categoryCounts[0] += 1
                    catZeroOneAbun.append(hits)
                    PAGeDict[gene][location] = (hits, 0)
                #f.write('%s\t%s\t%s\n' % (str(location), str(hits), category))

                # If the abundance already exists, append the location to the
                # location list in the tuple. Otherwise, create the entry and
                # add the category
                try:
                    hitsDict[hits][0].append(location)
                except:
                    hitsDict[hits] = ([location], category)

            #g.write(str(max(catZeroOneAbun) if catZeroOneAbun else 0) + ',' +
                #str(max(catTwoAbun) if catTwoAbun else 0) + ',' +
                #str(max(catThreeAbun) if catThreeAbun else 0) + '\n')

        elif(mode == 1):
            # Sort all locations in which a tag maps in that gene so they can
            # be listed in order.
            for location in sorted(PAGeDict[gene].keys(), key=int):
                hits = PAGeDict[gene][location]

                # Calculate category
                if(hits <= 2):
                    category = '4'
                    categoryCounts[4] += 1
                    PAGeDict[gene][location] = (hits, 4)
                elif(hits <= globalMedian):
                    category = '3'
                    categoryCounts[3] += 1
                    PAGeDict[gene][location] = (hits, 3)
                    catThreeAbun.append(hits)
                elif(hits > globalMedian and hits <= seventyFivePercentile):
                    category = '2'
                    categoryCounts[2] += 1
                    PAGeDict[gene][location] = (hits, 2)
                    catTwoAbun.append(hits)
                elif(hits > seventyFivePercentile and
                        hits <= ninetyPercentile):
                    category = '1'
                    categoryCounts[1] += 1
                    PAGeDict[gene][location] = (hits, 1)
                    catOneAbun.append(hits)
                else:
                    category = '0'
                    categoryCounts[0] += 1
                    PAGeDict[gene][location] = (hits, 0)
                    catZeroAbun.append(hits)
                #f.write('%s\t%s\t%s\n' % (str(location), str(hits), category))

                # If the abundance already exists, append the location to the
                # location list in the tuple. Otherwise, create the entry and
                # add the category
                try:
                    hitsDict[hits][0].append(location)
                except:
                    hitsDict[hits] = ([location], category)

        for hits in sorted(hitsDict.keys(), key=int):
            f.write('%s\t' % hits)
            f.write(','.join(sorted(hitsDict[hits][0], key=int)))
            f.write('\t%s\n' % hitsDict[hits][1])
            #f.write('%s\t%s\t%s\n' % (hits, sorted(hitsDict[hits][0], key=int), hitsDict[hits][1]))
            #for i in range(len(hitsDict[hits][0])):
            #    loc = sorted(hitsDict[hits][0])[i]
            #    if(i == len(hitsDict[hits][0]) - 1):
            #        print(i)
            #        f.write('%s\t%s;\n' % (loc, hitsDict[hits][1]))
            #    else:
            #        f.write('%s,' % loc)

        hitsDict.clear()
            #g.write(str(max(catZeroAbun) if catZeroAbun else 0) + ',' + 
                #str(max(catOneAbun) if catOneAbun else 0) + ',' +
                #str(max(catTwoAbun) if catTwoAbun else 0) + ',' + 
                #str(max(catThreeAbun) if catThreeAbun else 0) + '\n')

    f.write('# Transcriptome=%s\n' % transcriptomeFilename)
    f.write('# Genes=%s\n' % numGenes)
    f.write('# Uncorrected non-ambiguous bases=%s\n' % baseCounts)
    f.write('# Eligible bases for degradome-derived 5 prime ends=%s\n' %
        baseCountsOffTagLen)
    for i in range(len(categoryCounts)):
        f.write('# Category %s_bases=%s\n' % (i, categoryCounts[i]))
    for i in range(len(categoryCounts)):
        categoryList.append(categoryCounts[i] / baseCountsOffTagLen)
        f.write('# Category %s_fraction=%s\n' % (i, categoryCounts[i] / 
            baseCountsOffTagLen))
    
    f.close()
    #g.close()
    return(categoryList)

def writeValidatedTargetsFile(header, validatedTargets, outputFile):
    """Write validated targets to an output file

    Args:
        header: the stripped header from the target file
        validatedTargets: list of the validated targets
        outputFile: file to output validated targets to

    """

    # 
    pValueIndex = len(validatedTargets[0]) - 1
    categoryIndex = len(validatedTargets[0]) - 2
    windowRatioIndex = len(validatedTargets[0]) - 3
    
    #
    for i in range(len(validatedTargets)):
        validatedTargets[i].append("%.4f" % (float(validatedTargets[i][pValueIndex])/float(validatedTargets[i][windowRatioIndex])))

    #
    correctedPValueIndex = len(validatedTargets[0]) - 1

    f = open(outputFile,'w')

    f.write(header + ',cleavagePosition,PAREAbundance,10-nt window abundance,'\
        'PARE reads/window abundance,category,p-value,corrected p-'\
        'value\n')

    # noiseFilter requires that p value be <= .25 and window ratio be
    # >= .25 for a tag to be considered
    if(args.noiseFilter):
        for target in validatedTargets:
            target[windowRatioIndex] = '%.3f' % float(target[windowRatioIndex])
            target[correctedPValueIndex] = '%.4f' % float(
                target[correctedPValueIndex])
            # Include any target with pvalue <= .25 and with window ratio
            # >= .25
            if((float(target[correctedPValueIndex]) <= .25) and (float(target[
                    windowRatioIndex]) >= .25)):
                # If category 4s are not to be included and the category
                # of the current category is 4, then omit it.
                if(args.cat4Show and (int(target[categoryIndex]) == 4)):
                    pass
                else:
                    for j in range(len(target)):
                        f.write(str(target[j]))
                        # If at last element, write new line
                        if(j == correctedPValueIndex):
                            f.write('\n')
                        # Otherwise, write a comma
                        else:
                            f.write(',')

    # If noise filter is not on, then include any target with p value < .5 
    else:
        for target in validatedTargets:
            # Include any target with pvalue < .5
            if(float(target[correctedPValueIndex]) < .5):
                for j in range(len(target)):
                    f.write(str(target[j]))
                    # If at last element, write new line
                    if(j == correctedPValueIndex):
                        f.write('\n')
                    # Otherwise, write a comma
                    else:
                        f.write(',')

def resultUniq(filetag):
    
    """Read validated files for each library 
    and generate a single file with unique results"""

    print("\n#### Fn: resultUniq ####################################")

    # fls = [file for file in os.listdir('./output') if re.search(r'revmapped.csv', file)]
    fls = [file for file in os.listdir('./output') if file.endswith (filetag)]
    
    print ('Combining results from all the files to generate a single report\n')
    print ('--Files with lib-wise results:',fls)

    
    validatedComb   = './output/temp.csv'
    fh_out          = open(validatedComb ,'w')
    
    ## Combine files
    header = "" 
    for x in fls: 
        # print (x)
        afile   = open('./output/%s' % (x), 'r')
        header          = afile.readline().strip('\n') ## Use later
        data            = afile.read()
        afile.close()
        fh_out.write(data)
    fh_out.close()
    
    ## Sort combined result file:
    fh_in       = open(validatedComb, 'r') 
    parsed_in   = [line.strip('\n').split(',') for line in fh_in]
    parsed_in.sort(key=lambda k: (-int(k[9]) )) ## PARE and corrected p-value
    parsed_in.sort(key=lambda k: (float(k[14]) ))

    uniqRevmapped   = './output/All.libs.validated.uniq.csv'
    fh_output2      = open(uniqRevmapped, 'w')
    fh_output2.write("%s\n" % header)

    ## Uniq
    added_keys=set()## A set to store first 3 elements from input file: miRNA-DNA, chr# and cleavage site and than use it to compare further entries in file
    parsed_out_count=0## To keep count of unique entries
    for ent in parsed_in:
        #print(ent[0],ent[1],ent[15],ent[18])
        genename = ent[1] ## To avoid different variations of same gene to be counted as uniq
        lookup=tuple((ent[0],genename,ent[8]))## miR name + Target Gene + position of cleavage on gene
        if lookup not in added_keys:
            fh_output2.write('%s\n' % (','.join(ent)))
            parsed_out_count+=1
            added_keys.add(lookup)## Once a new entry is found it is recorded so as to compare and neglect further entries
        else:
            pass

    if os.path.isfile(validatedComb):
        os.remove(validatedComb)

    fh_in.close()
    fh_out.close()
    fh_output2.close()
    
    return uniqRevmapped

def genomicCoord(ent): ####

    '''
    Reverse maps entry
    '''
    
    ## Create a dictionary from list of coords to be searched later
    ## Gene_coords structure: 1, 'c','AT1G01020', 5928, 8737, protein_coding
    ## ent structure: >ath-miR401,AT2G06095,971-991,ACAGCCAGCTGTGGTCAAAGC,TGTCGATCGACACCAGTTTCG,1,5A15,21M,981,13,13,1,2,0.000217,0.000217
   
    
    # print (ent)
    # gene_name   = ent[1].split('_')[0] # for exmaple, AT2G06095_up - Changed to nect line in sPARTA1.17
    gene_name   = ent[1] # for exmaple, AT2G06095_up
    bind_site   = ent[2].split('-')
    cleave_site = int(ent[8])
    
    ## Reverse map co-ordinates ##########################################################
    # print ('**Reverse mapping of Co-ordinates will be performed**')
    if gene_name in coord_dict_wat:
        # print ('Entry: %s in positive strand: %s' % (ent[0:4],coord_dict_wat[gene_name]))
        geno_start = coord_dict_wat[gene_name][1]###Use for reverse mapping of postive genes

        #print('Looking for chr_id')
        chr_id  = coord_dict_wat[gene_name][0]
        #print('chr_id found')
        strand  = 'w' ## AVlilable in dictionary also coord_dict_crick[gene_name][1]
        gtype   =  coord_dict_wat[gene_name][2] ## Gene type
        new_cleave_site     = (int(geno_start)-1)+int(cleave_site)###1 is reduced to give correct positions
        new_bind_site_start = (int(geno_start)-1)+int(bind_site[0])
        new_bind_site_end   = (int(geno_start)-1)+int(bind_site[1])
        new_bind_site       = '%s-%s' % (new_bind_site_start,new_bind_site_end)
    else:
        # print ('Entry: %s in reverse strand: %s' % (ent[0:4],coord_dict_crick[gene_name]))
        geno_end    = coord_dict_crick[gene_name][2] ### Use for reverse mapping of negative genes
        #print('Looking for chr_id')
        chr_id      = coord_dict_crick[gene_name][0]
        #print('chr_id found')
        strand      = 'c' ## Available in dictionary also coord_dict_crick[gene_name][1]
        gtype       =  coord_dict_crick[gene_name][2] ##Gene type
        new_cleave_site     = (int(geno_end)+1)-int(cleave_site)###1 is added to give correct positions
        new_bind_site_end   = (int(geno_end)+1)-int(bind_site[0])###As the sequence was reversed before TF and CL, their binding start and end direction has also changed - Verified-OK
        new_bind_site_start = (int(geno_end)+1)-int(bind_site[1])
        new_bind_site       = '%s-%s' % (new_bind_site_start,new_bind_site_end)

    # print("Rev Mapped: %s,%s,%s,%s,%s" % (str(chr_id),strand,new_cleave_site,new_bind_site_start,new_bind_site_end))
    rev_mapped_entry = ("%s,%s,%s,%s,%s,%s" % (','.join(ent),str(chr_id),strand,new_cleave_site,new_bind_site_start,new_bind_site_end))
    
    return rev_mapped_entry

def ReverseMapping():
    '''
    Creates Coords dictionary, Reversemap coordinates, and output results in file
    '''

    print("#### Reverse Mapping Initiated ###########################################")
    
    ## Read strand from gene coords file becuase uniqRevmapped file doesn't have this info.
    print("Step 1/4: Reading coords file")
    fh_in       = open("coords","r")
    fileRead    = fh_in.readlines()
    coords      = [] ## List to store strand

    for i in fileRead:
        ent     = i.strip("\n").split(",")
        # print("This is read entry",ent)
        coords.append((ent))
    print("Coords read from file:%s entries" % (len(coords)))
    print('snippet of coords:',coords[1:10])
    print("Step 1/4: DONE!!\n\n")
    time.sleep(1)

    #### 2. PREPARE DICTIONARY OF COORDS ########################################
    print("Step 2/4: Preparing dictionary of coordinates")
    global coord_dict_wat, coord_dict_crick
    coord_dict_wat      = {} ## Dictionary of genes at watson strand
    coord_dict_crick    = {} ## Dictionary of genes at crick strand

    global nproc
    nproc ='1' ## Need better handling
       
    for i in coords:### gene_coords is a list in script, also written out as file of same name
        # print ("This is a coord:",i)
        strand      = i[1]
        # GeneName    = i[2].split('_')[0] # for exmaple, AT2G06095_up
        GeneName    = i[2] # for exmaple, AT2G06095_up
        # print('chr:%s | gene:%s | strand:%s' % (i[0],i[2],i[1]))
        if strand == 'c':### if entry in reverse strand
            atuple = (i[0],i[3],i[4],i[5])
            coord_dict_crick[GeneName] = atuple ## Gene name as key and chr_id,start, end and gene type as value
        elif strand == 'w':
            atuple = (i[0],i[3],i[4],i[5])
            coord_dict_wat[GeneName] = atuple   ## Gene name as key and chr_id,start, end and gene 
        else:
            print("Something went wrong with reverse mapping")
            print("Encountered unexpected charater for chromosome strand:%s" % (strand))
            print("Please check your annotation file, it seems to have unexpected characters for strand")
            sys.exit()
    
    print("Entries in watson dict:%s | crick dict:%s" % (len(coord_dict_wat),len(coord_dict_crick)))
    print("Strand dictionary made")
    print("Step 2/4: DONE!!\n\n")
    time.sleep(1)
    # sys.exit()

    ##### 3. Read the scoring input extend file and change coordinates ##########
    print("Step 3/4: Reading validated targets in list")
    ResFls = [file for file in os.listdir('./output') if file.endswith ('_validated')]
    for afile in ResFls:
        print("Caching results '%s' for reverse mapping"% (afile))
        fh_in       = open('./output/%s' % afile, 'r')### PARE VALIDATED results
        header      = fh_in.readline() ## Waste header
        ScoInpExt   = [] ## list of final results or parlist
        acount      = 0  ## All entries
        for res in fh_in:
            acount      +=1
            res_strp    = res.strip('\n')
            ent         = res_strp.split(',')
            ScoInpExt.append(ent)


        print("%s results in cached from validated file:%s" % (afile,acount))
        print("Step 3/4: DONE!!\n\n")
        time.sleep(1)

        ## Rev Map #######
        ##################
        print("Step 4/4: Reverse mapping using coords dict and targets list")
        
        ## # TEST- SINGLE CORE - For troubleshooting ----##
        print('\n***Entering genomicCoord- Serial***\n')
        ValidTarGeno = []
        for i in ScoInpExt:
            print("PARE lib:",afile)
            # print("\nEntry for reverse mapping:%s" % (i))
            z = genomicCoord(i)
            ValidTarGeno.append(z)
        
        ## NORMAL- PARALLEL MODE - Comment test above for normal use  
        # print('\n***Entering genomicCoord - parallel***\n')
        # print('**Reverse mapping initiated**\n\n')
        # ValidTarGeno = PPResults(genomicCoord, ScoInpExt) ## Results are in form of list
        
        print ('Reverse mapping complete for:%s\n\n\n' % (afile))
        # print("Step 4/4: Done!!\n\n")
        time.sleep(1)

        #### WRITE RESULTS ##############
        ################################

        print ("Writing Results")
        # os.remove(uniqRevmapped) #Remove previous "AllLibValidatedUniq_revmapped.csv"    
        revmapRes   = './output/%s_revmapped.csv' % (afile)
        fh_out      = open(revmapRes, 'w')
        fh_out.write('%s,Chr,Strand,genomeCleavePosition,genomeBindStart,GenomeBindEnd\n' % (header.strip('\n')))
        
        for i in ValidTarGeno: ## Write Results from list to file
            if i != 'E13-3-13':##Filter for error 13-13-13 that is small window abundance = 0 and ratio calculation error
                fh_out.write('%s\n' % (i))
            else:
                print (i)
        fh_in.close()
        fh_out.close()

        ## Remove non-revmapped file (_validated files)
        if os.path.isfile('./output/%s' % afile):
            # print("validated file removed:%s" % (afile))
            os.remove('./output/%s' % afile)

    return None

    
    ### END OF REVERSE MAP ##########

def dedup_process(alib):
    '''
    To parallelize the process
    '''
    print("\n#### Fn: De-duplicater #######################")

    afastaL     = dedup_fastatolist(alib)         ## Read
    acounter    = deduplicate(afastaL )            ## De-duplicate
    countFile   = dedup_writer(acounter,alib)   ## Write

    return countFile

def dedup_fastatolist(alib):
    '''
    New FASTA reader
    '''

    ### Sanity check
    try:
        f = open(alib,'r')
    except IOError:                    
        print ("The file, %s, does not exist" % (alib))
        return None


    ## Output 
    fastaL      = [] ## List that holds FASTA tags

    print("Reading FASTA file:%s" % (alib))
    read_start  = time.time()
    
    acount      = 0
    empty_count = 0
    for line in f:
        if line.startswith('>'):
            seq = ''
            pass
        else:
          seq = line.rstrip('\n')
          fastaL.append(seq)
          acount += 1

    read_end    = time.time()
    # print("-- Read time: %ss" % (str(round(read_end-read_start,2))))
    print("Cached file: %s | Tags: %s | Empty headers: %ss" % (alib,acount,empty_count)) 

    return fastaL
                   
def deduplicate(afastaL):
    '''
    De-duplicates tags using multiple threads and libraries using multiple cores
    '''
    dedup_start  = time.time()

    # deList = [] ## Hold deduplicated tags and their abudnaces in a tuple

    acounter    = collections.Counter(afastaL)

    dedup_end  = time.time()
    # print("-- dedup time: %ss" % (str(round(dedup_end-dedup_start,2))))

    return acounter 

def dedup_writer(acounter,alib):
    '''
    writes rtag count to a file
    '''

    print("Writing counts file for %s" % (alib))
    countFile   = "%s.fas" % alib.rpartition('.')[0]  ### Writing in de-duplicated FASTA format as required for phaster-core
    fh_out       = open(countFile,'w')

    acount      = 0
    seqcount    = 1 ## TO name seqeunces
    for i,j in acounter.items():
        # fh_out.write("%s\t%s\n" % (i,j))
        fh_out.write(">seq_%s|%s\n%s\n" % (seqcount,j,i))
        acount      += 1
        seqcount    += 1

    print("Total unique entries written for %s: %s" % (alib,acount))

    fh_out.close()

    return countFile

def coreReserve(cores):
    '''
    Decides the core pool for machine - written to make PHASworks comaptible with machines that 
    have less than 10 cores - Will be improved in future - Used as is from phasdetect (PHASworks)
    '''

    # if cores == 0:
    if cores == "Y":
        ## Automatic assignment of cores selected
        totalcores = int(multiprocessing.cpu_count())
        if totalcores   == 4: ## For quad core system
            nproc = 3
        elif totalcores == 6: ## For hexa core system
            nproc = 5
        elif totalcores > 6 and totalcores <= 10: ## For octa core system and those with less than 10 cores
            nproc = 7
        else:
            nproc = int(totalcores*0.85)
    else:
        ## Reserve user specifed cores
        nproc = int(cores)

    return nproc

##############################################################################################
#### MAIN FUNCTION ###########################################################################
def main():

    ## Pre-run check and imports ######################
    checkUser()
    checkLibs()
    checkfiles()
    global numpy,stats,scipy
    from scipy import stats
    import numpy,scipy
    # import rpy2.robjects as robjects
    # from rpy2.robjects.packages import importr
    # from rpy2.robjects.vectors import FloatVector
    # RStats = importr('stats')
    ###################################################


    if args.generateFasta:
        chromoDict,chromoLenDict          = genomeReader(args.genomeFile)
        # If the annotatyion type is a GFF file, run the GFF parser
        if(args.annoType == 'GFF'):
            genome_info,genome_info_inter = gffParser(args.annoFile)
        # If the annotatyion type is a GTF file, run the GTF parser
        elif(args.annoType == 'GTF'):
            genome_info,genome_info_inter = gtfParser(args.annoFile) 
        coords                            = extractFeatures(args.genomeFile,chromoDict,chromoLenDict,genome_info,genome_info_inter) ## Extracts Coords from GFF3
        fastaOut,fastaList              = getFASTA1(args.genomeFile,coords,chromoDict) ##Creates FASTA file
        unambiguousBaseCounter(fastaOut, args.minTagLen)
        print('This is the extracted file: %s' % (fastaOut))
    # 
    elif args.featureFile:
        print("\nThe input FASTA file is considered 'as is' for analysis\n")
        fastaOut    = args.featureFile ### Make it better
        fastaList   = fastaReader(fastaOut)
        unambiguousBaseCounter(fastaOut, args.minTagLen)
    else:
        print("Please provide input to '--featureFile' or '--genomeFile'")
        print("See sPARTA example commands: https://github.com/atulkakrana/sPARTA.github/tree/master/sparta")
        print("Script will exit now")
    
    ### FRAGMENTATION ###################
    ### Script Timer
    runLog      = 'runtime_%s' % datetime.datetime.now().strftime("%m_%d_%H_%M")
    fh_run      = open(runLog, 'w')
    fh_run.write('tarPred:%s | tarScore: %s | Uniq filter:%s\nLibs:%s\nGenomeFile:%s | GenomeFeature:%s' % (args.tarPred,args.tarScore,args.repeats, ','.join(args.libs),args.genomeFile,args.genomeFeature))
    fh_run.write ('\nLibs: %s' % (','.join(args.libs)))
    FragStart   = time.time()
    
    if args.fileFrag:
        start       = time.time()###time start
        fragList    = fragFASTA(fastaOut,fastaList)##
        end         = time.time()
        print ('fileFrag time: %s' % (round(end-start,2)))
    else:
        fragMem     = open('frag.mem','r')
        fragRead    = fragMem.readlines()
        print("Frags from earlier processed file: '%s' will be used in this run" % (fragRead[0].strip('\n').split('=')[1]))
        for i in fragRead:
            akey,aval = i.strip('\n').split('=')
            if akey == "genomeFeature":
                if int(aval) != args.genomeFeature:
                    print('\n------ERROR-------')
                    print("The earlier processed genome or feature set belongs to genomeFeature: %s" % (aval))
                    print("Your current genomeFeature input: %s" % (str(args.genomeFeature)))
                    print("Please input either correct genomeFeature value or re-run the analysis with '-featurefile' or '-genomeFile and -annoFile'")
                    print("\nSystem will exit now")
                    sys.exit(0)
                else:
                    pass
                
            if akey =="frags":
                fragList =  i.split('=')[1].split(',')
                #print('fragList from fragMem:', fragList)
        fragMem.close()
        #fragList = [file for file in os.listdir() if re.search(r'.*\.frag\.[0-9].*\.fa', file)] # deprecated v1.03
        print ('The fragments: %s' % (fragList))
        
    FragEnd = time.time()
    #print ('\n\nThe fragmentation time is %s\n\n' % (round(FragEnd-FragStart,2)))
    fh_run.write('Fragmentation time is : %s\n' % (round(FragEnd-FragStart,2)))

    #####################################
    
    ## TARGET PREDICTION ################### 
    TFStart = time.time()
    print('tarPred: %s | tarScore: %s | Uniq filter: %s' % (args.tarPred,args.tarScore,args.repeats))
    ## 
    if args.indexStep:
        shutil.rmtree('./index', ignore_errors=True)
        os.mkdir('./index')
    if args.tarPred and args.tarScore:
        shutil.rmtree('./predicted', ignore_errors=True)
        os.mkdir('./predicted')
        print('\nFragments to be indexed and used for TP: %s' % (fragList))

        # Bug found by Uciel Pablo. Moved this line from outside of if 
        # clause to within as this is only run if tarPred is on.
        miRs    = miRinput()
        start   = time.time()###time start
        
        ## Test mode - To find errors ###
        # for i in fragList:
        #    tarFind4(i)
        
        ## Parallel mode
        PP(tarFind4,fragList)
        end     = time.time()
        print ('Target Prediction time: %s' % (round(end-start,2)))
        
        targComb = FileCombine()
        #targComb = 'All.targs' ## Test - open line above when real
        
        start = time.time()###time start
        predTargets = tarParse4(targComb)
        end = time.time()
    
        #print ('Target Prediction time: %s' % (round(end-start,2)))
        
    elif not args.tarPred and args.tarScore:
        targComb = FileCombine()
        #targComb = 'All.targs' ## Test - open line above when real
        
        start = time.time()###time start
        predTargets = tarParse4(targComb)
        end = time.time()
        
        #print ('Target Scoring time: %s' % (round(end-start,2)))
    
    else: ## 
        print("!!Target prediction is OFF - Files in 'predicted' folder might be old!!")
        predTargets = './predicted/All.targs.parsed.csv'
    
    ###Timer
    TFEnd = time.time()
    print ('\n\nTarget prediction time is %s\n\n' % (round(TFEnd-TFStart,2)))
    fh_run.write('Target prediction time is : %s\n' % (round(TFEnd-TFStart,2)))
    #########################################
    
    ## PARE PROCESS AND MAP #################
    PAREStart = time.time()

    ## Deduplicatiom ########################
    ## Needs to be integrated and tag2FASTA will be done in this loop i.e. a collapsed FASTA file written with de-duplicated
    ## file or even if tag-count is upplied - Check how switches will work together
    # if libFormat    == "F":
    #     ### Convert FASTA to Tagcount
    #     ### Sanity check
    #     fh_in       = open(libs[0],'r')
    #     firstline   = fh_in.readline()
    #     if not firstline.startswith('>') and len(firstline.split('\t')) > 1:
    #         print("** File doesn't seems to be in FASTA format")
    #         print("** Please provide correct setting for @libFormat in 'phasworks.set' settings file")
    #         sys.exit()
    #     else:
    #         print("#### Converting FASTA format to counts #######")
    #         dedup_start     = time.time()
            
    #         ## TEST
    #         # newList = []
    #         # for alib in libs:
    #         #     aname = dedup_process(alib)
    #         #     newList.append(aname)
    #         # libs = newList

    #         libs            = PPResults(dedup_process,libs)
    #         # print('Converted libs: %s' % (libs))
    #         dedup_end       = time.time()
    #         fh_run.write("FASTA conversion time:%ss\n" % (round(dedup_end-dedup_start,2)))
        
    # elif libFormat  == "T": 
    #     ### Can be used as-is, check if it is really 
    #     ### Sanity check
    #     fh_in = open(libs[0],'r')
    #     firstline = fh_in.readline()
    #     if firstline.startswith('>'):
    #         print("** File seems tobe in FASTA format")
    #         print("** Please provide correct setting for @libFormat in 'phasworks.set' settings file")
    #         sys.exit()
    #     else:
    #         # print("File seems to be in correct format")
    #         pass

    # else:
    #     print("** Please provide correct setting for @libFormat in 'phasworks.set' settings file")
    #     print("** If sRNA data is in tag count format use 'T' and for FASTA format use 'F' ")
    #     sys.exit()
    
    if args.tag2FASTA:
        shutil.rmtree('./PARE',ignore_errors=True)
        os.mkdir('./PARE')
        PP(tag2FASTA2,args.libs)

    ###########################################
    
    indexFls = [file for file in os.listdir('./index') if file.endswith ('index.1.bt2')]
    print ('These are index files: ',indexFls)
    
    if args.map2DD:
        shutil.rmtree('./dd_map',ignore_errors=True)
        os.mkdir('./dd_map')
        global templib
        for templib in args.libs:
            ## Serial -Test
            #for dd in indexFls:
                #print ('Lib:%s mapped to index: %s' % (templib,dd))
                #mapdd2trans(dd)
            
            ###
            PP(mapdd2trans,indexFls)
    
    ##Timer
    PAREEnd = time.time()
    print ('\n\nPARE processing run time is %s\n\n' % (round(PAREEnd-PAREStart,2)))
    fh_run.write('PARE processing and mapping time is : %s\n' % (round(PAREEnd-PAREStart,2)))
    ## INDEXER ##############################
    ###
    PredStart = time.time()
    
    if args.validate:
        fastaOut = './genomic_seq.fa'
        predTargets = './predicted/All.targs.parsed.csv'
        shutil.rmtree('./PAGe',ignore_errors=True)
        os.mkdir('./PAGe')
        shutil.rmtree('./output', ignore_errors=True)
        os.mkdir('./output')
        validatedTargetsFilename = './output/MPPPValidated.csv'
    
        # 
        targetFinderFile = readFile(predTargets)
        # 
        header = targetFinderFile[0]
        del(targetFinderFile[0])
    
        # 
        global targetFinderList
        targetFinderList = createTargetFinderDataStructure(
            targetFinderFile)
    
        # 
        baseCountsFile = readFile('baseCounts.mem')
        baseCounts = int(baseCountsFile[0])
        baseCountsOffTagLen = int(baseCountsFile[1])
    
        for tagCountFilename in args.libs:
            PAGeIndexDict = {}
            # Variable holding all hits > 2
            allHits = []
            PAGeIndexList = []
            global tagCountFile
            tagCountFile = readFile(tagCountFilename)
            library = os.path.splitext(tagCountFilename)[0]
            PAGeOutputFilename = './PAGe/%s_PAGe' % library
            validatedTargetsFilename = './output/%s_validated' % library
    

            bowtieFiles = [file for file in os.listdir('dd_map') if file.startswith('%s' % library)]
            print("Creating PAGe Index dictionary for lib %s" % library)
            PAGeStart = time.time()
            
            # 
            PAGeIndexAndHits = PPResults(createPAGeIndex, bowtieFiles)
            # 
            for element in PAGeIndexAndHits:
                PAGeIndexDict.update(element[0])
                allHits.extend(element[1])
                PAGeIndexList.append(element[0])

            PAGeEnd = time.time()
            print("PAGe Indexing took %.2f seconds" % (PAGeEnd - PAGeStart))
            fh_run.write("PAGe Indexing took: %.3f seconds\n" % (PAGeEnd-PAGeStart))
    
            # 
            print("Writing PAGeIndex file...")
            PAGeWriteStart = time.time()
            global categoryList
            categoryList = writePAGeFile(PAGeIndexDict, args.genomeFeature,
                allHits, baseCounts, baseCountsOffTagLen, PAGeOutputFilename,
                fastaOut, library)
            PAGeWriteEnd = time.time()
            print("File written. Process took %.2f seconds" % (PAGeWriteEnd - PAGeWriteStart))
            fh_run.write("PAGe index file written. Process took %.2f seconds\n" % (PAGeWriteEnd - PAGeWriteStart))
    
            # 
            print("Finding the validated targets")
            validatedTargetsStart = time.time()
            
            # 
            # 
            validatedTargetsList = PPResults(validatedTargetsFinder, PAGeIndexList)
            validatedTargets = []
            for validatedTarget in validatedTargetsList:
                validatedTargets.extend(validatedTarget)

            validatedTargetsEnd = time.time()
            print("All validated targets found in %.2f seconds" % (validatedTargetsEnd - validatedTargetsStart))
            fh_run.write("All validated targets found in %.2f seconds\n" % (validatedTargetsEnd - validatedTargetsStart))
            
            # 
            # 
            if(validatedTargets):
                # print(validatedTargets)
                print("Writing validated targets")
                writeValidatedTargetsFile(header, validatedTargets,
                    validatedTargetsFilename)

            else:
                print("No targets could be validated for the set of miRNAs "\
                      "in lib %s." % library)

        ## Revmap results
        if args.featureFile: ## User used feature file no reverse mapping
            filetag         = '_validated'
            uniqRevmapped   = resultUniq(filetag) # uniqRevmapped='./output/AllLibValidatedUniq.csv'
        else:
            ReverseMapping()
            filetag         = 'revmapped.csv'
            uniqRevmapped   = resultUniq(filetag) # uniqRevmapped='./output/AllLibValidatedUniq.csv'


    PredEnd = time.time()
    print ('\n\nIndexing and Prediction run time is %s\n\n' % (round(PredEnd-PredStart,2)))
    fh_run.write('Indexing and Prediction run time is : %s seconds \n' % (round(PredEnd-PredStart,2)))
    fh_run.write('Script run time is : %s\n' % (round(PredEnd-FragStart,2)))
    fh_run.close()

#### RUN ##########################################

if __name__ == '__main__':
    nspread = 6
    args.accel = coreReserve(args.accel)
    # if args.accel == 'Y':
    #     args.accel = int(multiprocessing.cpu_count()*0.85)
    # else:
    #     args.accel = int(args.accel)

    
    start = time.time()
    main()
    end = time.time()
    print ("The complete 'sPARTA' run time is %s seconds" % (round(end-start,2)))
    print('The run has completed sucessfully.....CHEERS! - Exiting..\n')
    sys.exit(0)

## Change log
## Moved call to miRinput from outside of target prediction if clause to within to prevent sPARTA from looking for it when tarPred is off
## Unambiguous bases are now calculated during first run only and saved to file baseCounts.mem and pulled from there always durng validate step
## Changed last line to sys.exit(0)

## v1.02 -> v1.03
## Reverse fix in miRferno, miRNA input moved to TF function, implemented frag memory files

## v1.03 -> 1.04
## TarFind3 tweaked, --cat4Filter is now --cat4Show

## v1.04 -> 1.05
## Fixed compatibility to Bowtie 2.4.4 - Tested OK

## v1.06 -> v1.07
## A final report with unique targets from all libraries generated
## Fix: Fixes Bug: 010115 - GFF file lines start with ##
## extractFeature bug updated from MPPP - scaffold fix
## File check and message display to rename genome file without integers
## nthread changed to nspread

## v1.07 -> 1.10 [Major update]
## Two modules updated for improved target prediction and postprocessing to capture more secondary structures
## Sensitivity improved for both Heuristic and Exhaustive modes [But not implemented in this script - Need to be tested for sconfidence score]
## Added a filter whch throws our biologically invalid interactions

## v1.10 -> v1.11rc [stable]
## Fixed a bug introduced in v1.10 - 1) Fixed bugs with args in tarfind4 - Fixed 2) nspread got changed to nthred by mistake - Fixed

## v1.11 -> v1.12
## Fixed a bug introduced due to regresssin, effected fetching seqeunces for reverse strand

## v1.12 -> v1.13
## Fixed a bug where tags weren't being chopped to the specified tag length when creating PAGe index.
## Replaced tagLen variable with minTagLen and maxTagLen in order to add more flexbility with PARE tag lengths.

## v1.13 -> v1.14
## Fixed a bug in which the PARE file name was parsed based on the presence of
## an underscore. Causes problems when there are more than 1 underscore in the 
## filename

## v1.14 -> v1.15
## Fixed a bug causing minTagLen and maxTagLen to be stored as strings rather 
## than intergers if supplied by the user at the command line.

## v1.15 -> v1.16
## Added revFerno: Reverse maps sPARTA results to genome co-ordinates -AK/PP

## v1.16 -> v1.17b
## Modified extractFeature module for better compatibility wit scaffolded genome -AK
## Fixed minor bug where any input fasta file having number in its name would be deleted
## Fixed target prediction switch so that running w/o pare validation gives no error. I manually made chages to master, same as done by Reza in his branch - RH
## 'w' and 'c' set to lower case for all modules
## Cleaned/Organized the code and on-screen messages
## Added a new function to check for libraries before running sPARTA
## Added reverse mapping as default, done only when genomeFile is used to extract genomeFeature, noto in case of featureFile

## v1.17 -> v1.18
## Removed pyfasta requirement
## Fixed bug with revrsemapping for analyses with --featureFile, since it's not possible, reversemapping turned off

## v1.18 -> v1.19
## Removed rpy2 requirement
## Add GTF support, and new switch -annotype GTF/GFF3 -annofile FILENAME
#### Features extracted from GTF can have coordinates differet from GFF3 file, because in GTF we use exons to compute gene coordinates
#### This can lead to difference in number of validated targets (as some features were missing from GTF file) and difference in name of
#### intergenic region (due to absence of a few features the up and dpwn correspond to other gene)
## Added functionality to count the 'N's in fetched seqeunce, before wriring and write those seqeunces which are not purely 'N'. 
#### As these give error in bowtie2 while making index. Also to capture 'N' or 'n's all seqeunces are modified to upper case

## v1.19 - v1.20
## Added annoType annoFile argument in place of gffFile to allow for GTF files and adjusted corresponding argument requirements -RH
## Changed PARE reads/window abundance to only display up 3 decimal places -RH

## v1.20 -> 1.21
## Removed the stray rpy2 imports from main

## v1.21 -> 1.22
## Fixed reverse mapping error (same as revFerno) - No valid triggers were being found at 'c' strand - In reverse mapping 
#### function for crick strand genomic start site was beng used instead of genomic end to compute 
#### new coordinates - Fixed by replacing geno_end    = coord_dict_crick[gene_name][1] with
#### geno_end    = coord_dict_crick[gene_name][2] -  see the difference in coordinates
## Fixed precision bug, which was introduced in v1.21 - fixed by Reza
## Added GTF comaptibility, to make it comaptible with cufflinks experiments

## v1.22 -> v1.23 [Atul]
## Improved reverse mapping. If coords includes chromosome end i.e. '-' as end coordinate then revese mapping 
#### generated an error. Now, a dict of chromo lengths "chromoLenDict" is used to fetch chromosome ends
## Added a module to slightly better assign the cores - Used from PHAS detect

## v1.23 -> v1.25 [Atul] [major fix]
## GFF parser skips empty lines - resolves an issue fron user in Iowa
## Added error check for match/missing identifiers between GFF/GTF and FASTA file
## Fixed - Bug 05-12-17a - This lead to 2-nt shift in extracted seqeunces. The ise was that chromosome seqeunce was
#### being stored as list in the dictionary. In EXTRACTFASTA1, value fetched from disctionary was being converted to 
#### string in single command so the chromo seqeunce has extra " [' " at 5'-end and extra " '] " at 3'-end
## Fixed bug 05-12-17b - where first line of GFF was skipped using readline, sometimes skipping a real entry, lines
#### are now checked for '#' and skipped
## Added a function to check for all input files before executing analysis
## Added checkuser function to pass message
## Added prints and formatted command-line info



## FUTURE DEVELOPMENT [Pending
## Added a de-duplicator to sPARTA - No need to run tally - This needs to be integrated to sPARTA
## Fine tune paralleization in validation part
## Add degradome plots
## Add funtionality to include probabilty of accesibilty and RNA-RNA duplex formation, combine this -value with PARE based p-value
## fasta file is read once to fastaList. What does the 'unambiguousBaseCounter' function does with FASTA file? Can we use the fastaList to avoid file reading again?
## Which switches to use to not perform reverse mapping if feature file is supplied? Rename final uniq file to include 'revammpped' consitant with library specific files. 