#!/usr/local/bin/python3

## phasmerge    : Collapses library specific results to genome-level, generates a summary, 
##                matches with known annotations and compares PHAS summaries 
## updated      : version-1.27 04/04/17
## author       : kakrana@udel.edu, atulkakrana@gmail.com

## Copyright (c): 2016, by University of Delaware
##              Contributor : Atul Kakrana
##              Affilation  : Meyers Lab (Donald Danforth Plant Science Center, St. Louis, MO)
##              License copy: Included and found at https://opensource.org/licenses/Artistic-2.0

import os,glob,sys,difflib,time,shutil,argparse,math,operator
import operator,datetime,subprocess,multiprocessing,re,getpass
from multiprocessing import Process, Queue, Pool
from operator import itemgetter
from itertools import groupby

########### PHASER DEVELOPER SETTINGS ########

setFile         = "phasis.set"
memFile         = "phasmerge.mem"
res_folder      = "summary_%s"   % (datetime.datetime.now().strftime("%m_%d_%H_%M"))      ## Folder with all ther results
comp_folder     = "compare_%s"   % (datetime.datetime.now().strftime("%m_%d_%H_%M"))
cleanup         = 1
cores           = 0

### COLLAPSER DEFAULTS #############
fileType        = 'L'               ## 'L' = *cluster.boundary.without.PARE.validation.list file, N = *NO.by.PARE.file, Y = *YES.by.PARE.file, C: if formatted CSV for intercomparision
region          = 'I'               ## 'G': Genic and 'I': Intergenic - This affects overlap ration in main()

### SUMMARIZER DEFAULTS ############

fetchMax        = 1                                                 ## Fetch the max abundance phasi from tag position summmary, write in separate file
fetchLibAbun    = 0                                                 ## 0: Fetch libwise abundances for tags with 'phase' len required fpr heatmaps and other comparision between loci | 1: All tags from phaster, required to see abundance of phasiRNA tag and tags against all tags in these loci. The latter would be sum of abundances from all libs computed manualy in excel

head            = "Y"
coordsep        ='\t'                                               ## Seprator used in the file
phasedID        = 'N'                                               ## If given than phas ID is used to extract coords, if not than coords specified by user are used
namecol         = 1
pvalcol         = 2
chrcol          = 3                                                 ## As in excel format
startcol        = 4                                                 ## As in excel format
endcol          = 5                                                 ## As in excel format
libcol          = 7

phasiLenFilter  = 'N'                                               ## 'Y' then tags of phase length will be written from the cluster | 'N' - All tags will be written
minAbun         = 1                                                 ## Minimum abundance of tag to be written
matchThres      = 0.99                                              ## Ratio of match required by the cluster to phased loci | For transcripts, to capture biggest cluster like for mapping to IR, use a lower ratio so that longer transcript with smaller match ratio can be included
selfmergeratio  = 0.05                                              ## Default was 0.45. Use a higher value to get lots of the PHAS even if these are overlapping from different p-values.This is good for manual curation of a few loci. Otherwise use a lower value 0.25 to 0.45 to pick one of those overlapping.
startbuff       = 0                                                 ## While extracting sequence through coords2FASTA a buffer is added to start, 
                                                                    ## start position in phased ID has this buffer added, ## so minus this buffer 
                                                                    ## for better matching with real loci

maxTagRatioCut  = 0.90
phasCyclesCut   = 7
minAbunCut      = 320                                               ## Atleast 20 reads per position x 2 (strands) x 8 (pos) = 320


#### ANNOTATION DEFAULTS ############

overlapPerc     = 0.1      ## Minimum percentage of PHAS covered by exon. Keep 0 for finding overlapping genes and 0.20 or above for finding precursors
ntoverlapCutoff = 50       ## Values of minimum overlap in nts. This is required if PHAS overlap is checked for precursors and not overlap with any region of genes. FOr precursors one will expect an overlap of 4 phases or more
annomode        = 3        ## 1: PASA feature (GTF) file | 2: Trinity | 3: Rocket feature (GTF) file or generic gtf file with missing "transcript" entry | 4: Both rocket and trinity feature files | 5: PacBio GTF | 6: Comprehensive transciptome (from PASA) - NOTE: New mode should be registered in overlapchecker function


#### Command Line ##############################
################################################
parser = argparse.ArgumentParser()
# flags = parser.add_argument_group('required arguments') ## Add required arguments to this group and optional to parser


reqflags        = parser.add_argument_group("required arguments")
mergeflags      = parser.add_argument_group("Additional/Optional arguments required for --merge mode") ## Add required arguments to this group and optional to parser
compflags       = parser.add_argument_group("Additional/Optional arguments required for --compare mode") ## Add required arguments to this group and optional to parser

reqflags.add_argument('-mode', default='merge', help= 'merge: For summarizing results from'\
    ' sRNA libraries of choice | compare: To compare PHAS summaries from two different set of libraries.'\
    '[Default: merge mode]',required=True)
reqflags.add_argument('-dir',  default='None', type=str, help='directory from your'\
    ' phasdetect run which need to be summarized. [Compulsory]', required=True)

mergeflags.add_argument('-pval',  default='', type=str, help='pvalue cutoff to'\
    ' filter the phased siRNAs loci or transcipts. [Optional]', required=False)
mergeflags.add_argument('-gtf',  default='', type=str, help='GTF file from genome'\
    ' annotation or transcriptome mapped to genome. GTF file must be formatted using'\
    ' gffread (cufflinks) utility [Optional]. To convert your GFF or GTF file to '\
    ' required format: gffread <my.gff3> -T -o <my.gtf>', required=False)
mergeflags.add_argument('-safesearch',  default='T', type=str, help='Turning this ON'\
    ' filters weakly phased loci based on best k-value, maxtagratio and abundance'\
    ' T: For turned ON (default) | F: Turned OFF.[Optional]', required=False)
mergeflags.add_argument('-debug',  default='T', type=str, help='Turning this ON'\
    ' collects usage data from users directory and sends to the authors.'\
    ' T: For turned ON | F: Turned OFF (default).[Optional]', required=False)

compflags.add_argument('-dir2',  default='None', type=str, help='another directory from'\
    ' your phasmerge run, which need to be compared to first one. [Compulsory]', required=False)

args = parser.parse_args()

#### CHECKS 
if args.mode != "merge" and args.mode != "compare":
    print("Unknown value %s for '-mode' supplied" % (args.mode))
    # print("Please input correct mode 'merge' or 'compare'")
    parser.print_help()
    sys.exit()

if args.mode == "merge":
    if args.dir == None:
        print("\nPlease specify directory from phasmerge analysis using the '-dir' parameter")
        print("To see all requred parameters, run: python3 phasmerge -h\n")
        parser.print_help()
        sys.exit()
    elif args.gtf:
        import sqlite3
    else:
        pass

elif args.mode == "compare":
    if args.dir == None:
        print("\nPlease specify first directory from phasmerge analysis using the '-dir' parameter")
        print("To see all requred parameters, run: python3 collapser -h\n")
        parser.print_help()
        sys.exit()
    elif args.dir2 == None:
        print("\nPlease specify a second directory for the comparison between PHAS summaries, use the '-dir2' parameter")
        print("Like the first directory this must be a reesult from phasmerge analysis")
        print("To see all requred parameters, run: python3 phasmerge -h\n")
        parser.print_help()
        sys.exit()
    else:
        ## All required values supplied by user
        pass

else:
    print("Check your input for '-mode' parameter ")
    parser.print_help()
    sys.exit()

#### CLEANUPS
if args.dir.endswith("/"): ## If user add "/" at end of path, this is removed
    args.dir = args.dir[0:-1]
if args.dir2.endswith("/"):
    args.dir2 = args.dir2[0:-1]

#### COLLAPSER FUNCTIONS ####
#############################

def checkDependency():
    '''Checks for required components on user system'''

    print("\n#### Fn: checkLibs ###########################")
    
    goSignal  = True ### Signal to process is set to true 

    ### Check PYTHON version
    pythonver = sys.version_info[0]
    if int(pythonver) >= 3:
        print("--Python v3.0 or higher          : found")
        pass
    else:
        print("--Python v3.0 or higher          : missing")
        goSignal    = False
        # print("See README for how to INSTALL")

    if args.gtf:
        ### Check sqlite3 module for python3
        if 'sqlite3' in sys.modules:
            print("--sqlite3 (python)               : found")
            pass
        else:
            print("--sqlite3 (python)               : missing")
            goSignal    = False
            # print("See README for how to INSTALL")

        ### Check SQLite installation
        issqlite3 = shutil.which("sqlite3")
        if issqlite3:
            print("--SQLite (v3)                    : found")
            pass
        else:
            print("--SQLite (v3)                    : missing")
            goSignal    = False
            # print("See README for how to INSTALL")


    if goSignal == False:
        if args.gtf:
            print("\n** Please install the missing libraries required for matching the annotations")
            print("**'phasmerge' will work fine without the '-gtf' option")
            print("**'phasmerge' has unmet dependendies and will exit for now\n")
            sys.exit()

        else:
            print("\n** Please install the missing libraries before running the analyses")
            print("**'phasmerge' has unmet dependendies and will exit for now\n")
            sys.exit()

    return None

def readSet(setFile):
    '''
    Read and parse external settings file
    '''

    if os.path.isfile(setFile):
        pass
    else:
        print("---Settings file 'phasis.set' not found in current directory")
        print("---Please copy it to same directory as script and rerun")
        sys.exit()

    print("\n#### Fn: Settings Reader ####################")
    
    fh_in   = open(setFile, 'r')
    setFile = fh_in.readlines()
    fh_in.close()
    
    for line in setFile:
        if line: ## Not empty
            if line.startswith('@'):
                line = line.strip("\n")
                # print(line)
                akey,aval = line.split('=')
                param = akey.strip()
                value = aval.strip()
                # print(param,value)
                
                ## Extract values ######### 

                if param.strip() == '@runType':
                    global runType
                    runType = str(value.strip())
                    # print(runType)
                    if (runType != "G") and (runType != "T") and (runType != "S"):
                        print("Please input correct setting for '@runType' parameter in 'phasis.set' file")
                        print("Script will exit for now\n")
                        sys.exit()
                    else:
                        print('User Input runType               :',runType)

                elif param.strip() == '@index':
                    global index
                    index = str(value.strip())
                    if index:
                        print('User Input index location        :',index)
                    else:
                        print('User Input index location        : None')

                elif param.strip() == '@userLibs':
                    global libs
                    # libs = list(map(str,value.strip().split(',')))
                    libs     = [str(x) for x in value.strip().split(',') if x.strip() != '' ] ## This is my dope...
                    print('User Input Libs:                 :',",".join(libs))

                elif param.strip() == '@reference':
                    global reference
                    reference = str(value.strip())
                    print('User Input reference             :',reference)

                elif param.strip() == '@phase':
                    global phase
                    phase = int(value.strip())
                    print('User Input for phase length      :',phase)
                
                elif param.strip() == '@libFormat':
                    global libFormat
                    libFormat = str(value.strip())
                    if (libFormat != "T") and (libFormat != "F"):
                        print("Please input correct setting for '@libFormat' parameter in 'phasis.set' file")
                        print("Script will exit for now\n")
                        sys.exit()
                    else:
                        print('User Input for library type      :',libFormat) 

            else:
                #print("Missed line:",line)
                pass
    # sys.exit()
    return libs

def pvaluereader():
    '''
    Get the best matching p-value to user cutoff
    '''
    print("\n#### Fn: pvalue capture ###################")

    ### 
    if os.path.isdir(args.dir): ## Check to see its a file from bowtie and not tophat mapped folder - Untested
        # print("Directory from phaser analysis found")
        pass
    else:
        print("Specified directory not found: %s" % (args.dir))
        print("Please confirm the directory path, or")
        print("confirm 'phasdetect' run finished successfully")
        print("To see all required parameters, run: python3 phasmerge -h\n")
        print("Script exiting...\n")
        sys.exit()

    ### Read the list result files
    clustL      = [file for file in os.listdir("%s" % (args.dir)) if file.endswith ('.cluster')]
    pvalS       = set()

    # if not pvalS:
    #     print("No *.cluster files found in specified directory:%s" % (args.dir))
    #     print("please check")

    for aclust in clustL:
        # print(aclust,aclust.rsplit("_",4)[1])
        pval = float(aclust.rsplit("_",4)[1].replace('p',''))
        pvalS.add(pval)

    # print("--P-values cached",pvalS)

    ## Sanity check #######
    if not pvalS:
        print("No 'clusters' found - please confirm 'phasdetect' run finished successfully")
        print("Otherwise there are no PHAS clusters in your data")
        sys.exit()
    else:
        pass

    ### Sort pval list ####
    pval_sorted = sorted(pvalS,reverse=True)
    print("Total %s clusters cached at different p-values: %s" % (len(clustL),", ".join(str(x) for x in pval_sorted)))
    # bestpval = percentile(pval_sorted,25)
    # print("-- Best P-value: %s"% (bestpval))

    ### Use the specified or choose best bsed on confidence levels
    pcutoff = None ### Set to None for checks
    if args.pval:
        ### User specified a p-value see if its in list or other better one
        pcutoffL = [] ## List of valid cutoffs
        for apval in pval_sorted:
            # print(apval,float(args.pval))
            if apval            <= float(args.pval):
                pcutoff         = float(apval)
                print("User specified p-value being used:%s" % (pcutoff))
                break
            else:
                pass

        if pcutoff == None:
            ### Bad luck no matching or better p-value found
            print("\nNo PHAS cluster at or above user specified p-val cutoff %s" % (args.pval))
            print("Choose a lower value from these options: %s" % (", ".join(str(x) for x in pval_sorted)))
            print("Alternatively, you can run 'phasmerge' without the '-pval' switch, phasTER will try to find best pvalue for analysis")
            sys.exit()

        else:
            pass


    else:
        ### User did not specified optimize for best
        pcutoff = float(percentile(pval_sorted,25))

        if pcutoff <= float(1e-05):
            ### seems decent use this one
            print("Best guessed p-value being used: %s" % (pcutoff))
            print("User can specify a lower or higher p-val cutoff using '-pval' option")
            pass

        else:
            ### Optimized p-value is lower then expected cutoff
            pcutoff = float(percentile(pval_sorted,5))
            print("** PHAS cluster identified in current data are below recommended confidence threshold (p-val < 1e-05)")
            print("** To continue, please specify a lower p-value using '-pval' switch")
            print("** Possible values for '-pval' switch are : %s" % (", ".join(str(x) for x in pval_sorted)))
            print("** Recommended value for '-pval' switch is: %s" % (pcutoff))
            print("** Recommended command: python3 phasmerge -dir %s -pval %s" % (args.dir,pcutoff))
            sys.exit()

    # sys.exit()
    return pcutoff,pval_sorted

def prepare(pcutoff,libs,res_folder):

    print("\n#### Fn: prepare ############################")
    
    ### Strip libs extension for matching ##############
    libs_name = []
    for i in libs:
        alib = i.rsplit(".",1)[0]
        libs_name.append(alib)

    ### Make a new folder ##############################
    ####################################################
    temp_folder = "./%s/%s" % (res_folder,"temp")
    print("WARNING: Temporary files exits from earlier run, these will be deleted")
    shutil.rmtree("%s" % (res_folder),ignore_errors=True)
    os.mkdir("%s" % (res_folder))
    os.mkdir("%s" % (temp_folder))

   ### Get files for collapser #########################
    listL      = [file for file in os.listdir("%s" % (args.dir)) if file.endswith ('.list')]
    # print(listL)
    print("%s list files found" % (len(listL)))
    
    acount = 0 ### Count of library results copied to collapsed folder
    for alistF in listL:
        alib    = alistF.rsplit(".",7)[0]
        # print("alib:",alib)
        if alib in libs_name:
            # print("copied for summarization:%s" % (alistF))
            acount += 1
            afile   = "%s/%s" % (args.dir,alistF)
            shutil.copy(afile,temp_folder)
    # sys.exit()

    ### Get files for summarizer ######################
    ###################################################
    clustL      = [file for file in os.listdir("%s" % (args.dir)) if file.endswith ('.cluster')]
    combL       = []    ## List of file to be combined
    bcount      = 0     ## Counter for cluster files
    # print(clustL)
    print("%s clust files found" % (len(clustL)))
    for aclust in clustL:
        # print("Clust File:",aclust)
        info    = aclust.rsplit("_",4)
        # print("-Info:",info)
        
        if len(info)==5: ### If there is some other file or contaminant then it is not considered, need better way to filter those out
            alib    = info[0].rsplit('.',2)[0]
            # print("-alib:",alib)
            pval    = info[1].replace("p","")
            aphase  = info[3]

            # print("-%s,%s,%s" % (alib,pval,aphase))
            if (alib in libs_name) and (float(pval) <= float(pcutoff)) and (str(aphase) == str(phase)): ## phase added to ensure that if user runs both 21 and 24-nt PHAS in same directory, the files from different analysis are not picked up. 
            ### V1.14 and above all cluster file for p-values less then the used are picked up, because collapser uses
            ### files from all diffrent confodence levels. Also the matchThres for cluster and PHAS is increased to 0.99
            ### because now files for all confidence are available so cluster should match at high cutoff, reducing number
            ### of total clusters picked up
                # print("--Cluster file found:%s" % (alib ))
                # print("----Filename:%s | Lib:%s | pval = %s" % (aclust,libs,pval))
                afile = "%s/%s" % (args.dir,aclust)
                combL.append(afile)
                bcount +=1
    # print("--%s cluster files captured for concatanation" % (len(combL)))
    # print("There are files:",( [x.rpartition('/')[-1] for x in combL]))

    ### Concatanate these files #######################
    ###################################################
    aname       = "%s/ALL.%sPHAS_p%s_srna.cluster" % (res_folder,aphase,pcutoff) ### Out file name
    clustfile   = FileCombine(combL,aname)

    ### Sanity Check ##################################
    ###################################################
    if acount == 0 or bcount == 0 :
        print("** No phasdetect results detected")
        print("** Check if periodicity i.e. 21nt, 22nt or 24nt from phasdetect analysis")
        print("   and periodicity mentioned in 'phasis.set' is same")
        sys.exit()
    else:
        pass
    print("--%s files prepared for collapsing" % (acount))

    print("--Working folder:%s | Temporary Folder:%s\n" % (res_folder,temp_folder))
    # sys.exit()

    return temp_folder,clustfile

def collapser(temp_folder,fileType):
    '''
    Master function that parses results files, makes a dictionary, checks overlap 
    '''
    if fileType == 'L':
        print('\nList files selected for analysis - Converting them to readable format')
        fls     = glob.glob(r'./%s/*.PARE.validation.list' % (temp_folder))
        # print ('Here are the files that will be converted:',fls,'\n')
        print ('Total files to analyze: %s' % (len(fls)))

        ### Prepare first file for comaprision
        firstfile = fls[0]
        firstdict = dictmaker(firstfile)
        sys.exit()
        
        ### Start the comaprision
        compareflag == False
        for fl in fls[1:]:
            ## Prepare chr/scaffold/transcrpt specific dict of PHAS

            if compareflag == False:
                ## No cpmarision made yet this is the first and uses list directly from files
                adict       = dictmaker(fl)
                rawinputs   = inputmaker(firstdict,adict)


    return None

def mergePHAS(aninput):
    '''
    Takes two chr/scaffold and transcriptome specific lists and removes redundant PHAS - written to paralleized the process, and to reduce the numer of matchings required
    '''
    print ("\n#### Fn: PHAS collapser #####################")

    alist,blist = aninput
    # print("Alist",alist)
    # print("Blist",blist)
    main_dict   = {}    ## declare empty dictionary
    tmp_list    = []    ## List to hold all co-odinates before making a dictionary based on chromosme
    neg_list    = []    ## List to store keys that needs to be removed
       
    ### First List - Fill up the dictionary with entries from first file
    ####################################################################
    alib        = "nd" ## Not determined in latest version of collapser v1.13 onwards, kept for future addition
    blib        = "nd"
    
    acount      =1
    tmp_dict    = {}    ## Dictionary to store values for one file - recycled after every file
    for anent in alist:
        # print(anent)
        key     = 'a-%s-%s-%s' % (anent[2],anent[3],anent[4])    ### 'a' added to diffrentiate the key from PHAS in different file that have same coordinates. Key is 'a',chrid, start and end makes a key - 01/03/2017
        chrid   = anent[2]
        start   = int(anent[3])
        end     = int(anent[4]) ## 1 added because when opening a range using start and end, end number is not included in range - - Critical bug fixed in v4->v5 and later regressed/removed in v9->v10

        value = ((chrid,start,end),anent,alib,key) ## Added key in value so that we don't need to make it later
        # print("Key",key,"| Value",value)
        tmp_dict[key] = value
        acount+=1
    
    ### Prepare dict for comaprision
    main_dict.update(tmp_dict) ## Update the main dict if this is the first file
    # xdict   = main_dict.copy()
    # xcount  = len(main_dict)

    ### Start comparision with second file
    ######################################
    bcount      = 0     ## Count entries in bfile
    tmp_dict    = {}    ## Dictionary to store values for one file - recycled after every file
    if len(main_dict) != 0:
        ### This check is helpful for transcriptome where the PHAS might be absent in either library
        for ent_splt in blist:
            # print(ent_splt)
            bcount          +=1
            ratiodict       = {}                                ## Dictonary to hold ratios of comparision of this entry with all in dictionary

            #### Compare with dict entries and get ratio
            new_chrid   = ent_splt[2]                           ## Chromosome or transcript name
            new_start   = int(ent_splt[3])
            new_end     = int(ent_splt[4])                      ## 1 added because when opening a range using start and end, end number is not included in range - Critical bug fixed in v4->v5 and later regressed/removed in v9->v10
            newRegion   = list(range(new_start,new_end))
            newKey      = 'b-%s-%s-%s' % (ent_splt[2],ent_splt[3],ent_splt[4]) ### 'b' added to diffrentiate the key from PHAS in different file that have same coordinates. Key is 'b',chrid, start and end makes a key - 01/03/2017
            newValue    = ((new_chrid,new_start,new_end),ent_splt,blib,newKey) 
            
            #### Find a match for this PHAS in first file
            for i in main_dict.values():                         ## Compare with all dictionary values
                #print(i)
                existKey    = i[3]
                exist_chrid = i[0][0]
                exist_start = i[0][1]
                exist_end   = i[0][2]
                if new_chrid == exist_chrid:                     ## Check if chr or transcript is same
                    existRegion = list(range(exist_start,exist_end))
                    sm          = difflib.SequenceMatcher(None,existRegion,newRegion)
                    ratiodict[str(existKey)]=round(sm.ratio(),2) ## Make a dict of main dict entries and their comparision ratio with current new entry

                else:
                    ### None found
                    ratiodict[str(existKey)]=round(0.00,2)       ## None of the existing entry matches with the current ones
                    pass

            
            #### Decide if entry is different enough to be added
            ####################################################
            
            existKey = max(ratiodict,key=ratiodict.get)          ## Key from main_dict with max comparable ratio for current entry
            maxratio = ratiodict[existKey]                       ## Max ratio with the earlier phased locus
            # print("\n\n#####################")
            # print(existKey.strip(),maxratio)                     ## If maxratio is zero same entry will appear again here
            
            if maxratio <= overlapCutoff:
                ## Overlap is less then cutoff, treat as new loci, no need to delete any entry from existing PHAS
                tmp_dict[newKey]=newValue                       ## Life = one file
                # print('Adding new loci')

            elif maxratio > overlapCutoff: 
                ## Choose the longest loci between the matching two
                # print('Selecting longest loci')
                # print('Remade',existKey, 'New Key',newKey,)
                uniqid,achrid,astart,aend = existKey.rsplit("-",3)     ## Start and end extracted from existKey because exist_start and exist_end belonged to last ent in dictionary and not the one with max overlaping existKey
                existLen    = (int(aend)-int(astart)+1)
                newLen      = (int(new_end)-int(new_start)+1)
                # print ('Length of existing:%s | Length of new:%s' % (existLen,newLen) )
                
                if newLen > existLen : ### New Loci is longer
                    # print ('New phased loci is longer')
                    neg_list.append(existKey)
                    tmp_dict[newKey]=newValue
                
                else: ## The loci in dictionary is longer
                    # print('Existing phased loci is longer or equal to new - No updates made')
                    pass
            else:
                # print('Redundant')
                pass
    else:
        ## No PHAS in alist so all PHAS in blist should be added to dictionary
        for ent_splt in blist:
            # print(ent_splt)
            bcount          +=1
            new_chrid       = ent_splt[2]                           ## Chromosome or transcript name
            new_start       = int(ent_splt[3])
            new_end         = int(ent_splt[4])                      ## 1 added because when opening a range using start and end, end number is not included in range - Critical bug fixed in v4->v5 and later regressed/removed in v9->v10
            newKey          = 'b-%s-%s-%s' % (ent_splt[2],ent_splt[3],ent_splt[4]) ### 'b' added to diffrentiate the key from PHAS in different file that have same coordinates. Key is 'b',chrid, start and end makes a key - 01/03/2017
            newValue        = ((new_chrid,new_start,new_end),ent_splt,blib,newKey)
            tmp_dict[newKey]=newValue
            chrid          =new_chrid                              ## Alist was empty to chr_id to use below in print statement would be mepty, fill it here 

    
    main_dict.update(tmp_dict) ### Update the main dict with new non-overlapping PHAS
    # ydict   = main_dict.copy()
    # ycount  = len(main_dict)
        
    ######################## Test ####################
    ##print ('Dictionary')
    #fh_test = open('keysTest','w')
    #for i in main_dict.keys():
    #    #print (i)
    #    fh_test.write('%s\n' % i)
    #print ('\nLength of dictionary: %s' % (len(main_dict)))
    #
    #
    #fh_test2 = open('NegKeytest', 'w')
    #for i in neg_list:
    #    fh_test2.write('%s\n' % i)
    #print ('Length of negative list: %s' % (len(neg_list)))
    ################################################
    
    #### Remove keys in negative list
    # print("\n### Merging complete - Now removing the redundant entries")
    # print(blist)
    neg_set = set(neg_list) ## This will remove redudant entry picked from different confidence levels
    for akey in neg_set:
        # print ("To be removed:%s" % (akey))
        try:
            del main_dict[akey]
            # print (akey, '\nKey found in main dict and is being removed')
        except KeyError:
            # print (akey, '\nKey not found') ### This could be many if p-value cutoff includes different confidence levels (<=)
                                            ### Suppose bfile has same PHAS loci from different p-vals then the first instance in neg_list will seek and delete matching PHAS but other instances (of different p-val) will not find that key - tested OK
            pass


    print ("### Number of merged phased loci for %s: %s" % (chrid,len(main_dict)))

    # print("\n*** FIRST DICT", [x for x in xdict.keys()])
    # print("*** SECOND DICT",[y for y in ydict.keys()])
    # print("\n*** NEG SET",neg_set)
    # print("*** NEG LIST",neg_list)
    # print("\n###### First File: chr:%s | entries:%s" % (chrid,acount))
    # print("###### Second File: chr:%s | entries:%s" % (new_chrid,bcount))
    # print("###### Before collapsed entries:%s | Collapsed entries:%s | After removal of neg set:%s" % (xcount,ycount,len(main_dict)))

    return main_dict

def writer_collapse(collapsedL,pcutoff):
    '''
    Writes collapsed results - updated with v1.13 to handle chr/scaffold or trascript wise grouped results
    '''

    print ("\n#### Fn: Collapsed Writer ############")

    outfile1    = "./%s/%sPHAS_p%s_collapsed.txt" %     (res_folder,phase,pcutoff)
    fh_out1     = open(outfile1,'w')
    fh_out1.write('Name\tp-val\tChr\tStart\tEnd\tStrand\tLib\n')
    outfile2    = "./%s/%sPHAS_p%s_collapsed.list" %    (res_folder,phase,pcutoff)
    fh_out2     = open(outfile2,'w')### For our Genome viewer, No header required

   
    ### Default list mode - updated in v1.13 
    if fileType == "L":

        finalphasL      = []    ## List to store cleaner PHAS results
        acount          = 0     ## Count of merged PHAS
        for achr in collapsedL: ## achr is a dict with key: chr-start-end and values: [((chr,start,end),(phase,p-val,chr,start,end),filename,chr-start-end),(....)]
            if achr:                 ## Not empty
                phasL   = [x for x in achr.values()]
                for aphas in phasL:
                    # print(aphas)
                    apval       = aphas[1][1]
                    achr        = aphas[1][2]
                    astart      = aphas[1][3]
                    aend        = aphas[1][4]
                    astrand     = "NONE"
                    alib        = aphas[2]    
                    finalphasL.append((apval,achr,astart,aend,astrand,alib))
                    acount+=1

        # print("%s results prepared for writing" % (acount))

        ## Sort results on number of chromsomes/scaffolds or transcripts
        finalphasL_s    = sorted(finalphasL, key=lambda x: int(x[2]), reverse=False ) ## on start site
        finalphasL_s2   = sorted(finalphasL_s, key=lambda x: '{0:0>8}'.format(x[1]).lower()) ## sort on chr/scaffold or trans
        
        ## write results
        anum    = 1    ## For naming PHAS
        for ares in finalphasL_s2:
            # print(ares)
            fh_out1.write('Phas-%s\t%s\n' % (anum,"\t".join(ares)))
            fh_out2.write('%s.%s.%s\n' % (ares[1],ares[2],ares[3]))
            anum +=1

    
    ### Compare mode files not updated in v1.13
    elif fileType == 'C':
        fh_out1.write('Name\tp-val\tChr\tStart\tEnd\tStrand\tLib\tfileName\n')
        for value in main_dict.values():
            print(value)
            #print('Phas-%s\t%s\t%s\t%s\t%s\t%s\n' % (anum,value[1][1],value[1][2],value[1][3],value[1][4],value[1][6]))
            fh_out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (value[1][0],value[1][1],value[1][2],value[1][3],value[1][4],value[1][5],value[1][6],value[2]))
            fh_out2.write('%s.%s.%s\n' % (value[1][2],value[1][3],value[1][4]))
            # fh_out3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (value[1][0],value[1][1],value[1][2],value[1][3],value[1][4],value[1][5],value[1][6]))

    
    else:
        print("Unknown 'fileType' - Developer's error")
        print("Please report issue here: https://github.com/atulkakrana/phasTER/issues")
        sys.exit()

    print("Total collapsed PHAS after merging:%s" % (len(finalphasL)))
            
    fh_out1.close()
    fh_out2.close()
    # fh_out3.close()
    return outfile1,outfile2

def listConverter(afile,pcutoff):
    '''
    Parses phaster-core results file
    '''

    print("\n#### Fn: listConvertor ####################")

    fh_in   = open(afile,'r')
    fh_out  = open('%s_converted.list' % (afile.rpartition('.')[0]),'w')
    entries = fh_in.readlines()

    #### Prepare PHAS list
    alist   = []    ## List to hold PHAS entries
    for i in entries:
        if i.strip(): ## Remove an empty line from end file
            ent_splt         = i.strip('\n').split('=')
            pval,phase,trash = ent_splt[0].strip().split('|')
            chromo_start,end = ent_splt[1].strip().split('..')
            chromo,sep,start = chromo_start.rpartition(':')
            if float(pval) <= float(pcutoff):
                if runType      == 'G':
                    fh_out.write('%s\t%s\t%s\t%s\t%s\tNONE\tNONE\n' % (phase,pval,chromo.strip(),start,end)) ##Chromosome has space before it which later gives error while key matching
                    alist.append((phase,pval,chromo.strip(),start,end))
                elif runType    == 'T' or runType == 'S': ## Header has lots of stuff
                    fh_out.write('%s\t%s\t%s\t%s\t%s\tNONE\tNONE\n' % (phase,pval,chromo.strip(),start,end)) ##Chromosome has space before it which later gives error while key matching
                    alist.append((phase,pval,chromo.strip(),start,end))
                else:
                    pass
    print("%s elements in PHAS list from %s file" % (len(alist),afile))
    # print("%s list" % (afile), alist)
    
    fh_in.close()
    fh_out.close()
    
    return alist

def groupPHAS(PHASlist):
    '''
    Prepares a chr/scaffold and transcript specific dictionary from a PHAS list (mode 0) or PHAS merge results (mode 1)
    '''
    print("\n#### Fn: groupPHAS ####################")


    #### Prepare grouped list from a list of PHAS corresponding to the chr/scaffolda and transcript
    adict       = {} ### Grouped dict
    groupL      = [] ### List to store grouped PHAS
    PHASlist_s  = sorted(PHASlist,key = itemgetter(2))
    agrps       = groupby(PHASlist_s, itemgetter(2))

    acount      = 0 ### Counts uniqe groups i.e. chr/scaffolds/transcripts
    for (key, data) in agrps:
        # print("KEY:",key)
        # print("PHAS",data)
        dataL  = [] ### group specific PHAS
        for i in data:
            dataL.append(i)
        groupL.append((key,dataL))
        acount+=1

    ### For test
    # for i in groupL:
    #     print("KEY",i[0])
    #     print("DATA",i[1])
        
    print("Groups identified:%s" % (acount))
    # print("Grouped",groupL)


    return groupL

def selfMerge(dictitems):
    '''
    For each library remove redudant loci either from dfferent p-values or lengths. Input is chr/scaffold and transcript specific list
    '''
    print("\n#### Fn: selfMerge ###########################")
    chrkey,chrvalL  = dictitems ### Value is grouped list from python
    nonredundantL   = [] ### Final PHAS list
    processedL      = [] ### Keys for identifiers either processed or redundant

    print("++ Merging %s PHAS for chr:%s" % (len(chrvalL),chrkey))
    # print("and values",chrvalL)
    for aphas in chrvalL:
        aname,apval,achr,astart,aend = aphas
        akey    = "%s-%s-%s" % (astart,aend,apval) ### Assumption is that key will be unique due to pval, even if cordinates are same
        aregion = list(range(int(astart),int(aend)))
        tempD   = {} ### Stores key and full PHAS entry for overlapping PHAS to write results
        tempL   = [] ### Stores p-values and length for overlapping PHAS to make a decision
        
        if akey not in processedL:
            ### This PHAS is not compared yet or this didn't overlapped with anyone yet
            processedL.append(akey) ### MArk this processed
            # print("\n--To be merged:",aphas)

            ### Find overlapping PHAS
            ########################
            for bphas in chrvalL:
                bname,bpval,bchr,bstart,bend = bphas
                bregion = list(range(int(bstart),int(bend)))
                sm      = difflib.SequenceMatcher(None,aregion,bregion)
                aratio  = round(sm.ratio(),2)

                if aratio > selfmergeratio:
                    ## Some overlap found cache the PHAS and other parameters
                    bkey        = "%s-%s-%s" % (bstart,bend,bpval)
                    blength     = int(bend)-int(bstart)
                    tempD[bkey] = bphas
                    tempL.append((bkey,bpval,blength)) 


            ### Choose the best candidate
            #############################
            if tempL:
                # print("--Overlapped",tempL)
                ### There were a few overlapping PHAS
                # print("%s overlaps found" % (len(tempL)))
                # print("Overlaps",tempL)
                tempL_s = sorted(tempL, key=lambda x: float(x[1]), reverse=False)
                
                ### Choose on best p-value
                if float(tempL_s[0][1]) < float(tempL_s[-1][1]):
                    # print("Sorted on p-val",tempL_s)
                    # print("Best candidate selected on p-val (%s)" % ((tempL_s[0][1])))
                    ### There were more then one overlap and btter p-val exists
                    ### If the ent itself is best then it will be added to final 
                    ### output otherwise others will be chosen
                    finalkey = tempL_s[0][0]
                    finalval = tempD[finalkey]
                    nonredundantL.append(finalval)

                    ### Trash the rest of keys i.e. these are never checked again for overlap
                    for i in tempL_s[1:]:
                        processedL.append(i[0])

                else:
                    ### Choose the one with longest length. This works even if there is just one 
                    ### overlap with self
                    tempL_s = sorted(tempL, key=lambda x: int(x[2]), reverse=True)
                    # print("Sorted on lengths",tempL_s)
                    # print("Best candidate selected on length (%s nt)" % (tempL_s[0][2]))
                    finalkey = tempL_s[0][0]
                    finalval = tempD[finalkey]
                    nonredundantL.append(finalval)

                    ### Trash the rest of keys i.e. these are never checked again for overlap
                    for i in tempL_s[1:]:
                        processedL.append(i[0])

            else:
                ### No  overlap for query PHAS
                # print("No overlap - Unique PHAS")
                nonredundantL.append(aphas) ### Included in results
                continue
        else:
            # print("Compared already, best candidate cached")
            pass

    # print("#### Self Merged",nonredundantL)
    # print("Group:%s | Total PHAS:%s | Nonredundant PHAS:%s" % (chrkey,len(chrvalL),len(nonredundantL)))
    # sys.exit()

    return chrkey,nonredundantL

def listTocollapsed(alist):
    '''
    Converts non-redundant list of PHAS from first file to 
    a faux collapsed list for writing results in case there is just one file
    and no comaprision will be done
    '''
    
    collapsedL = [] ## List of dicts and format is key: chr-start-end and values: [((chr,start,end),(phase,p-val,chr,start,end),filename,chr-start-end),(....)]

    acount = 0      ## Total PHAS from all chr/scaffolds and trans
    for achr,avals in alist:
        ## Prepare a dict for this chr/scaffold or trans
        chrD = {} ## Hold resuts for specific chr/scaffold or transcriptome
        for aphas in avals:
            acount +=1
            achr    = aphas[2]
            astart  = aphas[3]
            aend    = aphas[4]
            afile   = "nd"
            akey    = "%s-%s-%s" % (achr,astart,aend)
            # print(akey)
            aval    = ((achr,astart,aend),(aphas),afile,akey)
            # print(aval)
            # collapsedL.append(aval)
            chrD[akey] = aval

        ## Add to a list as results from different chr/scaffolds and transcriptome are reported by PPresults
        collapsedL.append(chrD)

    print("Groups in List:%s | Groups in collapsed format list:%s" % (len(alist),len(collapsedL)))
    print("Total PHAS:%s" % (acount))

    return collapsedL

def collapsedToDict(alist):
    '''
    Takes a collapsed list and converts to a dict format as required for next comaprision
    '''

    print("\n#### Fn: collapsedToDict #####################")
    ### The list of results for different chr/scaffold an transcript has merged PHAS results in form of dict [(chr 1 results as dict),(chr 2 results as dict)....]
    ### Where key: chr-start-end and values: [((chr,start,end),(phase,p-val,chr,start,end),filename,chr-start-end),(....)]
    ### The dicts are already grouped based on chr/scaffold and transcript unlike mode=0
    acount  = 0     ### Counts uniqe groups i.e. chr/scaffolds/transcripts
    adict   = {}    ### Grouped dict like mode 0
    for achr in alist:
        ### achr is dict for specific chr/scafold and trans from PHAS merge
        if achr: ### It could be empty
            # print("## Collapsed list values")
            # print("\n",achr.values())
            phasL   = [x for x in achr.values()]
            # print("Len of phasL:%s " % len(phasL))
            key     = phasL[0][0][0]    ### Steal the chromosme/scaffold or transcriptome id from first entry
            groupL  = []                ### This will store extracted info from every PHAS entry
            for ent in phasL:
                aphas = ent[1]          ### ent[1] is ('24', '1e-05', '10', '11746264', '11746693')
                groupL.append(aphas)

            # print("\n Mode1 groupL:",groupL) ## TO check if structure is same as mode 0
            adict[key] = groupL
            acount+=1
    print("Length of collapsed dict:%s | and unique entities:%s" % (len(adict),acount))
    # print("Keys:",(adict.keys()))
    # print("Values:",(adict['1']))

    return adict

def inputMaker(adict,bdict):
    '''
    prepares a chr/scaffold and transcript-wise list from set of two collapsed dicts. Last entry of chr/scaffold and transcriptome speciifc list is the filename
    '''
    print ("\n#### Fn: inputMaker #########################")
    ### Get unique keys i.e. chr/scaffolds and transcripts from both lists
    rawinputs   = []    ## This will hold matching tuples i.e. for a unique chr/scaffold or trascriptome that will be used for comaprisions
    keyset      = set() ## a set of unique keys from both dicts   

    # print("Keys:",(adict.keys()))
    # print("Values:",(adict['1']))

    ## Trim filenames
    # # print(afile)
    # aname = afile.rpartition('/')[-1].rsplit(".",6)[0] ### Filename 2600_chopped.txt.cluster.boundary.without.PARE.validation.list
    # bname = bfile.rpartition('/')[-1].rsplit(".",6)[0]

    ### Get unique keys
    akeys = adict.keys()
    bkeys = bdict.keys()
    
    for x in akeys:
        keyset.add(x)
    for y in bkeys:
        keyset.add(y)

    print("Total unique keys from both lists:%s" % (len(keyset)))
    # print(keyset)

    acount = 0 ## Count of unique identifiers for which PHAS found
    bcount = 0 ## Count of unique identifiers for which PHAS found
    for key in keyset:
        # print ("Preparing inputs for merging for key:%s" % (key))
        
        #### Values from alist
        try:
            avals = adict[key]
            # print (key, '\nKey found in adict and its PHAS fetched')
            acount+=1
        except KeyError:
            # print (key, '\nKey not found')
            avals = []    ## Temp. lists for values from specific chr/scaffold or trascriptome 
            pass

        #### Values from blist
        try:
            bvals = bdict[key]
            # print (key, '\nKey found in bdict and its PHAS fetched')
            bcount+=1
        except KeyError:
            # print (key, '\nKey not found')
            bvals = []    ## Temp. lists for values from specific chr/scaffold or trascriptome
            pass

        ### Add these chr/scaffold or tramscriptome specific PHAS to rawinputs
        # avals.append(aname)
        # bvals.append(bname)
        rawinputs.append((avals,bvals))
        # print(rawinputs[0])
        # sys.exit()


    # print("Unique identifiers in rawinputs:%s" % (len(rawinputs)))
    print("Unique identifiers for PHAS in A:%s | B:%s" % (acount,bcount))

    return rawinputs

#### SUMMARIZER #############
#############################

def PHASreader(coordsfile):
    '''
    Reads coordinates file and return a list of PHASE ID and their coords values
    '''

    print("\n#### Fn: PHAS Reader #################")

    fh_in = open(coordsfile)
    if head == 'Y':
        phashead = fh_in.readline()

    phasCount   = 0       ## Count of entries read
    phasList    = []       ## List to store results 
    for ent in fh_in:

        # print('\n\nEntry from input file being matched: %s' % (ent))
        coords = ent.split(coordsep)
        #print('coords:',coords)
        if phasedID == 'Y':                                                         ## Extract coords from phased ID
            loci = coords[0].split('_')                                             ## Phas-10_w_3:27574117:27574772
            # phasID = loci[0]
            fusedcoords     = loci[2].split(':')                                    ## 3:27574117:27574772
            aname           = coords[namecol-1]
            apval           = coords[pvalcol-1]
            alib            = coords[libcol-1].strip('\n')
            get_chr_id      = fusedcoords[0]
            astart          = fusedcoords[1]
            aend            = fusedcoords[2]
            get_start       = int(astart)+startbuff                         ## It should be added to reduce length of loci
            get_end         = int(aend)+1 
            phasID          = '%s_%s_%s' % (get_chr_id,astart,aend)           ## 1 added because when opening a range using start and end, end number is not included in range
            phasCount       += 1
        
        elif phasedID == 'N': ## User specified coords
            # print("File with user specified columns will be used")
            aname           = coords[namecol-1]
            apval           = coords[pvalcol-1]
            alib            = coords[libcol-1].strip('\n')
            get_chr_id      = coords[chrcol-1]
            astart          = coords[startcol-1]
            aend            = coords[endcol-1] 
            get_start       = int(int(astart)+startbuff)
            get_end         = int(aend)+1                                       ## 1 added because when opening a range using start and end, end number is not included in range
            phasID          = '%s_%s_%s' % (get_chr_id,astart,aend)             ## Name to be used in output file
            # print("Phased Loci: %s #######################################################" % (phasID))
            phasCount       += 1

        else:
            print("Please input correct value for 'phasedID' or check your file")

        get_value  = (list(range(int(str(get_start)),int(str(get_end)))))
        phasList.append((aname,apval,alib,phasID,get_chr_id,get_start,get_end,get_value))
        
    print("Entries read:%s | Entries cached:%s" % (phasCount,len(phasList)))

    return phasList,phashead

def getClust(aninput):

    # print ("\n#### Fn: Cluster Search #########################")
    
    ent,clusters,nphas,totalphas   = aninput   ## Input for parallelized instance
    entL            = []        ## This will store three sublists as a single result for multiprocesing and decoded later
    # resList         = []        ## Store final results as (phas,[(phasiRNA),(PhasiRNA)],[extra info])
    # resList2        = []        ## Store phasiRNAs from all clusters as (phas,[(phasiRNA),(PhasiRNA)],[extra info])
    finalClustL     = []        ## Store matching clusters
    phasCount       = 0         ## Total phased loci in file
    uniqMatchCount  = 0         ## Atleast one cluster present for one phased loci
    allMatchCount   = 0         ## Total number of matched cluster
                                                            
    # print(ent)
    aname,apval,alib,phasID,get_chr_id,get_start,get_end,get_value = ent
    # print("This is the PhasId: %s | values:%s" % (phasID,get_value))
    # print("\nPhaseID being queried for cluster:%s ##############" % (phasID))
    phasCount +=1 
    print("%s/%s clusters cached" % (nphas,totalphas))

    ### Find matching cluster
    matchCount          = 0         ## Total maching clusters for a phased loci - if same cluster in multiple libraries
    finalMatchList      = []        ## Holds best cluster, from multiple libraries
    tempAllList         = []        ## Hold phasiRNAs from all matching clusters, of use for phased 
                                    ## transcripts to capture allphasiRNAs, must be used with low matchThres
    for aclust in clusters[1:]:
        tempMatchList   = []        ## To hold results of current matching cluster
        aclust_splt     = aclust.split('\n')
        header          = aclust_splt[0].split()
        clust_id        = header[2]
        chr_id          = header[6]
        start           = header[10]
        end             = int(header[12])+1 ##1 added because when opening a range using start and end, end number is not included in range
        value           = (list(range(int(str(start)),int(str(end)))))
        # print ('Cluster:', (value))

        if runType == 'G':                 ## Normal genomic coordinates with integer chr_id
            # print(chr_id)
            # get_chr_id  = str(get_chr_id)
            get_chr_id  = int(get_chr_id)
            chr_id      = int(chr_id)
        else:
            ## Chromosomes are transcript names  i.e. strings
            pass

        
        if get_chr_id == chr_id:
            sm          =   difflib.SequenceMatcher(None,get_value,value) ## The rratio corresponds to larger loci i.e. match/length of longer nucleotides
            # print("Get Value:",get_value)
            # print("Current value",value)
            aratio      = round(sm.ratio(),2)
            # print("Ratio:%s" % (aratio))
            
            if round(sm.ratio(),2) >= matchThres:
                ### Matched - phasiRNA from this cluster
                # print ('\nMatching cluster found:%s' % ''.join(header))
                # print("Allowed Ratio:%s | Current Ratio:%s" % (matchThres,aratio))
                matchCount  +=1
                finalClustL.append((aclust_splt))
                
                phasiCyc    = 0    ## Stores phasing cycles
                phasiSig    = 0    ## Stores total abundance of phase size sRNAs
                otherSig    = 0    ## Stores total abundance of other size sRNAs
                kvalsL      = []   ## List to store kvalues for each phasiRNAs
                
                for i in aclust_splt[1:-1]:## Because header was the first entry of block and not required here, Last entry is always empty
                    # print ("Matched Cluster:\n",i)
                    phasient    = i.split('\t')
                    phasiname   = phasient[4].replace("|","_")
                    phasiseq    = phasient[5]
                    phasilen    = int(phasient[6])
                    phasiabun   = int(phasient[7])
                    phasihits   = int(phasient[10].split("=")[1])
                    phasipos    = int(phasient[3])
                    phasistrand = phasient[2].translate(str.maketrans("+-","wc"))
                    phasipval   = phasient[12]
                    phasikval   = int(phasient[9].split('=')[1])
                    # print(phasipos)
                    # sys.exit()

                    # print(phasiname,phasiabun,phasiseq,phasilen)
                    kvalsL.append(phasikval)
                    tempMatchList.append((phasiname,phasiabun,phasiseq,phasilen,phasihits,phasipos,phasistrand,phasipval))
                    tempAllList.append((phasiname,phasiabun,phasiseq,phasilen,phasihits,phasipos,phasistrand,phasipval)) ## Records all phasiRNAs from all clusters

                    if int(phasilen) == phase:
                        phasiCyc +=1
                        phasiSig += phasiabun
                    else:
                        otherSig += phasiabun
                sizeRatio   = round(phasiSig/(phasiSig+otherSig),2) 
                bestkval    = max(kvalsL) ## Best k-value achieved by this cluster
                # print("Current Cycles:%s | Current sig. strength:%s" % (phasiCyc,phasiSig))
                # print("Current Cycles:%s | Current sig. strength:%s" % (bestkval,phasiSig))
                tempMatchList.append((bestkval,phasiSig,phasID,clust_id,sizeRatio))

                ## Decide the best and remove other from list ##############################
                ############################################################################
                if finalMatchList:
                    ## There exists a previosly matched cluster
                    exist_bestkval = finalMatchList[-1][0]
                    exist_phasiSig = finalMatchList[-1][1]
                    # print("Existing Cycles:%s | Existing sig. strength:%s" % (exist_bestkval,exist_phasiSig))

                    if bestkval > exist_bestkval: ## New cluster has more cycles
                        del finalMatchList[0:]
                        finalMatchList = list(tempMatchList)
                        # print("--- New cluster selected ---")

                    elif bestkval == exist_bestkval: ## Both have same cycles
                        if phasiSig > exist_phasiSig: ## New one has more total abundance of phased siRNAs
                            del finalMatchList[0:]
                            finalMatchList = list(tempMatchList)
                            # print("--- New cluster selected ---")
                    
                    else: ## Existing/old one was long i.e. had more cycles
                        # print("Earlier recorded cluster is retained")
                        pass

                else: ## This is the first cluster
                    finalMatchList  = list(tempMatchList)
                    allphasiList    = list(tempMatchList) 

                # print("\nFinal Match List:",finalMatchList)
        else:
            # print("No Match with this cluster")
            pass

    tempAllList.append((bestkval,phasiSig,phasID,clust_id,sizeRatio)) ## This list has phasiRNAs from all clusters but to keep the structure same as original resList, helpful while writing results, this info is added

    phasinfo = [aname,apval,get_chr_id,get_start,get_end,alib]
    # resList2.append((phasID,tempAllList,phasinfo))
    # resList.append((phasID,finalMatchList,phasinfo)) ## Add best matched cluster entry to the final list, there has to be one best matched entry per PHAS

    resList     = (phasID,tempAllList,phasinfo)
    resList2    = (phasID,finalMatchList,phasinfo)
    
    allMatchCount += matchCount ## Add the matched cluster for each entry
    if matchCount > 0 :
        uniqMatchCount+=1
        
    # print("\nTotal phas loci: %s | Matched: %s" % (len(phasList),len(resList)))
    # print ("SUMMARY: Phased loci in input file:%s | Loci match threshold: %s |Uniq matched cluster: %s | Total matched clusters found:%s" % (phasCount,matchThres,uniqMatchCount,allMatchCount))
    # print("NOTE: If matched clusters more then phased loci that means same cluster was present in different libs\n")
    # print("NOTE: Don't forget to uniq the miRNAs")

    entL = [resList,resList2,finalClustL]

    return entL

def allphasiWriter(clustfile,resList):
    '''
    This takes all the phasiRNAs on a transcript and return unique phasiRNAs file from all clusters - You need those file for study of precursors. phasiRNAs need to unique'd because in cases of foldback same phasiRNAs could match to 5'-w and 3'-c
    '''

    print("\n#### Function:allphasiWriter ##########")

    outfile         = "%s_thres%s_allphasi.csv" % (clustfile.rpartition(".")[0],matchThres)
    # outfile2        = outfile = clustfile+'thres_%s_allphasi.csv' % matchThres                       ## 2437.txt.score_p1e-07_sRNA_21_out.cluster
    fh_out          = open(outfile,'w')
    # print(outfile)
    # sys.exit()

    ## Find entries with unique seq
    for ent in resList: ## entry corresponds to one phased loci
        # print("\nEntry",ent)
        phasID      = ent[0]
        phasCycles  = ent[1][-1][0]
        phasSig     = ent[1][-1][1]
        clustID     = ent[1][-1][2]
        sizeRatio   = ent[1][-1][4]
        phasiList   = ent[1][0:-1]
        phasiList_s = sorted(phasiList,key=itemgetter(5)) ## Sort on position so that results from multiple clusters are at-least sorted in postions
        # print("Sorted all phasiRNAs:", phasiList_s)
        # sys.exit()
        # print("%s | phasCycles:%s | phasSig:%s" % (phasID,phasCycles,phasSig))

        tagset = set() ## To record only the uniq tags, in case of lilium many a times same tag is documented twice
        for i in phasiList:
            print("-Phasi",i)  ## Phasiname,phasiabun,phasiseq,phasilen,phasihits,phasipos
            tag = i[2]
            
            if tag not in tagset: ## Avoid recordinf same tag agains as observed in casne on lilium
                tagset.add(tag)
                
                ## Write phasiRNAs 
                if phasiLenFilter == 'Y': ## If tags filter is ON
                    if (int(i[3]) == int(phase)) and (int(i[1]) >= minAbun): ### Size specified in settings
                        # print(phasID,clustID,i[0],i[1],i[2])
                        fh_out.write('>%s_Clust%s_%s,%s,%s,%s,%s,%s,%s\n%s,%s,%s,%s,%s,%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],i[5],i[3],i[6],i[7],tag,i[1],i[4],i[5],i[3],i[6],i[7]) )  ## phasID,clust_id,phasiname,phasiabun,phasihits,phasiseq,phasiabun,phasihits,phasipos,phasistrand,phasipval
                        tagset.add(tag)

                else:
                    if int(i[1]) > minAbun:
                        # print(phasID,clustID,i[0],i[1],i[2])
                        fh_out.write('>%s_Clust%s_%s,%s,%s,%s,%s,%s,%s\n%s%s,%s,%s,%s,%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],i[5],i[3],i[6],i[7],tag,i[1],i[4],i[5],i[3],i[6],i[7]) )  ## phasID,clust_id,phasiname,phasiabun,phasihits,phasiseq,phasiabun,phasihits,phasipos,phasistrand,phasipval
                        tagset.add(tag)
                        pass


    fh_out.close()
    print("\n\nExiting :allphasiWriter\n")


    return outfile

def writer_summ(clustfile,resList,dictList,pcutoff):
    '''
    write the results
    '''
    print("\n#### Fn: Writer_summ #########################")

    outfile = "./%s_thres%s_phasi.csv" % (clustfile.rpartition(".")[0],matchThres)
    fh_out  = open(outfile,'w')

    if fetchMax == 1: ### Fetch max phasi for each loci 
        outfile2    = "./%s/%sPHAS_p%s_summary.txt" % (res_folder,phase,pcutoff)
        fh_out2     = open(outfile2,'w')

        libsHead    = libs
        fh_out2.write("Name\tP-val\tChr\tStart\tEnd\tIdentifier\tBest k-val\tPhasi ratio\t Max Tag Ratio\t%s\tTotal Phasi Abundance\tMost Abun Tag (MAT)\t MAT Abun\tMAT2\tMAT2 Abun\tBestLib\n" % ('\t'.join(x for x in libsHead)))
        # print(queryLibs) ## Libs whose abindance will be summed to give final abundance of tags
        abunList    = [] ### List to capture tags and abundance for each phased loci
        libAbunList = [] ## Lib-wise abundances of tag

    print("Writing:")
    print("--PhasiRNAs extended CSV file")
    print("--Summary file with lib-specific abundances")
    for ent in resList: ## entry corresponds to one phased loci
        phasID      = ent[0]
        phasCycles  = ent[1][-1][0]
        phasSig     = ent[1][-1][1]
        clustID     = ent[1][-1][2]
        sizeRatio   = ent[1][-1][4]
        phasiList   = ent[1][0:-1]
        phasinfo    = ent[2] ### The PHAs coordinates from phaser
        phastart    = int(phasinfo[3])
        phasend     = int(phasinfo[4])
        phaslen     = phasend-phastart
        # print("phasinfo:",phasinfo)
        # sys.exit()
        # print("%s | phasCycles:%s | phasSig:%s" % (phasID,phasCycles,phasSig))


        tagset = set() ## To record only uniq tags, in case of lilium many a times same tag is documented twice
        for i in phasiList:
            # print("-Phasi",i)  ## Phasiname,phasiabun,phasiseq,phasilen,phasihits,phasipos
            tag = i[2]
            
            if tag not in tagset: ## Avoid recording same tag agains as observed in casne on lilium
                tagset.add(tag)
                
                ## Write phasiRNAs 
                if phasiLenFilter == 'Y': ## If tags filter is ON
                    if (int(i[3]) == int(phase)) and (int(i[1]) >= minAbun): ### Size specified in settings
                        # print(phasID,clustID,i[0],i[1],i[2])
                        fh_out.write('>%s_Clust%s_%s,%s,%s,%s,%s,%s,%s\n%s,%s,%s,%s,%s,%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],i[5],i[3],i[6],i[7],tag,i[1],i[4],i[5],i[3],i[6],i[7]) )  ## phasID,clust_id,phasiname,phasiabun,phasihits,phasiseq,phasiabun,phasihits,phasipos,phasistrand,phasipval
                        tagset.add(tag)

                else:
                    if int(i[1]) > minAbun:
                        # print(phasID,clustID,i[0],i[1],i[2])
                        fh_out.write('>%s_Clust%s_%s,%s,%s,%s,%s,%s,%s\n%s%s,%s,%s,%s,%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],i[5],i[3],i[6],i[7],tag,i[1],i[4],i[5],i[3],i[6],i[7]) )  ## phasID,clust_id,phasiname,phasiabun,phasihits,phasiseq,phasiabun,phasihits,phasipos,phasistrand,phasipval
                        tagset.add(tag)
                        pass


                ## Get max abundance phasiRNA for specified phase
                if fetchMax == 1:

                    ## Get lib-wise abundaces mode
                    if fetchLibAbun == 0:
                        if len(tag) == int(phase):
                            atag,abun_sum,lib_abun = getAbundanceLocal(tag,dictList)
                            libAbunList.append((lib_abun))

                    elif fetchLibAbun == 1: ## All the tags
                        atag,abun_sum,lib_abun = getAbundanceLocal(tag,dictList)
                        libAbunList.append((lib_abun))
                    
                    else:
                        print("Libwise abundances won't be fetched")
                        pass

                    ## Tag specific abundances for fetching most abundant tag
                    if len(tag) == int(phase): ### Size specified in settings
                            abunList.append((atag,abun_sum))
            else:
                # print("Tag recorded once already#####################################\n")
                # sys.exit()
                pass


        # print("Elements in phasiList:%s" % (len(phasiList)))
        # print("Elements in libAbunList:%s" % (len(libAbunList)))


        ## Process Fetch max results for this phased loci
        #################
        if fetchMax == 1:
            
            abunList_sort = sorted(abunList, key=operator.itemgetter(1),reverse=True) ## Sort list on abundances to get max abundant phasiRNA
            # print("\nExample sorted values:%s" % (abunList_sort[0:10]))
            maxTag      = abunList_sort[0][0]   ## Max abundant phasiRNA sequence
            maxAbun     = abunList_sort[0][1]   ## Abundance of max abundant phasiRNA
            if len(abunList_sort) > 1: ## In one case no second tag was found
                maxTag2     = abunList_sort[1][0]   ## Max abundant phasiRNA sequence
                maxAbun2    = abunList_sort[1][1]   ## Abundance of max abundant phasiRNA
            else:
                maxTag2     = "na"
                maxAbun2    = "0"
                print("Prediction is of really low quality - Just one tag of phase size found")
                time.sleep(2)
                # sys.exit()
            totalAbun   = sum(int(i[1]) for i in abunList_sort) ## Phased abundance
            nminabun    = sum(int(i[1]) >= minAbunCut for i in abunList_sort)
            maxTagRatio = round(maxAbun/totalAbun,2)


            ## Sum Lib-wise abundances for all tags
            libAbunSum = [0]*len(libAbunList[0]) ## This will hold sum of all tags, intialized for number of libraries
            for tag in libAbunList:
                # print(tag)
                libAbunSum  = [sum(x) for x in zip(libAbunSum,tag)]
            
            # print("Tag:%s | maxPhasTag:%s | totalPhasAbun:%s" % (maxTag,maxAbun,totalAbun))
            # print("Libwise Abundances",libAbunSum)

            ## Write - Loci, most abundant tag, most abundant tag abun,total phased abun, number of phased tags, and lib-wise abundances
            if args.safesearch == "T":
                ## Apply filters
                # print("Total Abudnace:",totalAbun)
                if (float(maxTagRatio) <= maxTagRatioCut) and (int(phasCycles) >= phasCyclesCut) and (phaslen >= (phase*6)+3) and (nminabun >= 1): ##  3-nt addded to PHAS len to account for dicer offsets
                    fh_out2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('\t'.join(str(x) for x in phasinfo[:-1]),phasID,str(phasCycles),sizeRatio,maxTagRatio,'\t'.join(str(x) for x in libAbunSum),totalAbun,maxTag,maxAbun,maxTag2,maxAbun2,phasinfo[-1]))
                else:
                    ## Bad Quality prediction - Skip
                    pass

            else:
                ### No filters applied
                fh_out2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('\t'.join(str(x) for x in phasinfo[:-1]),phasID,str(phasCycles),sizeRatio,maxTagRatio,'\t'.join(str(x) for x in libAbunSum),totalAbun,maxTag,maxAbun,maxTag2,maxAbun2,phasinfo[-1]))
            
            abunList    = [] ## Empty before next phased loci
            libAbunList = [] ## Empty lib-iwse abundances of tag list before next entry


    fh_out.close()
    fh_out2.close()

    return outfile,outfile2

def clustWriter(finalClustL,pcutoff):
    '''
    writes the final cluster file
    '''

    outfile = "./%s/%sPHAS_p%s_clust.txt" % (res_folder,phase,pcutoff)
    fh_out  = open(outfile,'w')

    acount = 0
    for aclust in finalClustL:
        # print(aclust)
        fh_out.write("%s\n" % ("\n".join(aclust)))
        acount+=1

    print("Clusters written:%s" % (acount))

    fh_out.close()

    return outfile

def getAbundanceLocal(tag,dictList):
    """The core function for getting the abunance from a single tag
    from a series of libraries provided in tag count format. Makes
    the calls to create the data structures to get the abundance of
    a tag across several libraries

    Args:
        tag: Sequence to search the libraries for
        tagCountFilenamesList: The names of all libraries
    Returns:
        The input tag, abundance in all libraries, and the sum of 
        abundance across all libraries. 
    
    """
    # Initialize a list for the abundance results
    tagAbunList = []
    # Search for the tag in all libraries
    for dict in dictList:
        tagAbunList.append(getAbunFromDict(dict, tag))

    # The dictionaries may not be sorted in the given order, so sort them
    # according to the key (that being the index returned by getAbunFromDict)
    # print("\nTag:",tag)
    # print("Unsorted:",tagAbunList)
    sortedTagAbunList = sorted(tagAbunList, key=operator.itemgetter(0),
        reverse = False)
    # print("Sorted:",sortedTagAbunList)

    # Store just the abundances for return
    abunList = []
    for abun in sortedTagAbunList:
        abunList.append(abun[1])

    return(tag,sum(abunList),abunList)

def getAbundance(cur,tag,finalLibs):
    '''Input is tag for each loci and out put is tag with maximmum abumdance and sum of phasiRNAs - 
    rewritten in v1.0 for fetching tags from run master'''

    lib_abun = [] ## list to hold lib-wise abudnances
    
    for alib in finalLibs:
        # print("Lib:",alib)
    
        cur.execute("SELECT tag,norm FROM %s.run_master where tag = '%s' and lib_id = %s" % (db,tag,alib))### Convert intergenic to gene name so as to get strand
        info = cur.fetchall() ## Entries are redundant, one for every hit
        # print("Query fetched", info)

        if info:
            atag,norm_abun = info[0]
            lib_abun.append(norm_abun)
            # print("--Tag abundance:%s for lib:%s"% (tag,norm_abun))
        else:
            norm_abun = 0
            lib_abun.append(norm_abun)
            # print("--Tag abundance:%s for lib:%s"% (tag,norm_abun))


    abun_sum = sum(lib_abun)
    # print("--Lib-wise abundances",lib_abun)
    # print("--Sum of abundances:%s\n" % (abun_sum))

    # sys.exit()

    return tag,abun_sum,lib_abun
            
def prepareQuery(excludeLibs,cur):

    ### Prepare query of libs #################

    ### Get lib names
    columns = [] ## empty list
    cur.execute("SELECT DISTINCT(lib_id) FROM %s.run_master" % (db))
    info = cur.fetchall()
    libs = [x[0] for x in info]

    print("\nLibs:",libs)

    if excludeLibs:
        print("\n\nLibs specifed in excludeLibs %s will be skipped\n\n" % (excludeLibs))
        selectLibs      = [] ## Columns excluding the unwanted libraries
        excludeLibs_s   = [str(i) for i in excludeLibs] ## Converting all entries in exclude list to string for matching below
        
        for i in libs:

            ## Check if user made mistake in libType - as that will give results for all entries
            if type(i) is int and libType == 1:
                print("You seem to have input lib_id and chosen wrong libType")
                print("Check libType and excludeLibs match - Script will exit now")
                sys.exit()
            elif type(i) is str and libType == 0:
                print("You seem to have input lib_id and chosen wrong libType")
                print("Check libType and excludeLibs match - Script will exit now")
                sys.exit()
            else:
                print("All seems well")
                pass

            ### Filter libraries
            if str(i) not in excludeLibs_s: ## Tested OK
                selectLibs.append(i)
            else:
                print("excluded:",i)
                # sys.exit()
                pass
        
        finalLibs = selectLibs ## Norm_Sum and max_norm are not included

    else:
        finalLibs = libs ## ## Norm_Sum and max_norm are not included

    # print("finalLibs:%s" % (finalLibs))
    
    lib_col =",".join(str(x) for x in finalLibs)### Manually mentioned strating lib column - should work on all tag position summary tables
    
    print("\nLibrary Columns:",lib_col)
    # queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(")
    sumLibs = "%s" % lib_col.replace(",","+")
    queryLibs = "%s" % lib_col.replace(",",",")
    print("\nThese are sumLibs:",sumLibs)
    print("\nThis is query Libs:",queryLibs)
    # sys.exit()

    return queryLibs,sumLibs,finalLibs           

def FileCombine(alist,aname):

    # print ('--Files to concatanate:',alist)
    fh_out = open(aname ,'w')
    
    for x in alist:
        # print ("--Concatanating:%s"% (x))
        afile   = open('./%s' % (x), 'r')
        data    = afile.read()
        afile.close()
        fh_out.write(data)
    
    fh_out.close()
        
    return aname

def cleaner():
    '''
    Picks and cleans unwanted files
    '''
    print ("\n#### Garbage Cleaner #########################")
    garbage = [afile for afile in os.listdir('%s/' % (res_folder)) if afile.endswith (('.cluster','.zip'))] 
    print("Garbage:",",".join(garbage))

    for afile in garbage:
        apath =   "%s/%s" % (res_folder,afile)  
        if os.path.isfile(apath): ## Check to see its a file from bowtie and not tophat mapped folder - Untested
            print("Deleting:%s" % (apath))
            os.remove(apath)
        else:
            # print("Skipping cleanup, as its a directory %s" % (afile))
            pass

    return None

def readFileToDict(index):
    """This function simply reads all tags in a tag count file
       and stores them in a dictionary

    Args:
        filename: The name of an individual file to read
        index: A numberical index for the order that the results should
               eventually be provided
    Returns:
        Dictionary of the read file

    """

    filename = libs[index]
    print("Caching %s to fetch phasiRNA abundances" % (filename))

    # Initialize an empty dictionary for the results
    # Use the index as a key to the dictionary so that the file
    # can be associated with this dictionary later
    fileDict = {index: {}}

    if libFormat == "T":
    # Read through the file and store the tag count file into 
    # the dictionary
        fh_in       = open(filename, 'r')
        for line in fh_in:
            # tag and abundance should be separated by a tab. Strip
            # the new line at the end of each line
            tag, abun = line.strip('\n').split('\t')

            # Remove any white space that is inadvertently left
            tag = tag.strip()
            abun = abun.strip()

            # Store the 
            fileDict[index][tag] = int(abun)

        fh_in.close()

    elif libFormat == "F":
        
        ### Read files
        filename    = "%s.fas" % filename.rpartition(".")[0] 
        fh_in       = open(filename, 'r')
        fasta       = fh_in.read()
        fasta_splt  = fasta.split('>')
        
        for i in fasta_splt[1:]:
            ent     = i.split('\n')
            abun    = ent[0].split("|")[1].strip()
            tag     = ent[1].strip()
            # print("Ent:",ent)
            # print("Tag:%s | Abun:%s" % (tag,abun))

            fileDict[index][tag] = int(abun)

        fh_in.close()

    else:
        ## @libFormat not accepted
        sys.exit()

    return(fileDict)

def getAbunFromDict(dict, tag):
    """Get the abundance of the tag from a dictionary

    Args:
        tag: The tag to be queried from the dictionary
        dict: The dictionary of tags and their abundances in ech library
        libList: The list of libraries as provided by the user
    Returns:
        Three separate elements. The first is the tag, the second
        is the summed abundance across all libraries, and the final
        is a list of the individual library abundances

    """

    # Pull the library index from the dictionary
    index = list(dict.keys())[0]

    # Try to abundance of tag. If it doesn't exist in this library,
    # return 0 for the abundance
    try:
        return([index, dict[index][tag]])

    # Tag doesn't exist in this dictionary so return 0
    except KeyError:
        return([index,0])

def PPResults(module,alist):
    npool   = Pool(int(nproc))
    res     = npool.map_async(module, alist)
    results = (res.get())
    npool.close()

    return results

def percentile(data, percentile):
    size = len(data)
    return sorted(data)[int(math.ceil((size * percentile) / 100)) - 1]

#############################
#### ANNOTATE ###############

def gtfParser(afeatureFile):

    '''Parses PASA gtf file into featurename,chr,start,end,strand,feature'''

    print("\nFunction: gtfParser")
    
    with open(afeatureFile) as fh_in:
        lines = (line.rstrip() for line in fh_in) 
        gtfRead = list(line for line in lines if line) # Non-blank lines in a list
    fh_in.close()

    gtfList = [] ## List to hold parsed GTF entries

    for i in gtfRead:
        # print(i)
        ent = i.split("\t")
        # print("\nEnt",ent)
        if ent[2] == "transcript" or ent[2] == "exon":
            # print(ent)
            gchr    = ent[0].replace("chr","")
            gtype   = ent[2]
            gstart  = ent[3]
            gend    = ent[4]
            gstrand = ent[6].translate(str.maketrans("+-","wc"))
            aflag   = "P" ## PASA 
            info    = ent[8].strip("\n").split(";")
            # print(info,len(info))
            if len(info) == 3: ## With last one empty 
                ## Protein coding gene with a version number
                gid     = info[0].split()[1].replace('"','') ## Gene ID
                tid     = info[1].split()[1].replace('"','') ## Transcript ID
                # print(gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))
            # if mode == 6 and len(info) == 4: ## With last one empty 
            #     ## Protein coding gene with a version number
            #     gid     = info[0].split()[1].replace('"','') ## Gene ID
            #     tid     = info[1].split()[1].replace('"','') ## Transcript ID
            #     # print(gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
            #     gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))
            else:
                print("This entry has more than expected info")
                print(info)
                print("Debug!!")
                sys.exit()

        else:
            print("We don't need this info for current script")
            pass

    if len(gtfList) == 0:
        print("Check if the feature used for extraction i.e. gene, transcript, exon is correct")
        print("Debug!!")
        sys.exit()

    print("First 10 entries of gtfList list: %s" % (gtfList[:10]))
    print("Total entries fetched from GTF file:%s" % (str(len(gtfList))))

    print ("Exiting function - gtfParser\n")
    time.sleep(1)
    
    return gtfList

def gtfParser2(afeatureFile):

    '''This function parses Trinity and Rocket GTF file
    to give a trascript entry and entry for all exons - Basically this is the parser for gffread geenrated files'''

    print("\n#### Fn: gtfParser ###########################")

    #### Check File
    if os.path.isfile(afeatureFile):
        pass
    else:
        print("---GTF file could not found")
        print("---Please check if it exists")
        sys.exit()
    
    #### Read file
    with open(afeatureFile) as fh_in:
        lines = (line.rstrip() for line in fh_in) 
        gtfRead = list(line for line in lines if line) # Non-blank lines in a list
    fh_in.close()

    #### Parse File
    print("Parsing '%s' gtf file" % (afeatureFile))
    gtfList     = [] ## List to hold parsed GTF entries
    tempName    = [] ## Stores current trans name
    tempCoords  = [] ## Temp coords
    for i in gtfRead:
        # print(i)
        ent = i.split("\t")
        # print("\nEnt",ent)
        gScore  = ent[5] ## Score provided for mapping accuracy in Trinity GTF from GMAP 0 to 100
        gtype   = ent[2]
        if gtype == "exon" and (gScore == '.' or float(gScore) > 90.0):
            # print("\nExon:",ent)
            gchr    = re.sub("[^0-9]", "", ent[0]).lstrip('0')
            # gchr    = ent[0]        ## To include scaffolds and fragments
            gstart  = int(ent[3])
            gend    = int(ent[4])
            gstrand = ent[6].translate(str.maketrans("+-","wc"))
            info    = ent[8].strip("\n").split(";")
            # print(info)
            
            ## Parse the info and add the exon entries  #################
            #############################################################

            if annomode == 2 or (annomode == 4 and len(info)==4): ## This is  trinity GTF info
                ## Protein coding gene with a version number
                tid     = info[0].split()[1].replace('"','').split(".")[0] ## Transcript ID
                gid     = info[2].split()[1].replace('"','').rpartition("_")[0] ## Gene ID
                aflag   = 'T' ## Trinity
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))

            elif annomode == 3 or (annomode ==4 and len(info) >= 7): ## This is rocket GTF info, with and w/o p_id "P1" in the end
                tid     = info[0].split()[1].replace('"','') ## Transcript ID
                gid     = info[1].split()[1].replace('"','') ## Gene ID
                aflag   = 'R' ## Rocket
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))
            
            elif annomode == 5:
                tid     = info[0].split()[1].replace('"','').split(".")[0] ## Transcript ID
                gid     = info[1].split()[1].replace('"','').split(".")[0] ## Gene ID
                aflag   = 'PB' ## PacBio
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))

            elif annomode ==6:
                gid     = info[0].split()[1].replace('"','') ## Gene ID
                tid     = info[1].split()[1].replace('"','') ## Transcript ID
                aflag   = 'P' ## Pasa
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))


            else:
                print("-This entry has more than expected info")
                print("-Info block",info)
                print("-Info block len:%s" % len(info))
                print("-Debug!!")
                sys.exit()

            ## Check trascript change #########
            ###################################

            ## Check if transcript entry needs to be added i.e. all the exons from previous one has been captured
            if not tempName or tempName[0] == tid:
                tempName.append(tid)
                tempCoords.append(gstart)
                tempCoords.append(gend)
                tempStrand  = gstrand   ## strand of last transcript
                tempgid     = gid       ## gene name of last transcript
                tempChr     = gchr      ## chr of last transcript
                # print(tempCoords)
            elif tempName[0] != tid: ## Check if the new transcript is being read. if yes then add transcript entry for old exons using tempName and tempCoords
                # print("-New transcript read - summarizing transcript entry for earlier one")
                tstart      = min(tempCoords)
                tend        = max(tempCoords)
                ttid         = tempName[0]
                ttype       = 'transcript'
                # print('-',tempgid,ttid,tempChr,tstart,tend,tempStrand,ttype,aflag)
                gtfList.append((tempgid,ttid,tempChr,tstart,tend,tempStrand,ttype,aflag))
                # sys.exit()
                
                ## Empty lists and fill with current exon from new transcript
                tempName    = [] ## Empty trans name
                tempName.append(tid)
                tempCoords  = []
                tempCoords.append(gstart)
                tempCoords.append(gend) ## Empty trans coords
                tempStrand  = gstrand   ## strand of last transcript
                tempgid     = gid       ## gene name of last transcript
                tempChr     = gchr      ## chr of last transcript
                # sys.exit()
            else:
                print("-Unforseen scenario encountered")
                sys.exit()

        else:
            # print("We don't need this info for current script") ## CDS or GENE
            pass

    if len(gtfList) == 0:
        print("Check if the feature used for extraction i.e. gene, transcript, exon is correct")
        print("Debug!!")
        sys.exit()
    else:
        # print("First 10 entries of gtfList list: %s" % (gtfList[:10]))
        pass

    print("Total entries fetched from GTF file:%s" % (str(len(gtfList))))
    time.sleep(1)
    
    return gtfList

def overlapChecker(phasList,gtfList,pcutoff):
    '''Checks for overlap between genomic PHAS and transcipts frpm GTF file'''
    
    print("\n#### Fn: Annotate ############################")
    
    ## Prepare output file ####################################
    ###########################################################
    outfile = "./%s/%sPHAS_p%s_annotation.txt" % (res_folder,phase,pcutoff)
    fh_out  = open(outfile,'w')
    fh_out.write("Name\tP-val\tphasChr\tphasStart\tphasEnd\tidentifier\tOverlapping Transcript\toverlapping trans NTs\tPercentage Overlap to Trans\toverlapping Exon NTs\tPercentage overlap to Exon\tOverlapping Exons\tTrascript Strand\tTranscript Length\n")
    
    ## Prepare DB and featureTables ##########################
    ##########################################################
    DB = 'tempdb'
    try:
        os.remove(DB)
    except OSError:
        pass
    ## Prepare features table - This might include entries from two GTF files (mode-2), identified by the flags
    conn            = sqlite3.connect(DB)
    featureTable    = tableMaker(gtfList,conn) ## TURN ON ****
    
    # # Test Query
    # cur = conn.cursor()
    # cur.execute("SELECT * FROM %s where chr = 1 AND strand = 'w' limit 10" % (featureTable))
    # test = cur.fetchall()
    # print("Test results:",test)
    # # sys.exit()

    ## Flags i.e. data to query - This allows combining results from Trinity and rocket to one sheet
    flags       = [] ## Datatype to query
    if annomode     == 1 or annomode == 6:
        flags.append(('P'))
    elif annomode   == 2:
        flags.append(('T'))
    elif annomode   == 3:
        flags.append(('R'))
    elif annomode   == 4:
        flags.append(('R'))
        flags.append(('T'))
    elif annomode   == 5:
        flags.append(('PB'))
    else:
        print("Please check the mode selected for analysis")
        pass

    ## Find transcript overlapping PHAS ########################
    ############################################################
    print("Analysis for overlapping trans initialized")
    
    nonOverlapList  = []    ## Temp list so that these non-overlapping ones are just added once,as there is no uniq transcipt attached to them
    overlapNameList = []    ## Just the naemes of phased entries that has overlapping results
    for aflag in flags:     ## For provided datatypes
        # print("\nFlag:",aflag)
        
        for i in phasList: ## FOr every PHAS
            # print("\n###################Entry:",ent)
            aphasID     = i[0]
            aphasCycles = i[1][-1][0]
            aphasSig    = i[1][-1][1]
            aclustID    = i[1][-1][2]
            asizeRatio  = i[1][-1][4]
            aphasiList  = i[1][0:-1]
            aent        = i[2][:-1]              ## aname, apval, achr, astart, aend, trash
            aname,apval,achr,astart,aend = aent
            transList   = overlapTrans(aent,conn,featureTable,aflag)
            matchflag   = False

            ## Check overlap with exons of every overlapping transcript
            for trans in transList:             ## For every overlapping trasc
                atrans,tstrand,tlen,toverlap = trans
                print("-Computing overlap for transcript: %s" % (atrans))
                exonsOverlap,nexons     = overlapExons(aent,conn,featureTable,atrans,aflag)
                tperc                   = round((toverlap/(aend-astart+1)),3)
                eperc                   = round((exonsOverlap/(aend-astart+1)),3)

                ## Filters and output results ####
                ##################################
                if tperc >= overlapPerc and toverlap >= ntoverlapCutoff:
                    matchflag   = True
                    overlapNameList.append(aphasID)
                    fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('\t'.join(str(x) for x in aent), aphasID, atrans, str(toverlap), str(tperc), str(exonsOverlap), str(eperc), str(nexons), tstrand, str(tlen)))
                else:
                    # nonOverlapList.append((aent,aphasID)) ## This in current form wont catch those entries for which there is no overlap with transcripts i.e. PHAS with no overlap to any gene. The filter is designed to include those cases for which there is an overlap to transcripts but not with exons
                    pass

            ## Record non-overlapping
            if matchflag == False:
                nonOverlapList.append((aent,aphasID))

    #### Write results for which no overlapping transcript was found
    writtenList = [] ## Temp list to check the result for entry has been written to avoid duplicate for PHAS which has no results
    # print(overlapNameList)
    for bent in nonOverlapList:
        # print(bent)
        bphasID     = bent[1]
        if (bphasID) not in overlapNameList and (bphasID not in writtenList):
            fh_out.write("%s\t%s\tx\tx\tx\tx\tx\tx\tx\tx\n" % ('\t'.join(str(x) for x in bent[0]),bphasID,))
            writtenList.append(bphasID)
     
    fh_out.close()

    return outfile

def tableMaker(alist,conn):
    '''makes SQLlite table for query'''

    print("Preparing SQL table with GTF features")

    cur = conn.cursor()

    ## Make table
    featuretable = "tempTable"
    cur.execute('''DROP TABLE IF EXISTS %s''' % (featuretable)) ### Drop Old table - while testing    
    conn.commit()
    
    try:
        cur.execute('''CREATE TABLE %s (gene varchar(255),trans varchar(255),chr varchar(255), start integer, end integer, strand varchar(10), type varchar(255),flag varchar(255))''' % (featuretable))
        conn.commit()
        acount = 0
        for i in alist:
            # print(i)
            gid,tid,gchr,gstart,gend,gstrand,gtype,aflag = i
            if gchr:
                ## Check if there is a chr assigned otherwise comaprsion cant be done
                acount+=1
                # print ('Gene:',gid,'Trans:',tid,'| Chr:',gchr,'| Start:',gstart,'| End:',gend, '| Strand:',gstrand,'| Type:',gtype,'| Flag:',aflag)
                cur.execute("INSERT INTO %s VALUES ('%s','%s','%s',%d,%d,'%s','%s','%s')" % (featuretable,str(gid),str(tid),str(gchr),int(gstart),int(gend),str(gstrand),str(gtype),str(aflag))) ### Integer made string to accomodate scaffolds or fragments
                # featureList.append((str(aname),int(achr),int(astart),int(aend),str(astrand)))
                #conn.commit()
    
    except sqlite3.Error:
        print('ERROR:',Error)
        sys.exit()

    print("Feature table made with %s entries" % (acount))

    return featuretable

def overlapTrans(ent,conn,featureTable,aflag):
    '''This function returns a list of transcipts that overlaps with a PHAS'''
    cur         = conn.cursor()
    transList   = [] ## Store the transcipts that overlap
    # print("++++",ent)
    aname,apval,achr,astart,aend = ent

    # print("-Flag being queried:%s" % (aflag))
    # print("-Query entry:",ent)
    
    ## trans flanking phas or enclaved in phas
    cur.execute("SELECT * FROM %s where chr = %s AND flag = '%s' AND ((end between %s and %s) or (start between %s and %s)) AND type = 'transcript'" % (featureTable,achr,aflag,astart,aend,astart,aend))
    flankTrans = cur.fetchall()

    ## PHAS enclaved in trans
    cur.execute("SELECT * FROM %s where chr = %s AND flag = '%s' AND (%s between start and end) AND (%s between start and end) AND type = 'transcript'" % (featureTable,achr,aflag,astart,aend))
    bigTrans = cur.fetchall()


    # print("flanking Trans:",flankTrans)
    # print("Long Trans:",bigTrans)

    ## Combine both lists and report overlapping trans
    allTrans    = flankTrans + bigTrans
    for i in allTrans:
        # print(i)
        atrans  = i[1]
        astrand = i[5] 
        alen    = i[4]-i[3]
        cur.execute("SELECT * FROM %s where flag = '%s' AND trans = '%s' AND type = 'transcript'" % (featureTable,aflag,atrans))
        transinfo = cur.fetchall() 
        # print(transinfo)

        toverlap = 0
        tstart  = transinfo[0][3]
        tend    = transinfo[0][4]
        
        ## Trans is enclaved
        if astart <= tstart and aend >= tend:
            # print("-exon enclaved")
            toverlap = tend - tstart + 1

        ## PHAS is enclaved
        elif tstart <= astart and tend >= aend:
            # print("-PHAS enclaved")
            toverlap = aend - astart + 1
        
        ## Overlap at 5' of PHAS, use phas start as reference
        elif tstart <= astart and tend >= astart:
            # print("- 5' flank")
            toverlap = tend - astart + 1
        
        ## Overlap at 3' of PHAS, use PHAS end as reference
        elif tstart <= aend and tend >= aend:
            # print("- 3' flank")
            toverlap = aend - tstart + 1

        else:
            # print("-Unexpected overlap encountered - Transcript Overlaps but not Exons?")
            pass

        transList.append((atrans,astrand,alen,toverlap))

    # print("-Number of overlapping trans:%s" % (len(transList)))
    # print("-transList",transList)

    # if len(transList) > 0:
    #     sys.exit()

    return transList

def overlapExons(ent,conn,featureTable,atrans,aflag):
    '''This function will compute overlaps wit exons of overlapping transcripts'''

    aname,apval,achr,astart,aend = ent
    
    cur = conn.cursor()
    cur.execute("SELECT * FROM %s where flag = '%s' AND trans = '%s' AND type = 'exon'" % (featureTable,aflag,atrans))
    exons = cur.fetchall()
    # print("-These are the exons",exons)

    ## compute overlap with exons of this transcript
    nexons          = len(exons)
    exonsOverlap    = 0
    for aexon in exons:
        # print("-Checking exons:",aexon)
        xoverlap    = 0
        xstart      = aexon[3]
        xend        = aexon[4]

        ## Exon is enclaved
        if astart <= xstart and aend >= xend:
            # print("-exon enclaved")
            xoverlap = xend - xstart + 1

        ## PHAS is enclaved
        elif xstart <= astart and xend >= aend:
            # print("-PHAS enclaved")
            xoverlap = aend - astart + 1
        
        ## Overlap at 5' of PHAS, use phas start as reference
        elif xstart <= astart and xend >= astart:
            # print("- 5' flank")
            xoverlap = xend - astart + 1
        
        ## Overlap at 3' of PHAS, use PHAS end as reference
        elif xstart <= aend and xend >= aend:
            # print("- 3' flank")
            xoverlap = aend - xstart + 1

        else:
            # print("-Unexpected overlap encountered - Transcript Overlaps but not Exons?")
            pass
            # sys.exit()

        exonsOverlap+=xoverlap

    # print("-Overlap found:%s" % exonsOverlap)
    return exonsOverlap,nexons

#### COMPARE ################
#############################

def summparse(adir):
    '''
    parses phaster file
    '''
    ### Check Folder
    if os.path.isdir(adir): ## Check to see its a file from bowtie and not tophat mapped folder - Untested
        # print("Directory from phaser analysis found")
        pass
    else:
        print("Specified directory not found: %s" % (adir))
        print("Please confirm the directory path, or")
        print("confirm 'phasdetect' run finished successfully")
        print("To see all required parameters, run: python3 phasdetect -h\n")
        print("Script exiting...\n")
        sys.exit()

    summfileL   = [file for file in os.listdir("%s" % (adir)) if file.endswith ('_summary.txt')]


    ### Check File
    if len(summfileL) == 1:
        summfile    = "%s/%s" % (adir,summfileL[0])
    elif len(summfileL) == 0:
        print("No summary file found in %s directory" % (adir))
        print("Please make sure that the orginal file exists")
        print("Did you removed it from folder? Copy it back in original format and rerun script")
        sys.exit()
    elif len(summfileL) > 1:
        print("Multiple summary files found in %s directory" % (adir))
        print("Please make sure that the orginal file exists")
        print("Also remove any files that you made during your analysis from directory and rerun script")
        sys.exit()

    if os.path.isfile(summfile):
        pass
    else:
        print("---Summary file %s not found in %s/ directory" % (summfile.rpartition("/")[-1],adir))
        print("---Please check if file exists and rerun the script")
        sys.exit()

    print("\nFn: Read Summary #############################")
    summL   = []
    summD   = {}

    fh_in   = open(summfile,'r')
    fh_in.readline() ## Remove Header
    afile   = fh_in.readlines()
    fh_in.close()

    acount = 0
    for i in afile:
        ent     = i.strip('\n').split('\t') ## Name, p-val, chr, start, end, identifier, best k-val, phasi ratio, max tag ratio
        # print(ent)
        aname   = ent[0]
        apval   = ent[1]
        achr    = ent[2]
        astart  = ent[3]
        aend    = ent[4]
        aid     = ent[5]
        akval   = ent[6]
        aratio  = ent[7]
        amaxtag = ent[8]
        avalue  = (aname,achr,astart,aend,apval,aid,akval,aratio,amaxtag)
        # print("Phaster parsed:",avalue)
        summL.append(avalue)
        summD[aname] = avalue
        acount       +=1

    print("PHAS in file:%s| PHAS cached:%s" % (acount,len(summL)))

    return summL,summD,summfile

def compare(summL1,summD1,summL2,summD2,log_file):
    '''
    Compares the PHAS loci between two sumamries
    '''

    resList     = []    ## List to store final results
    negSet      = set() ## Nt matching from summary-2
    emptyfill   = ('x','x','x','x','x','x','x','x','x') ### Lazy wasy of filling empty cells
    acount      = 0 ## Count of matched from List1
    bcount      = 0 ## Count of matched from List1
    ccount      = 0 ## Count of unmatched from List1
    dcount      = 0 ## Count of unmatched from List2
    ecount      = 0 ## All the matches between List1 and List2 (redundant matches)
    #### Find Matched PHAS
    for aphas in summL1:
        # print(aphas)
        aname   = aphas[0]
        achrid  = aphas[1]
        astart  = int(aphas[2])
        aend    = int(aphas[3])
        aid     = aphas[5]
        akey    = 'a-%s-%s-%s' % (achrid,str(astart),str(aend)) ## 'a' added to diffrentiate the key from PHAS in different file that have same coordinates.
        aregion = list(range(astart,aend))
        matflag = False     ## Tracks the match stats for this aphas

        for bphas in summL2:
            bname   = bphas[0]
            bchrid  = bphas[1]
            bstart  = int(bphas[2])
            bend    = int(bphas[3])
            bid     = bphas[5]
            bkey    = 'b-%s-%s-%s' % (bchrid,str(bstart),str(bend)) ## 'a' added to diffrentiate the key from PHAS in different file that have same coordinates.
            bregion = list(range(bstart,bend))

            if achrid == bchrid:
                sm1         = difflib.SequenceMatcher(None,aregion,bregion) ## Mapping bphas over aphas
                matchratio1 = round(sm1.ratio(),5)

                if (matchratio1 >= 0.25):
                    # print("aphas",aphas)
                    # print("bphas",bphas)
                    # print("Match Ratio:%s" % (matchratio1))
                    matblocks   = sm1.get_matching_blocks() ### [(0,0,200),(),...] A list of all matches with tuple
                                                            ### corresponding to one path cof match. In tuple first
                                                            ### two elements are coordinates for seqeunce1 and sequence2 
                                                            ### followed my length of match. Expect just one tuple as there
                                                            ### is just one patch of match.
                    # print(matblocks)
                    amatcoord   = astart+int(matblocks[0][0])
                    bmatcoord   = bstart+int(matblocks[0][1])
                    matcoords   = "%s:%s" % (amatcoord,bmatcoord)
                    # print("astart:%s | bstart:%s | amatcoord:%s | bmatcoord:%s" % (astart,bstart,amatcoord,bmatcoord))
                    matlen      = int(matblocks[0][2])
                    abphas      = aphas+bphas
                    resList.append((abphas,matchratio1,matcoords,matlen))
                    negSet.add(bid)
                    
                    matflag = True
                    ecount  +=1

                
                else:
                    ## Different PHAS loci on same chromsome and transcript
                    pass
            
            else:
                ## Different chr or trascripts
                pass

        #### aphas unmatched?
        #####################
        if matflag == False:
            ## No match found for this aphas entry
            print("No Match",aphas)
            abphas      = aphas+emptyfill
            matchratio  = 0 
            matcoords   = "none"
            matlen      = "none"
            resList.append((abphas,matchratio,matcoords,matlen))
            ccount      +=1
        else:
            ## This aphas is matched well
            acount      +=1


    #### Include unmatched bphas
    ############################
    for bphas in summL2:
        bid     = bphas[5]
        ccount  = 0
        if bid not in negSet:
            ## This bphas had no match with aphas
            abphas      = emptyfill+bphas
            matchratio  = 0
            matcoords   = "none"
            matlen      = "none"
            resList.append((abphas,matchratio,matcoords,matlen))
            dcount      +=1
        else:
            ## This bphas was matched to at-least one aphas
            bcount      +=1


    print("PHAS in summary1: %s | Matched:%s | Unmatched:%s" % (len(summL1),acount,ccount))
    print("PHAS in summary2: %s | Matched:%s | Unmatched:%s" % (len(summL2),bcount,dcount))
    print("Matched %s:%s (summary1:summary2)\n" % (acount,ecount))

    with open(log_file, "a") as fh_log:
        fh_log.write("\nPHAS in summary1: %s | Matched:%s | Unmatched:%s\n" % (len(summL1),acount,ccount))
        fh_log.write("PHAS in summary2: %s | Matched:%s | Unmatched:%s\n" % (len(summL2),bcount,dcount))
        fh_log.write("Matched %s:%s (summary1:summary2)\n" % (acount,ecount))


    return resList

def compare_writer(resList):
    '''
    writes compared results to a file
    '''

    outfile = "./%s/phasis.compare.txt" % (comp_folder)
    fh_out  = open(outfile,'w')
    fh_out.write("Name.a\tChr.a\tStart.a\tEnd.a\tP-val.a\tIdentifier.a\tBest k-val.a\tPhasi ratio.a\tMax Tag Ratio.a\tName.b\tChr.b\tStart.b\tEnd.b\tP-val.b\tIdentifier.b\tBest k-val.b\tPhasi ratio.b\tMax Tag Ratio.b\tMatch Ratio\tMatch Start coords [a:b]\tMatch Length\n")

    acount = 0
    for ent in resList:
        # print(ent)
        phasinfo    = [str(x) for x in ent[0]]
        matratio    = str(ent[1])
        matcoords   = str(ent[2])
        matlen      = str(ent[3])
        fh_out.write("%s\t%s\t%s\t%s\n" % ("\t".join(phasinfo),matratio,matcoords,matlen))
        acount+=1

    fh_out.close()

    print("Entries written: %s" % (acount))


    return outfile

def usagedata():
    '''
    collects usage data for improments to script
    '''

    print("\n#### Fn: usagedata ###########################")
    # Open a plain text file for reading.  For this example, assume that
    # the text file contains only ASCII characters.

    ### Header
    auser           = getpass.getuser()
    msg             = MIMEMultipart()
    msg['Subject']  = "phasis | '%s' | runid - %s" % (auser,datetime.datetime.now().strftime("%m_%d_%H_%M"))
    msg['From']     = "phasworks@gmail.com"
    msg['To']       = "kakrana@udel.edu"
    msg.preamble    = "Usage data for developer"

    attachlist = ["phasis.set","listdir.txt"]
    if os.path.isfile(setFile):
        pass
    else:
        print("---Settings file 'phasis.set' not found in current directory")
        print("---Please copy it to same directory as script and rerun")
        sys.exit()

    ### List of files
    # print("Attaching list of directory")
    alist   = os.listdir(os.getcwd())
    afile   = "listdir.txt"
    fh_out  = open(afile,'w')
    fh_out.write("\n\n$ls\n")
    for i in alist:
        fh_out.write("%s\n" % (i))
    fh_out.close()

    ### Attachements 
    for afile in attachlist:
        # print("Attaching file:%s" % (afile))
        with open(afile) as fp:
            atext = MIMEText(fp.read())
        msg.attach(atext)

    ### Send 
    # Send the message via our own SMTP server.
    # print("Sending settings file")
    s = smtplib.SMTP('localhost')
    s.send_message(msg)
    s.quit()

    return None

#### MAIN ###################
#############################
def main():
    
    checkDependency()
    if args.debug == "T":
        global smtplib
        global MIMEText
        global MIMEMultipart
        import smtplib 
        from email.mime.text import MIMEText
        from email.mime.multipart import MIMEMultipart
        usagedata()

    if args.mode == "merge": ## Default mode
        ### Collapser #########################################
        libs                    = readSet(setFile)
        pcutoff,pval_sorted     = pvaluereader()
        temp_folder,clustfile   = prepare(pcutoff,libs,res_folder)

        #### Write to memeory
        fh_mem = open("%s/%s" % (res_folder,memFile),'w')
        fh_mem.write("@phase:%s\n" % (phase))
        fh_mem.write("@pval:%s\n" % (pcutoff))
        
        global overlapCutoff
        if runType == 'G' or runType == "S":
            overlapCutoff = 0.05 ## = 0.25 for genomic and 0.50 for ncRNAs
        else:
            overlapCutoff = 0.10 ## = 0.25 for genomic and 0.50 for ncRNAs

        if fileType == 'L':
            print('\nList files selected for analysis - Converting them to readable format')
            fls     = glob.glob(r'%s/*.PARE.validation.list' % (temp_folder))
            # print ('Here are the files that will be converted:',fls,'\n')
            print ('Total files to analyze: %s' % (len(fls)))

            ### Prepare first file for comaprision
            firstfile   = fls[0]
            firstlist   = listConverter(firstfile,pcutoff)                          ## List of PHAS from file
            firstgrps   = groupPHAS(firstlist)                              ## PHAS list grouped on chr/scaffold or transcripts

            #### Test - serial selfMerge function
            # firstmergeL = []
            # for i in firstgrps:
            #     # print("## Group",i)
            #     akey,aval = selfMerge(i) ### Provides a merged list of PHAS for chr/scaffold and Trans
            #     firstmergeL.append((akey,aval))

            ### Parallel mode
            firstmergeL= PPResults(selfMerge,firstgrps)                     ## PHAS list made non-redundant from diff conf. levels
            firstmergeD= dict((i[0], i[1]) for i in firstmergeL)            ## Dict. of PHAS list based on chr/scaffold and trans
            collapsedL = listTocollapsed(firstmergeL)                             ## Collapsed list format incase there is just one library

            ### Start the comparision
            compareflag = False     ## Flag used to decide if it comaprision between first two files or afile and collapsed list
            totalfls    = len(fls)
            flcount     = 2         ## Two files are comapred in first loop so count starts from 2 
            for fl in fls[1:]:
                ## Prepare chr/scaffold/transcrpt specific dict of PHAS
                print("\n#### Comparing %s/%s file" % (flcount,totalfls))
                print("#### File being compared:%s" % (fl))
                fllist          = listConverter(fl,pcutoff)                         ## List of PHAS from file
                flgrps          = groupPHAS(fllist)                      ## PHAS list grouped on chr/scaffold or transcripts
                flmergeL        = PPResults(selfMerge,flgrps)            ## PHAS list made non-redundant from diff conf. levels
                flmergeD        = dict((i[0], i[1]) for i in flmergeL)   ## Dict. of PHAS list based on chr/scaffold and trans

                if compareflag == False:
                    ## No comparision made yet this is the first and uses list directly from files
                    rawinputs       = inputMaker(firstmergeD,flmergeD) ### File names are encoded in sub-lists
                    
                    ## Test - serial mergePHAS function
                    # collapsedL = [] ### List to hold resultsfrom all chr/scaffold and trascriptome
                    # for aninput in rawinputs:
                    #     ares = mergePHAS(aninput)
                    #     collapsedL.append(ares)
                    #     # print(collapsedL)

                    collapsedL      = PPResults(mergePHAS,rawinputs)## Compare PHAS from file with PHAS from collapsed list
                    compareflag     = True
                    flcount         += 1

                else:                
                    xD              = {}                            ## Set to empty before updation, just to be sure
                    xD              = collapsedToDict(collapsedL)   ## Reformat the collapsed results to dict from files
                    rawinputs       = inputMaker(xD,flmergeD)       ## Prepare for comaprision between both
                    collapsedL      = []                           ## Set to empty before updation, just to be sure
                    
                    ## Test - serial mergePHAS function
                    # for aninput in rawinputs:
                    #     ares = mergePHAS(aninput)
                    #     collapsedL.append(ares)
                    #     # print(collapsedL)

                    collapsedL      = PPResults(mergePHAS,rawinputs)## Compare PHAS from file with PHAS from collapsed list
                    flcount         += 1

        elif fileType == 'C':
            ### Function for future updates
            listConverter(temp_folder,fileType)
            main_dict = compare(temp_folder,fileType)
        else:
            ### Will not be user ever
            print("Unknown filetype - developers setting is off ")
            print("Download the author's orginal version from: https://github.com/atulkakrana/phasTER/releases")
            sys.exit()
            # main_dict = removeRedundant(temp_folder,pcutoff,fileType,overlapCutoff)
        
        collapsedfile,collapsedLfile = writer_collapse(collapsedL,pcutoff)
        fh_mem.write("@collapsedfile:%s\n" % (collapsedfile))
        fh_mem.write("@collapsedlist:%s\n" % (collapsedLfile))


        ### Summarizer ###########################################
        
        ## Read phasiFile
        phasList,phashead = PHASreader(collapsedfile)
        time.sleep(1)

        ## Sanity check
        if not phasList:
            print("** No PHAS loci or transcripts found at p-value %s" % (pcutoff))
            print("** Please use a lower p-value, here are some valid inputs: %s" % (", ".join(str(x) for x in pval_sorted)))
            sys.exit()

        ####################
        ## Get Clusters ####
        fh_in           = open(clustfile,'r')
        clusters        = fh_in.read().split('>')
        totalphas       = len(phasList)
        nphasL          = [i+1 for i in range(len(phasList))] ## an array of numbers of rdifferent PHAS for print statement
        rawinputs       = [(ent,clusters,nphas,totalphas) for ent,nphas in zip(phasList,nphasL)]
        fh_in.close()
        
        #### Serial - Test
        # masterL         = []
        # for i in rawinputs:
        #     entL = getClust(i)
        #     masterL.append(entL)

        #### Parallel - Original
        masterL = PPResults(getClust,rawinputs)

        # print("Clusters fetched for %s PHAS" % (len(masterL)))

        ### Decode the inputs to three separate lists
        resList     = []
        resList2    = []
        finalClustL = []
        for aninput in masterL:
            tempList,tempList2,tempClustL = aninput
            resList.append(tempList)   ## Unlisting the entries i.e. removing square brackets [()]
            resList2.append(tempList2) ## Unlisting the entries i.e. removing square brackets [()]
            for i in tempClustL:
                finalClustL.append(i)

        # print("Length of resList:%s | resList2:%s | finalClustL:%s" % (len(resList),len(resList2),len(finalClustL)))
        print("Total phas loci: %s | Clusters fetched: %s\n" % (len(phasList),len(resList)))
        ####################
        ####################


        if runType == 'T' or runType == 'S':
            allphasiFile = allphasiWriter(clustfile,resList2)
        else:
            # print("File with all PHAS will not be generated for this 'runType'")
            pass

        ## Prepare dictionary of tag count files for abundance queries
        dictList    = []
        indexList   = list(range(len(libs)))

        ### Serial - Use if you get error making dictionary from huge file > 3.5 GB (can be automated)
        # for anindex in indexList:
        #     adict = readFileToDict(anindex)
        #     dictList.append(adict)

        #### Parallel
        dictList    = PPResults(readFileToDict, indexList)

        ## Write the summary
        clustfile               = clustWriter(finalClustL,pcutoff)
        phasifile,summaryfile   = writer_summ(clustfile,resList,dictList,pcutoff)
        fh_mem.write("@summaryfile:%s\n"    % (summaryfile))
        fh_mem.write("@phasifile:%s\n"      % (phasifile))
        fh_mem.close()


        ##################
        #### ANNOTATE ####

        if args.gtf:
            #### Parse GTF ###
            featureFile     = args.gtf
            if annomode     == 1:   ## PASA GTF for overlapping transcripts
                gtfList     = gtfParser(featureFile)
            elif annomode   == 2 or annomode == 3 or annomode == 5 or annomode == 6: ## Trinity or Rocket GTF for overlapping transcipts
                gtfList     = gtfParser2(featureFile)
            elif annomode   == 4: ## Trinity and Rocket both for overlapping transcipts
                gtfList1    = gtfParser2(featureFile)  ## Both lists can be identified with aflag feature
                gtfList2    = gtfParser2(featureFile2) ## Both lists can be identified with aflag feature
                gtfList     = gtfList1 + gtfList2
            else:
                print("Please input correct mode in user settings - script will exit now")
                sys.exit()
            
            resFile     = overlapChecker(resList,gtfList,pcutoff)

        ### Prepare for revFerno and cleanup unwanted files
        # shutil.copy("./%s" % (setFile), "./%s" % (res_folder))
        if cleanup == 1:
            cleaner()

        print ("\n####'phasmerge' took", round(time.time()-start,2),"seconds")
        print ("#### Please see '%s' folder for results\n" % (res_folder))
        print("#### 'phastrigs' can be run by command: python3 phastrigs -mode auto -dir %s -mir your-miRNA-filename\n" % (res_folder))

    elif args.mode == "compare":

        ##################
        #### Prepare
        shutil.rmtree("%s" % (comp_folder),ignore_errors=True)
        os.mkdir("%s" % (comp_folder))

        log_file    = "%s/compare.log" % (comp_folder) 
        fh_log      = open(log_file,'w')
        fh_log.write("Date:%s\n" % ((datetime.datetime.now().strftime("%m_%d")) ))
        fh_log.write("Time:%s\n" % ((datetime.datetime.now().strftime("%H_%M")) ))


        ###################
        #### Read Summaries
        summL1,summD1,summfile1 = summparse(args.dir)
        summL2,summD2,summfile2 = summparse(args.dir2)
        fh_log.write("Summary1:%s\n" % (summfile1))
        fh_log.write("Summary2:%s\n" % (summfile2))
        fh_log.close()

        ##################
        #### Compare
        resList = compare(summL1,summD1,summL2,summD2,log_file)

        ##################
        #### Write Results
        compfile = compare_writer(resList)

        print ("\n####'phasmerge' took", round(time.time()-start,2),"seconds")
        print ("#### Please see '%s' folder for results\n" % (comp_folder))
        # print("#### 'phastrigs' can be run by command: python3 phastrigs -mode auto -dir %s -mir your-miRNA-filename\n" % (res_folder))

    else:
        sys.exit()

if __name__ == '__main__': 
    if cores == 0 :
        nproc = int(multiprocessing.cpu_count()*0.90)
    else:## As mannually entered by the user
        nproc = int(cores)
    start = time.time()
    main()
    sys.exit()


#### Change Log #############
#############################

## v01 -> v03
##Added functionality to compare between similar loci and select the longest one

## v03 -> v04 (Critical bug fixed)
## Corrected bug with 'key' - Range was being made one integer less as last integer is never counted and with every run (file) there is change in key (one integer less) it is
#### Therefore key from first run not matched woth earlier and all loci being retained

## v04-> v05 (Critical bug fixed)
## Bug with main key remade fixed

## v05 -> v06
### 1. If a lower cutoff is used than only entries pertaining to that cutoff should be analyzed, as phas are redundant between different cutoffs
###Problem observed was that 'Key not found' and reason for that is because Key has already been removed during first instance and later
###Instances of same loci but lower cutoff find main key missing from main dictionary thats why error.
###Only entries for input cutoff will be analyzed
### 2. capture loci that are different in two given files

## v06 -> V07
### Ping format added to compare files from different tissues
### Added library name in log file
### Added new files for each librray with clusters added - good to identify unique clusters i.e PARE validated and non-valiadted

###V08 -> v09 (Feb 2014)
##Added functionality to work on three file types and distinguish on the basis of switch - Working OK

## v09 -> v10 - Critical speed bug fixed
### The bug was related to appending chr to the start and end abd comparing the range. In cases when there is change in number of digits between start and end. The range for matching was generated incorrectly long and leads to freezing of script.
## In this version chromosome is not appended to start and end instead a check is doen for same chromosme and than range is generated, folllowed bu comparing ranges for overlap
## REgression - In this version I was getting error related to key-not found which is kind of historical error. Maybe in past (v04 -> v05) it was addressed by adding 1 to the end co-ordinates as while opening a range last number is not counted
## So after remaking the mainkey the last number as absent is not included in the key or something like that - For future 'Key not found' errors please check this first that
## Added p-val in name

## v1.0 -> v1.01
## Introduced functionality to combine results from csv files - this will help comparing old files with new files to identify new phased loci that needs manual curation
## Removed header as output of list converter as it's not generated always

## v1.0 -> v1.02
## Added overLap ratio for different regions in main()
## Removed chr_id as integer type as ncRNA names are not integers
## Modified main_key remade with additional replace of replace("'","") to make it functional for string chromosome/transcript IDs

## v1.03 -> 1.03
## Added a temp folder to write all "Uniq" files not working directory is not filled up with all less used files

## v1.03 -> v1.04
## Changed the RemoveRedundant module to be compatible with transcripts with colplex header like PacBio
## Changed the RemoveRedundant module so that no mainkey needs to be remade, instead it's stored as last value 

## v1.04 -> v1.05 [stable]
## "compare" module needs updated with changes made in v1.04 

## v1.06 -> v1.06 [stable]
## Fixed a syntax waring about gloabal variable 'overlapCutOff' declaration

## v1.06 -> collapser
## Find matching p-value less then user specified, select the highest size file name, extract uniq p-values, abd identify the one below specified cutoff
## Read settings file to get the library names, make working directory, and provide name of this folder (collapsed)
## Command line for cutoff with recomended cutoff

## v1.09 -> v1.10
## Added Reza modules to fetch library wise abundances
## Tested Reza module using the tagcount files downloaded from server with norm abudnaces and they match with the internal version

## v1.10 -> v1.11
## Input files with *.fas extension had *.txt extension in list files so these were not matching with the "libs". In prepare function
#### a new list of libs is made w/o file extensions and this is used to match with list and cluster files with just the library name

## v1.10 -> v1.12 [major]
## Updated runTYPEs - G,T and S modes
## Added a memory file which is read by revFerno
## overlapCutoff updated to use runType flags

## v1.12 -> v1.14 [major]
## Added faster algorithm for removing redudannt loci at both self (within the file) and global level (when comapring with others) -faster
## Added uniq alphabet to keys made for comarision between two lists so that PHAS with same coords between both files doesn't 
#### end up having same key which add common key to negtive list and updated the dict with same key - Such PHAS was being deleted
## New algorithm paralellized - faster
## Now includes all PHAS loci less then cutoff  
## Loci from every file, of different confidence levels are now uniques before comparing to other PHAS - Best loci and coords
## Fixed comarision based on length it was reversed before - best loci and coords
## This version will give best PHAS coordinates, compltely unique PHAS in final list
## all cluster file for p-values less then the used are picked up, because collapser usesfiles from 
#### all diffrent confidence levels. Also the matchThres for cluster and PHAS is increased to 0.99 because
#### now files for all confidence are available so cluster should match at high cutoff, reducing number the
#### total number of clusters picked up.
## closed lots of print statements in

## v1.14 -> v1.2
## Fixed error writing phasiRNAs for transcriptome mode and one more minor error
## CLuster files for all different confidence levels are included to fetch phasiRNAs abudnaces
#### It's because of this matchThres has been raised to 0.99 - To find exact same cluster 
#### for PHAS and reducing the number of cluster to scan. This also reduces noise in phasiRNA
#### abudnaces added due to presence of phasiRNAs from several clusters

## v1.2 -> v1.21
## Added checks for available p-value cutoffs in case they are lesser then hard default
## Added method find best pcutoff from available values if user doesn't specified
#### If this pcutoff is below recommneded confidence level then best available p-value is recommended
## MAde changes to accomodate pcutoff instead of default argument (1e-05) which is now empty 

## v1.21 -> v1.22 [major]
## Fixed error where comma-seprated @userlibs has a an empty entry, like a stray comma in end. Script doesn't ends in that case 
#### I though its because of many libraries being used and results not reported back i.e. multiprocessing issue but I was wrong
## Added a cluster writer, which writes matching clusters from all libs
## Fixed a major issue where in main(), for second, third and so on .... files the listconverted results from first file
### were being used instead of from the corresponding files
## Fixed an issue if transcipt PHAS was not present for both list to be compared. If alist is empty then all element ob blist are 
### added to merged list. And if blist is empty the for loop to comapre doen't even used.

## v1.21 -> v1.23 [major]
## "collapser" renamed to "phasmerge"
## Renames phaser.set and collapser.mem to phasworks.set and phasmerge.mem respectively
## "GetClust" functiona parallelized
## Added "compare" mode - Added additional match coords and lengths
## Added "-gtf" match functions
## Added dependency checks
## "Safe search" implemented, default cutoff are in developer settings area

## v1.23 -> v1.24
## Collect usage data for impromenets to scripts and for better debugging of reported issues
## Updated safesearch filters
## Moved sqlite3 import to args.gtf

## v1.24 -> v1.25 
## Fixed - Annotate functions being called even without a GTF file
## Reduced continuous phase positions cutoff to 10, and added phaslen filter to safesearch
## Closed a few print statements
## Added final prints for compare mode and added summary to comapre log file

## v1.25 -> v1.26 [stable]
## Optimized safesearch parameters by adding libwise abundance cutoff
## Added code to parse FASTA file from phasdetect with abudnaces in dictList module
## Changed the way files were handled in earlier function for tagcount file in dictList module
## Added sanity to checks to readset function
## Closed print statments for copying results, cluster and merging of these
## 'overlapCutoff' variable was named redundantly changed one for annotations to 'ntoverlapCutoff'

## v1.26 -> v1.27 [stable]
## Renamed to PHASIS

########################################
## PUBLIC RELEASE
## Turn OFF debug mode

#######################################
#### POSSIBLE PENDING ISSUES###########
## None