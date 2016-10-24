#!/usr/local/bin/python3

## This script uses the phasiRNA scripts from SVN So, make sure your svn has Pingchuan phasi-prediction scripts
## Author: Atul Kakrana kakrana@udel.edu

#### FUNCTIONS ###########################################

import os,sys,subprocess,multiprocessing,time,getpass,shutil,hashlib,datetime
from multiprocessing import Process, Queue, Pool
import os.path
# import mysql.connector as sql

#### USER SETTINGS ########################################


### Settings file
setFile     = "phaser.set"
memFile     = "phaser.mem"
res_folder  = "phased"
# allowedHost = ['tarkan.ddpsc.org']
# allowedUser = ['kakrana','suresh','kun','patel','rkweku']


## Genome - Mandatory
# runType            = 'N'                               ## Run on whole genome (Y) or transcriptome file (N) - If (N) then provide index of transcriptome
# index      = "../../index/merged.index"        ## If geno = 'Y' then index path from $ALLDATA; if geno = 'N' then prepare index of your file 
#                                                     ## using bowtie-build command like bowtie-build GENOMEFILE.fa NAME_FOR_INDEX

# ## sRNA - Mandatory
# db              = 'RICE_sbsQIFA_sRNA'               ## sRNA DB
# fetchLib     = 'N'                               ## (Y): Get IDs for all libs in DB (N): If you want to run on specific libs than 'N' and specifiy libs below
# # userLibs = [(2435,),(2436,),(2437,),(2438,),(2439,),(2440,),(2441,),(2495,),(2496,),(2497,),(2498,),(2499,),(2500,),(2501,)] ## Used only if fetchLib == 'N'
# userLibs = [(2509,),(2510,),(2511,),(2512,),(2513,),(2514,),(2515,),(2516,),(2517,),(2518,),(2519,),(2520,),(2521,),(2522,),(2596,),(2597,),(2598,),(2599,),(2600,),(2601,),(2602,)] ## Used only if fetchLib == 'N'
# phase           = 21                                ## Phase to use for prediction

## Degradome - Optional ####################################
deg             = 'N'                               ## Use Degradome validation, IF yes enter PARE db in line below
PARE            = 'GuturGu'                         ## If deg = 'Y' then File for degradome analysis

## ADVANCED SETTINGS #######################################
cores           = 0                                 ## 0: Most cores considered as processor pool | 1-INTEGER: Cores to be considered for pool
# nthread         = 3                               ## Threads perprocess
# server          = "tarkan.ddpsc.org"                ## Server to use to fetch library information and smallRNA libraries
# perl            = "/usr/local/bin/perl_5.18"        ## Josh updated the perl on Tarkan and its not ready yet for PHAS script FORK is missing and somemore modules -Check with Pingchuan help
perl            = "perl"
Local           = 3                                 ## [0]: Files in directory [2]: Get the libs from $ALLDATA with raw reads 
                                                    ## [3] Get library from srna db with reads filtered on number of hits

#############################################################
#############################################################

def checkuser():
    '''
    Checks if user is authorized to use script
    '''
    print ("#### Verifying User Authorization ################")
    auser = getpass.getuser()
    print("Hello '%s' - If you face any issues, then please report: https://github.com/atulkakrana/phasTER/issues \n" % (auser))
    # if auser in allowedUser:
    #     print("Hello '%s' - Issues need to be reproted: https://github.com/atulkakrana/phasTER/issues \n" % (auser))
    # else:
    #     print("YOU ARE NOT AUTHORIZED TO USE DEVELOPMENTAL VERSION OF 'PHASER'")
    #     print("Contact 'Atul Kakrana' at kakrana@gmail.com for permission\n")
    #     sys.exit()
    
    return None

def checkHost(allowedHost):

    '''
    Checks if Phster is allowed at this server
    '''
    print ("#### Pre-run checks #########################")
    f = subprocess.Popen("hostname", stdout=subprocess.PIPE,shell= True)
    output,err = f.communicate()
    #print (output.decode("ascii"))
    
    host = output.decode("ascii")
    print ('--Current host:',host.strip('\n'))
    
    ## DO not turn OFF this 'for' loop as that given an error while matching current host with allowedHost - Reason Unknown
    # print ('Allowed Hosts:')
    # for host in allowedHost:
    #     print (host)
    print("--Allowed hosts: %s" % (','.join(x for x in allowedHost)))
    
    
    if str(host.strip('\n')) in allowedHost:
        print("--phasTER is supported on this server - good to go!!!\n")
        pass
    
    else:
        print("--phasTER is not tested on this server")
        print("--Run your analysis at any of these servers:%s" % (','.join(x for x in allowedHost)))
        print("--Script will exit now\n")
        sys.exit()

    return None

def readSet(setFile):
    '''
    Read and parse external settings file
    '''

    if os.path.isfile(setFile):
        pass
    else:
        print("---Settings file 'phaser.set' not found in current directory")
        print("---Please copy it to same directory as script and rerun")
        sys.exit()

    print("#### Fn: Settings Reader ####################")
    
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
                
                ### Extract values ######### 
                if param.strip() == '@runType':
                    global runType
                    runType = str(value.strip())
                    print('User Input runType:              ',runType)

                elif param.strip() == '@reference':
                    global reference
                    reference = str(value.strip())
                    print('User Input reference location:   ',reference)
                
                elif param.strip() == '@index':
                    global index
                    index = str(value.strip())
                    print('User Input index location:       ',index)

                elif param.strip() == '@userLibs':
                    global userLibs
                    libs = list(map(str,value.strip().split(',')))
                    print('User Input Libs:                 ',libs)

                elif param.strip() == '@libFormat':
                    global libFormat
                    libFormat = str(value.strip())
                    print('User Input to auto fetch libs:   ',libFormat)

                elif param.strip() == '@genoFile':
                    global genoFile
                    genoFile = str(value.strip())
                    print('User Input genoFile:             ',genoFile)

                elif param.strip() == '@phase':
                    global phase
                    phase = int(value.strip())
                    print('User Input for phase length:     ',phase)
                
                elif param.strip() == '@path_prepro_git':
                    global phaster_path
                    phaster_path = str(value.strip()).rstrip("/")+"/phaster"
                    # phaster_path = str(value.strip()).rstrip("/")+"/core"
                    print('User Input for phaster path:     ',phaster_path) 

            else:
                #print("Missed line:",line)
                pass
    
    # sys.exit()
    return libs

def ConnectToDB(server):
    
    ##infile values are '0' when you dont want to pulaod data from local file and '1' when you wish to upload data by local file
    ##EX:con=sql.connect(host= server, user='kakrana', passwd='livetheday', local_infile = infile)
    ##Now later in script you can
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ','")
    
    print ('\nTrying to connect to mySQL server on %s' % (server))
    # Try to connect to the database
    try:
        con=sql.connect(host= server, user='kakrana', passwd='xxxxxxx')###local_infile = 1 not supported yet so a table has to be updated on row basis
        print ('Connection Established\n')

    # If we cannot connect to the database, send an error to the user and exit the program.
    except sql.Error:
        print ("Error %d: %s" % (sql.Error.args[0],sql.Error.args[1]))
        sys.exit(1)

    return con

## Get lib ids in the sRNA DB
def GetLibs(con,db):
    ##Function that gets just the list names, required to run script in parts.
    cur = con.cursor()
    cur.execute('select distinct(lib_id) from %s.library' % (db))
    libs = cur.fetchall()
    #print (libs)
    print ('\nTotal number of sRNA libraries found: %s\n' % (len(libs)))
    
    return libs###

## Deprecated
def PHASBatch(con,libs,runType,index,deg):
    
    #os.mkdir('./%s' % (lib))
    #output_path = './%s' % (lib)
    
    for lib in libs:
        print (lib)
        cur = con.cursor()
        cur.execute('SELECT processed_path FROM master.library_info where lib_id = %s' % (lib))
        path = cur.fetchall()
        #print(path[0][0])
        
        pro_file = path[0][0].replace('$ALLDATA', '/alldata')###Processed sRNA file
        out_file = '%s.txt' % (lib)
        rl = str(phase)
        nproc2 = str(nproc)
        sRNAratio = str(75)
        print (pro_file)
        
        if runType == 'Y':###Uses Whole genome as input
            if deg == 'Y':
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/kakrana/phasiRNA_prediction_pipeline.ver.genome.pl", "-i", pro_file, "-q", PARE, "-f", "-t", sRNAratio, "-d", index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
            else:
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/kakrana/phasiRNA_prediction_pipeline.ver.genome.pl", "-i", pro_file,"-f", "-t", sRNAratio, "-d", index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
        
        else: ### Uses FASTA file of genes as input         
            #pipe =subprocess.Popen(["perl5.18", "-v"])
            if deg == 'Y':
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/kakrana/phasiRNA_prediction_pipeline.ver.MUL.pl", "-i", pro_file, "-q", PARE, "-f", "-t", sRNAratio, "-d", index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
            else:
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/kakrana/phasiRNA_prediction_pipeline.ver.MUL.pl", "-i", pro_file, "-f", "-t", sRNAratio, "-d", index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
                    
        
        if retcode == 0:
            pass
        else:
            print("Problem with Phasing script - Return code not 0")
            sys.exit()
        
    return lib

### sRNA Libraries are fetched from server
def TagAbundanceFile(con,db,libs):
    
        for alib in libs:##For all the libraries
            
            ## Check if file already exsits in directory - This saves a lot of time downloading the same file
            filePath = '%s.fas' % (alib)
            if os.path.isfile(filePath) == False:
                print ('\nPreparing sRNA reads file for library: %s' % (alib[0]))
                #print (lib[0])
                #print ('Caching tag and count information from server for PARE alib %s' % (alib[0]) )
                cur = con.cursor()
                cur.execute("SELECT tag,norm from %s.run_master where lib_id = %s AND (hits between 0 and 20)" % (db,alib[0]))
                lib_info = cur.fetchall()
                #print('These are the tags:',lib_info[:10])
                
                fh_out = open('%s.fas' % (alib), 'w')##Naming file with lib_ids name
                print ('Library cached, writing abundance file')
                tag_num = 1
                for ent in lib_info:## All the PARE tags in a library
                    #print (ent)
                    fh_out.write('%s\t%s\n' % (ent[0],ent[1]))
                    tag_num += 1
                    
                fh_out.close()
            else:
                print('tag abundance file exists for library: %s' % (alib))
                pass

## Phasing anlysis - New
def PHASBatch2(aninput):

    print ("\n#### Fn: phaser #########################")
    # print("\naninput\n",aninput)
    lib,runType,index,deg,nthread = aninput
    
    ### Sanity check #####################
    if not os.path.isfile(lib):
        print("** %s - sRNA library file not found" % (lib))
        print("** Please check the library- Is it in specified directory? Did you input wrong name?")
        print("** Script will exit for now\n")
        sys.exit()
    else:
        print("sRNA library located - Running phasing analysis")
        pass
    #####################################


    pro_file = lib ### sRNA input file 
    out_file = './%s/%s.txt' % (res_folder,lib.rpartition(".")[0]) ## Output file suffix
    
    rl          = str(phase)
    # nproc2 = str(nproc)
    nthread     = str(nthread)
    sRNAratio   = str(75)
    print(pro_file)

    if runType == 'Y':### Uses Whole genome as input
        if deg == 'Y':
            retcode = subprocess.call([perl, "%s/phasiRNA_prediction_pipeline.ver.genome.pl" % (phaster_path), "-i", pro_file, "-q", PARE, "-f", "-t", sRNAratio, "-d", index, "-px", out_file, "-rl", rl, "-cpu", nthread])
        else:
            if libFormat == "T":
                aformat = "t"
                retcode = subprocess.call([perl, "%s/phasiRNA_prediction_pipeline.ver.genome.pl" % (phaster_path), "-i", pro_file,"-f", aformat, "-t", sRNAratio, "-d", index, "-px", out_file, "-rl", rl, "-cpu", nthread])
            elif libFormat == "F":
                aformat = "f"
                retcode = subprocess.call([perl, "%s/phasiRNA_prediction_pipeline.ver.genome.pl" % (phaster_path), "-i", pro_file,"-f", aformat, "-t", sRNAratio, "-d", index, "-px", out_file, "-rl", rl, "-cpu", nthread])
            else:
                print("** Invalid '@libFormat' parameter value")
                print("** Please check the '@libFormat' parameter value in setting file")
                print("** F for FASTA format | T for tag-count format are the only acceptable values")
                print("** Script will exit now")
                sys.exit()
    
    else: ### Uses FASTA file of genes as input         
        if deg == 'Y':
            retcode = subprocess.call([perl, "%s/phasiRNA_prediction_pipeline.ver.MUL.pl" % (phaster_path), "-i", pro_file, "-q", PARE, "-f", "-t", sRNAratio, "-d", index, "-px", out_file, "-rl", rl, "-cpu", nthread])
        else:
            if libFormat == "T":
                aformat = "t"
                retcode = subprocess.call([perl, "%s/phasiRNA_prediction_pipeline.ver.MUL.pl" % (phaster_path), "-i", pro_file, "-f", aformat, "-t", sRNAratio, "-d", index, "-px", out_file, "-rl", rl, "-cpu", nthread])
            elif libFormat == "F":
                aformat = "f"
                retcode = subprocess.call([perl, "%s/phasiRNA_prediction_pipeline.ver.genome.pl" % (phaster_path), "-i", pro_file,"-f", aformat, "-t", sRNAratio, "-d", index, "-px", out_file, "-rl", rl, "-cpu", nthread])
            else:
                print("** Invalid '@libFormat' parameter value")
                print("** Please check the '@libFormat' parameter value in setting file")
                print("** F for FASTA format | T for tag-count format are the only acceptable values")
                print("** Script will exit now")
                sys.exit()

    if retcode == 0:
        pass
    else:
        print("** Problem with Phasing script - Return code not 0")
        sys.exit()
        
    return None

## Balance process according to core pool
def PPBalance(module,alist):
    #print('***********Parallel instance of %s is being executed*********' % (module))
    start       = time.time()
    ##PP is being used for Bowtie mappings - This will avoid overflooding of processes to server
    nprocPP     = round((nproc/int(nthread))) 
    if nprocPP  < 1:
        nprocPP = 1 ## 1 here so as to avoid 0 processor being allocated in serial mode
    else:
        pass

    print("nprocPP:%s" % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)

def optimize(nproc):
    '''
    dirty optimization of threads per library
    '''

    nlibs       = len(libs)
    ninstances  = int(nproc/nlibs) ### Number of parallel instances to use

    if ninstances > 3:
        nthread = ninstances
    else:
        nthread = 3

    print("\n#### %s cores reserved for analysis #########" % (str(nproc)))
    print("#### %s threads assigned to one lib #########\n" % (str(nthread)))
    time.sleep(1)

    return nthread 

## Generate rawInputs fpr PP balance
def inputList(libs,runType,index,deg,nthread):
    '''generate raw inputs for parallel processing'''

    rawInputs = [] ## An empty list to store inputs for PP
    for alib in libs:
        rawInputs.append((alib,runType,index,deg,nthread))

    # print("These are rawInputs:",rawInputs)

    return rawInputs

def indexBuilder(genoFile):
    
    
    print ("#### Fn: indexBuilder #########################")
    ### Sanity check #####################
    if not os.path.isfile(genoFile):
        print("'%s' reference file not found" % (genoFile))
        print("Please check the genomeFile - Is it in specified directory? Did you input wrong name?")
        print("Script will exit for now\n")
        sys.exit()
    else:
        print("Reference file located - Preparing to create index")
        pass
    #####################################

    ### Clean reference ################
    fastaclean,fastasumm = FASTAClean(genoFile,0)

    ### Prepare Index ##################
    print ("**Deleting old index 'folder' !!!!!!!!!!!**")
    print("If its a mistake cancel now by pressing ctrl+D and continue from index step by turning off earlier steps- You have 5 seconds")
    time.sleep(5)
    shutil.rmtree('./index', ignore_errors=True)
    os.mkdir('./index')
    
    genoIndex   = '%s/index/%s' % (os.getcwd(),fastaclean.rpartition('/')[-1].rpartition('.')[0]) ## Can be merged with genoIndex from earlier part if we use bowtie2 earlier
    print('Creating index of cDNA/genomic sequences:%s**\n' % (genoIndex))
    retcode     = subprocess.call(["bowtie-build", fastaclean, genoIndex])
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
        # print("Reference index prepared sucessfully")
        pass
    else:
        print("There is some problem preparing index of reference '%s'" %  (reference))
        print("Is 'Bowtie' installed? And added to environment variable?")
        print("Script will exit now")
        sys.exit()
    #####################################

    ### Make a memory ###################
    fh_out      = open(memFile,'w')
    refHash     = (hashlib.md5(open('%s' % (reference),'rb').read()).hexdigest()) ### reference hash used instead of cleaned FASTA because while comparing only the user input reference is available
    indexHash   = (hashlib.md5(open('%s.1.ebwt' % (genoIndex),'rb').read()).hexdigest())
    print("\n@genomehash:%s | @indexhash:%s" % (refHash, indexHash) )
    fh_out.write("@timestamp:%s\n" % (datetime.datetime.now().strftime("%m_%d_%H_%M")))
    fh_out.write("@genomehash:%s\n" % (refHash))
    fh_out.write("@index:%s\n" % (genoIndex))
    fh_out.write("@indexhash:%s\n" % (indexHash))

    print("Index prepared:%s\n" % (genoIndex))

    # sys.exit()
    
    return genoIndex

def FASTAClean(filename,mode):
    
    '''Cleans FASTA file - multi-line fasta to single line, header clean, empty lines removal'''

    ## Read seqeunce file
    fh_in       = open(filename, 'r')
    print ("PHASER uses FASTA header as key for identifying the phased loci")
    print ("Cleaning '%s' reference FASTA file" % (filename))
    
    ## Write file
    if mode == 0:
        fastaclean = ('%s/%s.clean.fa' % (os.getcwd(),filename.rpartition('/')[-1].rpartition('.')[0])) ## os.getcwd(),fastaclean.rpartition('/')[-1].rpartition('.')[0]
    elif mode == 1:
        fastaclean = ('%s/%s.clean.rev.fa' % (os.getcwd(),filename.rpartition('/')[-1].rpartition('.')[0]))
    elif mode == 2:
        fastaclean = ('%s/%s.clean.revcomp.fa' % (os.getcwd(),filename.rpartition('/')[-1].rpartition('.')[0]))
    elif mode == 3:
        fastaclean = ('%s/%s.clean.comp.fa' % (os.getcwd(),filename.rpartition('/')[-1].rpartition('.')[0]))
    elif mode == 4:
        fastaclean = ('%s/%s.clean.rna.fa' % (os.getcwd(),filename.rpartition('/')[-1].rpartition('.')[0]))
    elif mode == 5:
        fastaclean = ('%s/%s.clean.dna.fa' % (os.getcwd(),filename.rpartition('/')[-1].rpartition('.')[0]))
    else:
        print("Input correct mode- 0: Normal | 1: Seqeunces reversed | 2: Seqeunces reverse complemented | 3: Seqeunces complemented only")
        print("USAGE: cleanFasta.v.x.x.py FASTAFILE MODE")
        sys.exit()

    fh_out1     = open(fastaclean, 'w')

    fastasumm   = ('%s.summ.txt' % (filename.split('.')[0]))
    fh_out2     = open(fastasumm, 'w')
    fh_out2.write("Name\tLen\n")
    
    fasta       = fh_in.read()
    fasta_splt  = fasta.split('>')
    acount      = 0 ## count the number of entries
    empty_count = 0
    for i in fasta_splt[1:]:
        ent     = i.split('\n')
        name    = ent[0].split()[0].strip()
        seq     = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
        alen    = len(seq)

        if mode == 0:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        elif mode == 1:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq[::-1]))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        elif mode == 2:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq[::-1].translate(str.maketrans("TAGC","ATCG"))))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        elif mode == 3:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq.translate(str.maketrans("TAGC","ATCG"))))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        elif mode == 4:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq.translate(str.maketrans("TAGC","UAGC"))))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        elif mode == 5:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq.translate(str.maketrans("UAGC","TAGC"))))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        else:
            print("Please enter correct mode")
            pass

        acount+=1
    
    fh_in.close()
    fh_out1.close()
    fh_out2.close() 

    print("Fasta file with reduced header: '%s' with total entries %s is prepared" % (fastaclean, acount))
    print("There were %s entries found with empty sequences and were removed\n" % (empty_count))
    
    return fastaclean,fastasumm

def readMem(memFile):
    '''
    Reads memory file and gives global variables
    '''
    print ("#### Fn: memReader #########################")
    fh_in   = open(memFile,'r')
    memRead = fh_in.readlines()
    fh_in.close()

    for line in memRead:
        if line: ## Not empty
            if line.startswith('@'):
                line        = line.strip("\n")
                # print(line)
                akey,aval   = line.split(':')
                param       = akey.strip()
                value       = aval.strip()
                # print(param,value)

                if param == '@genomehash':
                    global existRefHash
                    existRefHash = str(value)
                    print('Existing reference hash:              ',existRefHash)

                elif param == '@indexhash':
                    global existIndexHash
                    existIndexHash = str(value)
                    print('Existing index hash:                  ',existIndexHash)

                elif param == '@index':
                    global index
                    index = str(value)
                    print('Existing index location:               ',index)
                
                else:
                    pass

    return None

def main(libs):
    # global con
    # con = ConnectToDB(server)

    ### Open the runlog
    runLog          = 'runtime_%s' % datetime.datetime.now().strftime("%m_%d_%H_%M")
    fh_run          = open(runLog, 'w')
    phaser_start    = time.time()

    ## Did user provided its index? If Yes Skip making memory files
    if not index:
        ### Check genome file and index
        if not os.path.isfile(memFile):
            print("This is first run - create index")
            tstart      = time.time()
            genoIndex   = indexBuilder(reference)
            tend        = time.time()
            fh_run.write("Indexing Time:%ss\n" % (round(tend-tstart,2)))
        else:
            currentRefHash = hashlib.md5(open('%s' % (reference),'rb').read()).hexdigest()
            readMem(memFile)
            if currentRefHash == existRefHash:
                print("Existing index matches the specified reference and will be used")
                genoIndex = index
                fh_run.write("Indexing Time: 0s\n")
            else:
                print("Existing index does not matches specified genome - It will be recreated")
                tstart      = time.time()
                genoIndex   = indexBuilder(reference)
                tend        = time.time()
                fh_run.write("Indexing Time:%ss\n" % (round(tend-tstart,2)))
    else:        
        genoIndex = index
        if not os.path.isfile("%s.1.ebwt" % (genoIndex)):
            print("** %s - User specified index not found" % (genoIndex))
            print("** Please check the value for @index parameter in settings file")
            print("** Is it in specified directory? Did you input wrong name?")
            print("** Script will exit for now\n")
            sys.exit()
        else:
            print("User specified genoIndex located and will be used")
            fh_run.write("Indexing Time: 0s\n")
            pass

    ### Start the phasing analysis ##########

    ### 1. Make Folders
    shutil.rmtree("%s" % (res_folder),ignore_errors=True)
    os.mkdir("%s" % (res_folder))

    ### 2. Run Phaser
    print('These are the libs: %s' % (libs))
    rawInputs = inputList(libs,runType,genoIndex,deg,nthread)

    #### Test - Serial Mode
    # for aninput in rawInputs:
    #     PHASBatch2(aninput)

    #### Original - Parallel mode
    PPBalance(PHASBatch2,rawInputs)

    #### close runLog
    phaser_end = time.time()
    fh_run.write("Analysis Time:%ss\n" % (round(phaser_end-phaser_start,2)))
    fh_run.close()

if __name__ == '__main__':

    ###Processors to use####
    if cores == 0 :## Use default 95% of processors as pool
        nproc = int(multiprocessing.cpu_count()*0.90)
    else:## As mannually entered by the user
        nproc = int(cores)

    ###############
    checkuser()
    # checkHost(allowedHost)
    libs        = readSet(setFile)
    nthread     = optimize(nproc)
    main(libs)
    print ('\n\nPhasing Analysis finished successfully')
    sys.exit()

########### CHANGE LOG ########
###############################

### Version 01 -> v02
### Added PARE switch
### Added sRNA ratio option
### Added option to specify libs

## v02 -> v03
## Added option to get libs from the server with hits filter
## COrrected bug in main(), repaced libs with userlibs for specific librarues part
## Perl location added as variable

## v03 -> v04
## Changed order of user settings to make them more clear
## Added functionality to check if the abundance file for library already exists in folder - Saves a lot of time

## v04 -> v05
## Added local mode to run on smallRNA files specified by user and present in current directory unlike fetching from DB or ALLDATA
## Simplfied user settings

## v05 -> v06
## Added extension to sRNA library files for easy copying

## v06 -> v065
## Fixed regresion introduced in v06

## v065 -> v07 [major][stable]
## Includes fix for iltering out tags with no hits, these are inlcuded now for libraries that have no genomes

## v070 -> v075 [major]
## Paralelization schema improved - Now paralelized three libraries together, only the analysis part is parallelized and not downloading part

## v075 -> v080
## Changed nProc input from "Y" to 0
## Fixed bug from v075 if fethLibs = 'Y', then also use libs were being used for raw inputs

## v08 -> v085
## moved scripts to svn/users/kakrana and updated the paths

## v085 -> v090
## Script run on local libraries
## Localization complete
## Added sanity checks
## Index made if change in genome detected, and reused if genome/referenc eis not changed


## TO-DO
## Add automatic index resolution
## Add functionality to share library folder, for 21-,24- and 22nt analysis
