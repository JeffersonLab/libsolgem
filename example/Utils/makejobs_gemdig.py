#!/usr/bin/python

import sys
import os
import datetime
import math

if len(sys.argv) != 7:
    print "Supply macro and number of jobs;"
    sys.exit()
#    print " you may supply a preinit macro as a last argument"

#macro = sys.argv[1][:-4]
cwd = os.getcwd() 

#if len(sys.argv) == 4:
#    macro_pre = sys.argv[3][:-4]

fspec = int(sys.argv[1])
nsig = int(sys.argv[2])
infile_sig_prefix = sys.argv[3]
nmax = int(sys.argv[4])
ninsplit = int(sys.argv[5])
nbacktoadd = int(sys.argv[6])

#kludge... dunno how to not hardcode this yet... :/
infile_bkgd_prefix = ["/work/halla/sbs/efuchey/gep12_beam_bkgd_20181008", "/work/halla/sbs/efuchey/gep12_beam_bkgd_20181012", "/work/halla/sbs/efuchey/gep12_beam_bkgd_20181018", "/work/halla/sbs/efuchey/gep12_beam_bkgd_20181022", "/volatile/halla/sbs/efuchey/gep12_beam_bkgd_20181028"]

bkgd_str = "nobkgd"
if nbacktoadd>0:
    bkgd_str = "bkgd"
    
det_suf = ""
if fspec == 1:
    det_suf = "bbgem"

if fspec == 3:
    det_suf = "ft"

if fspec == 4:
    det_suf = "fpp"



runfiletxt="""#!/bin/csh
"""

print( str(nmax)+" events to split by bunch of "+str(ninsplit) )
nsplits = int(math.ceil(nmax/ninsplit))
print( "file to split into "+str(nsplits)+" splits" )

mytime = datetime.datetime.today()
nbg = 0
fbgf = 0
for i in range(0, nsplits) :
    suffix = "dig_"+det_suf+"_"+bkgd_str+'_'+mytime.strftime("%Y%m%d_%H")+'_'+str(i)
    suffix2 = "dig_"+det_suf+"_"+bkgd_str+'_'+mytime.strftime("%Y%m%d_%H")
    suffix3 = str(i)
    
    n0 = i*ninsplit
    n1 = (i+1)*ninsplit
    if fbgf+ninsplit*nbacktoadd/2 > 10000:
        nbg += 1
        fbgf = 0
    
    print( fbgf )
    if nbg>4:
        nbg = 0
    
    print 'Creating job ', suffix

    runfile = open("runjob_"+suffix+".sh", 'w')
    runfile.write(runfiletxt)
    runfile.write("\n")
    runfile.write("source .cshrc\n")
    runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/db_run.dat .\n")
    runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/db_ratedig.dat .\n")
    runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/db_generalinfo_"+det_suf+".dat .\n")
    runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/db_g4sbs_"+det_suf+".dat .\n")
    runfile.write("time ./DigiPass "+str(fspec)+" "+str(nsig)+" "+infile_sig_prefix+" "+str(n0)+" "+str(n1)+" "+str(nbacktoadd)+" "+infile_bkgd_prefix[nbg]+" "+str(fbgf)+" > run.out\n")
    runfile.write("mkdir /volatile/halla/sbs/efuchey/gemdigi_"+suffix2+"\n")
    #runfile.write("mkdir /work/halla/sbs/efuchey/"+suffix2+"\n")
    runfile.write("mv digitized_"+det_suf+"_"+bkgd_str+".root /volatile/halla/sbs/efuchey/gemdigi_"+suffix2+"/digitized_"+det_suf+"_"+bkgd_str+"_"+suffix3+".root \n")
    #runfile.write("mv beam_bkgd.root /work/halla/sbs/efuchey/"+suffix2+"/beam_bkgd_"+suffix3+".root\n")
    #runfile.write("jput beam_bkgd.root /mss/home/efuchey/"+suffix2+"/beam_bkgd_"+suffix3+".root\n")
    runfile.write("mv run.out /volatile/halla/sbs/efuchey/run_"+suffix+".out\n")
    runfile.close()

    subfile = open("gemdig_"+suffix+".jsub", 'w')
    
    subfile.write("PROJECT: sbs\n")
    #subfile.write("NODE_TAG: farm14\n")
    subfile.write("TRACK: simulation\n")
    #subfile.write("TRACK: debug\n")
    subfile.write("TIME: 2160\n")
    subfile.write("OS: centos7\n")
    subfile.write("JOBNAME: gemdig_"+suffix+"\n")
    subfile.write("MEMORY: 3000 MB\n")
    subfile.write("COMMAND: csh runjob_"+suffix+".sh\n")
    subfile.write("OTHER_FILES:\n")
    subfile.write("/u/home/efuchey/.cshrc\n")
    subfile.write("/u/home/efuchey/g4work/Tracking/libsolgem/test/runjob_"+suffix+".sh\n")
    subfile.write("/u/home/efuchey/g4work/Tracking/libsolgem/test/DigiPass\n")
    #subfile.write(cwd+'/'+macro+".mac\n")
    #if len(sys.argv) == 4:
    #   subfile.write(cwd+'/'+macro_pre+".mac\n")
    subfile.close()

    #os.system("/site/scicomp/auger-slurm/bin/jsub gemdig_"+suffix+".jsub")
    os.system("jsub gemdig_"+suffix+".jsub")

    fbgf+= ninsplit*nbacktoadd/2



#Ex: To send of farm (conveniently)
# python2 makejobs_gemdig.py 3 1 /volatile/halla/sbs/efuchey/gep12_elastic_sig_20180920_21 40000 10000 0
# python2 makejobs_gemdig.py 3 1 /volatile/halla/sbs/efuchey/gep12_elastic_sig_20180920_21 10000 1000 10






