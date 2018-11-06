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

filename = sys.argv[1]
print ( filename )
det_suf = sys.argv[2]
print ( det_suf )
bkgd = int(sys.argv[3])
nevent = int(sys.argv[4])
path_inputfile = sys.argv[5]
print ( path_inputfile )
nfiles = int(sys.argv[6])

nseg = 0
do_cuts = 0

bkgd_str = "nobkgd"
if bkgd>0:
    bkgd_str = "bkgd"

runfiletxt="""#!/bin/csh
"""
mytime = datetime.datetime.today()
for i in range(0, nfiles) :
    suffix = "ana_"+det_suf+"_"+bkgd_str+'_'+mytime.strftime("%Y%m%d_%H")+'_'+str(i)
    suffix2 = "ana_"+det_suf+"_"+bkgd_str+'_'+mytime.strftime("%Y%m%d_%H")
    suffix3 = str(i)
    
    print 'Creating job ', suffix

    runfile = open("runjob_"+suffix+".sh", 'w')
    runfile.write(runfiletxt)
    runfile.write("\n")
    runfile.write("source .cshrc\n")
    runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/db_run.dat .\n")
    #runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/db_ratedig.dat .\n")
    runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/db_generalinfo_"+det_suf+".dat .\n")
    runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/db_g4sbs_"+det_suf+".dat .\n")
    runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/db_sbs_"+det_suf+".tracker.dat .\n")
    runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/db_sbssim_cratemap.dat .\n")
    runfile.write("ln -s /u/home/efuchey/g4work/Tracking/libsolgem/test/sbssim.odef .\n")
    runfile.write("ln -s "+path_inputfile+"/"+filename+"_"+det_suf+"_"+bkgd_str+"_"+str(i)+".root "+filename+"_"+det_suf+"_"+bkgd_str+".root\n")
    print("time ./ReplayMCDig "+filename+" "+det_suf+" "+str(bkgd)+" "+str(nevent)+" "+str(nseg)+" "+str(do_cuts)+" > run.out\n")
    runfile.write("time ./ReplayMCDig "+filename+" "+det_suf+" "+str(bkgd)+" "+str(nevent)+" "+str(nseg)+" "+str(do_cuts)+" > run.out\n")
    runfile.write("mkdir /volatile/halla/sbs/efuchey/digana_"+suffix2+"\n")
    #runfile.write("mkdir /work/halla/sbs/efuchey/"+suffix2+"\n")
    runfile.write("mv "+filename+"_"+det_suf+"_"+bkgd_str+"_replayed_new.root /volatile/halla/sbs/efuchey/digana_"+suffix2+"/"+filename+"_"+det_suf+"_"+bkgd_str+"_replayed_"+suffix3+".root \n")
    runfile.write("cp db_sbs_"+det_suf+".tracker.dat /volatile/halla/sbs/efuchey/digana_"+suffix2+"/. \n")
    #runfile.write("mv beam_bkgd.root /work/halla/sbs/efuchey/"+suffix2+"/beam_bkgd_"+suffix3+".root\n")
    #runfile.write("jput beam_bkgd.root /mss/home/efuchey/"+suffix2+"/beam_bkgd_"+suffix3+".root\n")
    runfile.write("mv run.out /volatile/halla/sbs/efuchey/run_"+suffix+".out\n")
    runfile.close()

    subfile = open("digana_"+suffix+".jsub", 'w')
    
    subfile.write("PROJECT: sbs\n")
    subfile.write("NODE_TAG: farm14\n")
    subfile.write("TRACK: simulation\n")
    #subfile.write("TRACK: debug\n")
    subfile.write("TIME: 2160\n")
    subfile.write("OS: centos7\n")
    subfile.write("JOBNAME: gemdig_"+suffix+"\n")
    subfile.write("MEMORY: 1000 MB\n")
    subfile.write("COMMAND: csh runjob_"+suffix+".sh\n")
    subfile.write("OTHER_FILES:\n")
    subfile.write("/u/home/efuchey/.cshrc\n")
    subfile.write("/u/home/efuchey/g4work/Tracking/libsolgem/test/runjob_"+suffix+".sh\n")
    subfile.write("/u/home/efuchey/g4work/Tracking/libsolgem/test/ReplayMCDig\n")
    #subfile.write(cwd+'/'+macro+".mac\n")
    #if len(sys.argv) == 4:
    #   subfile.write(cwd+'/'+macro_pre+".mac\n")
    subfile.close()

    #os.system("/site/scicomp/auger-slurm/bin/jsub gemdig_"+suffix+".jsub")
    os.system("jsub digana_"+suffix+".jsub")



#Ex: To send of farm (conveniently)
# python2 makejobs_digana.py digitized ft 0 -1 /volatile/halla/sbs/efuchey/gemdigi_dig_ft_nobkgd_20181101_00 4
# python2 makejobs_digana.py digitized ft 1 -1 /volatile/halla/sbs/efuchey/gemdigi_dig_ft_bkgd_20181101_21 10













