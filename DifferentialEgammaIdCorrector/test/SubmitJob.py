import os
from os import system, environ

submitFile="""
universe              = vanilla
Executable            = prepare_job.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT

Transfer_Input_Files  = Setup.tar.gz, DifferentialEgammaIdCorrector_cfg.py

#x509userproxy         = x509up_u110926

+maxWallTime          = 300
RequestMemory         = 8000
RequestDisk           = 15000

"""
fileParts = [submitFile]
files = open("DoubleEG_16H_run.log","r")
#files = open("test_sample.txt","r")
count = 0
for ij in files:
    count += 1
    fileParts.append("error     = Condor_jobs/logs/MCjob16{}_$(Cluster)_$(Process).stderr\n".format(count))
    fileParts.append("Log       = Condor_jobs/logs/MCjob16{}_$(Cluster)_$(Process).log\n".format(count))
    fileParts.append("Output    = Condor_jobs/logs/MCjob16{}_$(Cluster)_$(Process).stdout\n".format(count))
    fileParts.append("Arguments ={} {}\n".format(count, ij.strip()))
    fileParts.append("Queue\n\n")

    fout = open("Condor_jobs/condor_sub_Data16H_{}.sub".format(count),"w")
    fout.write(''.join(fileParts))
    fout.close()
    fileParts.pop(-1)
    fileParts.pop(-1)
    system('condor_submit Condor_jobs/condor_sub_Data16H_%i.sub' % count)

files.close()

