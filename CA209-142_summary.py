# Works! for TP, TH, MSI, MutSig, manifest (from paired file of Sujaya) and from the original manifest
# succesfully processed CA209-142, 141 trials
# Note, this is developed only for trails data that was processed before SB processing

MasterDict={}
import os, tarfile 
import zipfile

home=(raw_input("    The path to the trial folder is "))

# home='/Users/kalinavn/Desktop/temp_data_BMS/input'
# os.chdir(home)

# first I am creating a folder that will host all the files that I unzip to read
# after the final table is generated i can delete this folder

if os.path.exists(home+'/Delete_myself')==False:    
    os.makedirs(home+'/Delete_myself')
os.chdir(home)
for filename in os.listdir('.'):
    if filename.endswith('.tar.gz'):
            tar = tarfile.open(filename, 'r')
            N1=filename.strip().split('.')
            L1=N1[0].split('.')        
            folder_name=str(L1[0])                                        
            tar.extractall(path='./Delete_myself/'+folder_name)
            tar.close() 
    if filename.endswith('.zip'):
            N1=filename.strip().split('.')
            L1=N1[0].split('.')        
            folder_name=str(L1[0])                                                                                
            zip_file = zipfile.ZipFile(filename, 'r')
            zip_file.extractall(path='./Delete_myself/'+folder_name)                       
            zip_file.close()

# setting the name of the trial                        
os.chdir(home+'/Delete_myself')
for f1 in os.listdir('.'):
    F1=f1.split('_')
    if F1[-1]=='MutSig':
        trial=F1[0]
        MasterDict[trial]={}                
    
###  here I need to change to the manifest and populate the general info correctly from Sujaya's combined_ss.txt file
os.chdir(home)
infile=open('combined_ss.txt', 'rU')
for L1 in infile:
    M1=L1.strip().split('\t')
    subjectID=str(M1[0])
    gender=str(M1[1])
    tumor=str(M1[2])
    normal=str(M1[5])
    for t1 in MasterDict:                               
        if tumor!='NA' and normal!='NA':
            MasterDict[t1][tumor]={'PN':subjectID,
                                    'normal':normal,
                                    'gender':gender, 
                                    'type':'TN'}
        if tumor!='NA' and normal=='NA':
            MasterDict[t1][tumor]={'PN':subjectID,
                                    'normal':normal,
                                    'gender':gender, 
                                    'type':'TO'}                            
                                    
infile.close()

### here I am colecting additional information from another manifest, which carry info for WES and RNA-seq. 
### I need to verify that the EA IDs are correct that mathc to the PN and populate the data for 
## ACCNUM (6), BARCODE (7), EA (10)+number(13) = EA ID, BATCHID (17), WELLID (18)

manifest={}
import csv

# uncomment this bellow:
manifest_file=(raw_input("    what is the name of the manifest file, including the full path?  "))
f=open(manifest_file, 'rU')

#f=open(home+'/BMS-Manifest-P-20170209-0007-CA209-142-EA-20170801-1.csv', 'rU')
infile=csv.reader(f)
c=0
for M1 in infile:
    c+=1    
    if c>1:                  
        EA_ID=str(M1[10]+M1[13])

        manifest[EA_ID]={'STUDYID':str(M1[1]), # trial to verify 
                                   'USUBJID':str(M1[2]), # to verify
                                   'BMSPROJECTID':str(M1[0]),
                                   'ACCNUM':str(M1[6]),
                                   'BARCODE':str(M1[7]),
                                   'BATCHID':str(M1[17]),
                                   'WELLID':str(M1[18])
                                 }
                                                                                  
f.close()
for p1 in MasterDict:
    for p2 in MasterDict[p1]:
        MasterDict[p1][p2]['BMSPROJECTID']='NA'
        MasterDict[p1][p2]['ACCNUM']='NA'
        MasterDict[p1][p2]['BARCODE']='NA'
        MasterDict[p1][p2]['BATCHID']='NA'
        MasterDict[p1][p2]['WELLID']='NA'
                            
for p1 in MasterDict:
    for p2 in manifest:
        if MasterDict[p1].has_key(p2):
            MasterDict[p1][p2]['BMSPROJECTID']=manifest[p2]['BMSPROJECTID']
            MasterDict[p1][p2]['ACCNUM']=manifest[p2]['ACCNUM']
            MasterDict[p1][p2]['BARCODE']=manifest[p2]['BARCODE']
            MasterDict[p1][p2]['BATCHID']=manifest[p2]['BATCHID']
            MasterDict[p1][p2]['WELLID']=manifest[p2]['WELLID']
            
            if p1!=manifest[p2]['STUDYID']:
                print 'incorrect trial info '
            if MasterDict[p1][p2]['PN']!=manifest[p2]['USUBJID']:
                print 'incorrect patient number info'
                     
        

#################################################### end of manifest loop  ##########################
                            
####################### Start  #######################################
########## Tumor purity (TP) and Tumor heterogenuity (TH)  ##########
                                                        
# This is adding new keys into the master dictionary for Tumor Purity (TP) and Tumor Heteregenuity (TH)
cluster_list=['1','2','3','4','5'] 
for k1 in MasterDict:
    for k2 in MasterDict[k1]:
        MasterDict[k1][k2]['TumorPurity-Strelka']='NA' 
        MasterDict[k1][k2]['TumorPurity-Tnsnv']='NA'
        MasterDict[k1][k2]['TumorPurity-VarDict']='NA' 
        MasterDict[k1][k2]['TH_strelka_cluster_count']='NA'
        MasterDict[k1][k2]['TH_Tnsnv_cluster_count']='NA'
        MasterDict[k1][k2]['TH_VarDict_cluster_count']='NA'
                                    
        for i1 in cluster_list:
            MasterDict[k1][k2]['TH_strelka_cluster'+i1+'_percent']='NA'
            MasterDict[k1][k2]['TH_strelka_cluster'+i1+'_number_of_mut']='NA' 
            MasterDict[k1][k2]['TH_Tnsnv_cluster'+i1+'_percent']='NA'
            MasterDict[k1][k2]['TH_Tnsnv_cluster'+i1+'_number_of_mut']='NA'
            MasterDict[k1][k2]['TH_VarDict_cluster'+i1+'_percent']='NA'
            MasterDict[k1][k2]['TH_VarDict_cluster'+i1+'_number_of_mut']='NA'
                                      
# here I am unpackaging all the folders and files for strelka, tnsnv or VarDict files, so that I can read each file
# note that the directory for each folder must be done because the names for strelka and Tnsnv are identical, so they need to be placed in the corresponding folders         
for k1 in MasterDict:   
    os.chdir(home+'/Delete_myself/'+k1+'-Tumor_heterogeneity_results')
# uncomment these below if using it for the first time, once the files are unzipped, but you are still working on the script, comment them out    
    if os.path.exists(home+'/Delete_myself/'+k1+'-Tumor_heterogeneity_results/Strelka')==False:
        dir_strelka=os.mkdir(home+'/Delete_myself/'+k1+'-Tumor_heterogeneity_results/Strelka') 
    if os.path.exists(home+'/Delete_myself/'+k1+'-Tumor_heterogeneity_results/Tnsnv')==False:
        dir_tnsnv=os.mkdir(home+'/Delete_myself/'+k1+'-Tumor_heterogeneity_results/Tnsnv')
    if os.path.exists(home+'/Delete_myself/'+k1+'-Tumor_heterogeneity_results/VarDict')==False:
        dir_VarDict=os.mkdir(home+'/Delete_myself/'+k1+'-Tumor_heterogeneity_results/VarDict')
    for filename in os.listdir('.'):   
        if filename.endswith('.tar'):
            tar = tarfile.open(filename, 'r')
            N1=filename.strip().split('.')
            L1=N1[0].split('_')        
            if len(L1)>2:  
                if L1[4]=='Strelka':                                       
                    tar.extractall(path='./Strelka')                    
                if L1[4]=='Tnsnv':
                    tar.extractall(path='./Tnsnv')                       
            else: 
                tar.extractall(path='./VarDict')
            tar.close()
            
##### same to TP_Strelka
strelka={}                        
for p1 in MasterDict:
    os.chdir(home+'/Delete_myself/'+p1+'-Tumor_heterogeneity_results/Strelka')
    strelka[p1]={}
    for filename in os.listdir('.'):
        c=0 
        if filename.endswith('tumor_purity.txt'):
            c+=1
            B1=filename.strip().split('.')
            name1=str(B1[0])
            infile=open(filename, 'rU') 
            for L1 in infile:
                M1=L1.strip('\n').split(' ')
                strelka[p1][name1]=M1[-1]
            infile.close()
for P1 in MasterDict:         
    for name1 in strelka[P1]:
        if MasterDict[P1].has_key(name1):
            MasterDict[P1][name1]['TumorPurity-Strelka']=strelka[P1][name1]
        else:        
            for t1 in  MasterDict[P1]:               
                if name1==MasterDict[P1][t1]['normal']:
                    MasterDict[P1][t1]['TumorPurity-Strelka']=strelka[P1][name1]
                                
##### same to TP_Tnsnv
Tnsnv={}                        
for p1 in MasterDict:
    os.chdir(home+'/Delete_myself/'+p1+'-Tumor_heterogeneity_results/Tnsnv')
    Tnsnv[p1]={}
    for filename in os.listdir('.'):
        c=0 
        if filename.endswith('tumor_purity.txt'):
            c+=1
            B1=filename.strip().split('.')
            name1=str(B1[0])
            infile=open(filename, 'rU') 
            for L1 in infile:
                M1=L1.strip('\n').split(' ')
                Tnsnv[p1][name1]=M1[-1]
            infile.close()
for P1 in MasterDict:         
    for name1 in Tnsnv[P1]:
        if MasterDict[P1].has_key(name1):
            MasterDict[P1][name1]['TumorPurity-Tnsnv']=Tnsnv[P1][name1]
        else:        
            for t1 in  MasterDict[P1]:               
                if name1==MasterDict[P1][t1]['normal']:
                    MasterDict[P1][t1]['TumorPurity-Tnsnv']=Tnsnv[P1][name1]

##### same for TP_VarDict
VarDict={}                        
for p1 in MasterDict:
    os.chdir(home+'/Delete_myself/'+p1+'-Tumor_heterogeneity_results/VarDict')
    VarDict[p1]={}
    for filename in os.listdir('.'):
        c=0 
        if filename.endswith('tumor_purity.txt'):
            c+=1
            B1=filename.strip().split('.')
            name1=str(B1[0])
            infile=open(filename, 'rU') 
            for L1 in infile:
                M1=L1.strip('\n').split(' ')
                VarDict[p1][name1]=M1[-1]
            infile.close()
for P1 in MasterDict:         
    for name1 in VarDict[P1]:
        if MasterDict[P1].has_key(name1):
            MasterDict[P1][name1]['TumorPurity-VarDict']=VarDict[P1][name1]
        else:        
            for t1 in  MasterDict[P1]:               
                if name1==MasterDict[P1][t1]['normal']:
                    MasterDict[P1][t1]['TumorPurity-VarDict']=VarDict[P1][name1]  

################ end of Tumor Purity  ###########################                                      

######## To collect TH data  ##########

# collect the total number of clones 
# the percent for the top 5 clones 
# the number of mutations for the top 5 clones
# this is from the same file: *tumor_purity.txt

#### TH for Strelka
strelka_cluster={}
for P1 in MasterDict:
    os.chdir(home+'/Delete_myself/'+P1+'-Tumor_heterogeneity_results/Strelka')
    strelka_cluster[P1]={}
    for filename in os.listdir('.'):             
        if filename.endswith('tumor_heterogeneity_report.txt'):
            B1=filename.strip().split('.')
            name1=str(B1[0])
            strelka_cluster[P1][name1]={'count':0}
            for i1 in cluster_list:
                strelka_cluster[P1][name1]['cluster'+i1+'_percent']='NA'
                strelka_cluster[P1][name1]['cluster'+i1+'_number_of_mut']='NA' 
                                    
            infile=open(filename, 'rU') 
            b=0
            for L1 in infile:
                b+=1
                if L1.startswith('cluster'):
                    M1=L1.strip('\n').split('\t')
                    strelka_cluster[P1][name1]['count']+=1
                    for i1 in cluster_list:
                        if str(M1[0])=='cluster'+i1:
                            strelka_cluster[P1][name1]['cluster'+i1+'_percent']=str(M1[1])
                            strelka_cluster[P1][name1]['cluster'+i1+'_number_of_mut']=str(M1[2])                
            infile.close()
for P1 in MasterDict:
    for name1 in strelka_cluster[P1]:     
        if MasterDict[P1].has_key(name1):
            MasterDict[P1][name1]['TH_strelka_cluster_count']=strelka_cluster[P1][name1]['count']
            for i1 in cluster_list:
                MasterDict[P1][name1]['TH_strelka_cluster'+i1+'_percent']=strelka_cluster[P1][name1]['cluster'+i1+'_percent']
                MasterDict[P1][name1]['TH_strelka_cluster'+i1+'_number_of_mut']=strelka_cluster[P1][name1]['cluster'+i1+'_number_of_mut']
        else:
            for t2 in MasterDict[P1]:    
                if name1==MasterDict[P1][t2]['normal']:
                    MasterDict[P1][t2]['TH_strelka_cluster_count']=strelka_cluster[P1][name1]['count']
                    for i1 in cluster_list:
                        MasterDict[P1][t2]['TH_strelka_cluster'+i1+'_percent']=strelka_cluster[P1][name1]['cluster'+i1+'_percent']
                        MasterDict[P1][t2]['TH_strelka_cluster'+i1+'_number_of_mut']=strelka_cluster[P1][name1]['cluster'+i1+'_number_of_mut']


#### TH for Tnsnv
Tnsnv_cluster={}
for P1 in MasterDict:
    os.chdir(home+'/Delete_myself/'+P1+'-Tumor_heterogeneity_results/Tnsnv')
    Tnsnv_cluster[P1]={}
    for filename in os.listdir('.'):             
        if filename.endswith('tumor_heterogeneity_report.txt'):
            B1=filename.strip().split('.')
            name1=str(B1[0])
            Tnsnv_cluster[P1][name1]={'count':0}
            for i1 in cluster_list:
                Tnsnv_cluster[P1][name1]['cluster'+i1+'_percent']='NA'
                Tnsnv_cluster[P1][name1]['cluster'+i1+'_number_of_mut']='NA' 
                                    
         
            infile=open(filename, 'rU') 
            b=0
            for L1 in infile:
                b+=1
                if L1.startswith('cluster'):
                    M1=L1.strip('\n').split('\t')
                    Tnsnv_cluster[P1][name1]['count']+=1
                    for i1 in cluster_list:
                        if str(M1[0])=='cluster'+i1:
                            Tnsnv_cluster[P1][name1]['cluster'+i1+'_percent']=str(M1[1])
                            Tnsnv_cluster[P1][name1]['cluster'+i1+'_number_of_mut']=str(M1[2])                
            infile.close()
for P1 in MasterDict:
    for name1 in Tnsnv_cluster[P1]:     
        if MasterDict[P1].has_key(name1):
            MasterDict[P1][name1]['TH_Tnsnv_cluster_count']=Tnsnv_cluster[P1][name1]['count']
            for i1 in cluster_list:
                MasterDict[P1][name1]['TH_Tnsnv_cluster'+i1+'_percent']=Tnsnv_cluster[P1][name1]['cluster'+i1+'_percent']
                MasterDict[P1][name1]['TH_Tnsnv_cluster'+i1+'_number_of_mut']=Tnsnv_cluster[P1][name1]['cluster'+i1+'_number_of_mut']
        else:
            for t2 in MasterDict[P1]:    
                if name1==MasterDict[P1][t2]['normal']:
                    MasterDict[P1][t2]['TH_Tnsnv_cluster_count']=Tnsnv_cluster[P1][name1]['count']
                    for i1 in cluster_list:
                        MasterDict[P1][t2]['TH_Tnsnv_cluster'+i1+'_percent']=Tnsnv_cluster[P1][name1]['cluster'+i1+'_percent']
                        MasterDict[P1][t2]['TH_Tnsnv_cluster'+i1+'_number_of_mut']=Tnsnv_cluster[P1][name1]['cluster'+i1+'_number_of_mut']


#### TH for VarDict
VarDict_cluster={}
for P1 in MasterDict:
    os.chdir(home+'/Delete_myself/'+P1+'-Tumor_heterogeneity_results/VarDict')
    VarDict_cluster[P1]={}
    for filename in os.listdir('.'):             
        if filename.endswith('tumor_heterogeneity_report.txt'):
            B1=filename.strip().split('.')
            name1=str(B1[0])
            VarDict_cluster[P1][name1]={'count':0}
            for i1 in cluster_list:
                VarDict_cluster[P1][name1]['cluster'+i1+'_percent']='NA'
                VarDict_cluster[P1][name1]['cluster'+i1+'_number_of_mut']='NA' 
                                    
            infile=open(filename, 'rU') 
            b=0
            for L1 in infile:
                b+=1
                if L1.startswith('cluster'):
                    M1=L1.strip('\n').split('\t')
                    VarDict_cluster[P1][name1]['count']+=1
                    for i1 in cluster_list:
                        if str(M1[0])=='cluster'+i1:
                            VarDict_cluster[P1][name1]['cluster'+i1+'_percent']=str(M1[1])
                            VarDict_cluster[P1][name1]['cluster'+i1+'_number_of_mut']=str(M1[2])                
            infile.close()
for P1 in MasterDict:
    for name1 in VarDict_cluster[P1]:     
        if MasterDict[P1].has_key(name1):
            MasterDict[P1][name1]['TH_VarDict_cluster_count']=VarDict_cluster[P1][name1]['count']
            for i1 in cluster_list:
                MasterDict[P1][name1]['TH_VarDict_cluster'+i1+'_percent']=VarDict_cluster[P1][name1]['cluster'+i1+'_percent']
                MasterDict[P1][name1]['TH_VarDict_cluster'+i1+'_number_of_mut']=VarDict_cluster[P1][name1]['cluster'+i1+'_number_of_mut']
        else:
            for t2 in MasterDict[P1]:    
                if name1==MasterDict[P1][t2]['normal']:
                    MasterDict[P1][t2]['TH_VarDict_cluster_count']=VarDict_cluster[P1][name1]['count']
                    for i1 in cluster_list:
                        MasterDict[P1][t2]['TH_VarDict_cluster'+i1+'_percent']=VarDict_cluster[P1][name1]['cluster'+i1+'_percent']
                        MasterDict[P1][t2]['TH_VarDict_cluster'+i1+'_number_of_mut']=VarDict_cluster[P1][name1]['cluster'+i1+'_number_of_mut']

################### End of TH data ####################################################


####################### Start of MicroSatilte Instability (MSI) data #######################################

### first need to unzip all the MSI files, since each file has a unique name, all can be unziped in the same dir

for a1 in MasterDict:    
    for b1 in MasterDict[a1]:
        MasterDict[a1][b1]['MSI_Total_Number_of_Sites']='NA'
        MasterDict[a1][b1]['MSI_Number_of_Somatic_Sites']='NA'
        MasterDict[a1][b1]['MSI_%']='NA'
        
MSI={}        
for k1 in MasterDict:   
    os.chdir(home+'/Delete_myself/BMS_'+k1+'_MSI')
    MSI[k1]={}
    for filename in os.listdir('.'):  ## note that this loop I did twice. For some reason if I run it once, I cannot read the input files. Untar twice fixed the problem
        if filename.endswith('.tar.gz'):
            tar = tarfile.open(filename, 'r')
            tar.extractall(path='.')                       
            tar.close()
    for filename in os.listdir('.'):   
        if filename.endswith('.tar.gz'):
            tar = tarfile.open(filename, 'r')
            tar.extractall(path='.')                       
            tar.close()        
        if filename.endswith('msi'):
            L1=filename.strip().split('_')
            name1=L1[0]
            MSI[k1][name1]={}
            infile=open(filename, 'rU')
            
            c=0
            for N1 in infile:
                c+=1
                if c==1:
                    M1=N1.strip().split()
                    MSI[k1][name1]['MSI_Total_Number_of_Sites']='NA'
                    MSI[k1][name1]['MSI_Number_of_Somatic_Sites']='NA'
                    MSI[k1][name1]['MSI_%']='NA'
                if c==2:
                    M1=N1.strip().split()
                    MSI[k1][name1]['MSI_Total_Number_of_Sites']=M1[0]
                    MSI[k1][name1]['MSI_Number_of_Somatic_Sites']=M1[1]
                    MSI[k1][name1]['MSI_%']=M1[2]                 
            infile.close()
for P1 in MasterDict:
    for name1 in MSI[P1]:
        if MasterDict[P1].has_key(name1):
            MasterDict[P1][name1]['MSI_Total_Number_of_Sites']=MSI[P1][name1]['MSI_Total_Number_of_Sites']
            MasterDict[P1][name1]['MSI_Number_of_Somatic_Sites']=MSI[P1][name1]['MSI_Number_of_Somatic_Sites']
            MasterDict[P1][name1]['MSI_%']=MSI[P1][name1]['MSI_%']
                
        else:        
            for t1 in  MasterDict[P1]:               
                if name1==MasterDict[P1][t1]['normal']:
                    MasterDict[P1][t1]['MSI_Total_Number_of_Sites']=MSI[P1][name1]['MSI_Total_Number_of_Sites']
                    MasterDict[P1][t1]['MSI_Number_of_Somatic_Sites']=MSI[P1][name1]['MSI_Number_of_Somatic_Sites']
                    MasterDict[P1][t1]['MSI_%']=MSI[P1][name1]['MSI_%'] 

##################  End of MSI data  #########################################
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               

######################### Start Mutation Signature (MutSig)  ############################

### First need to unzip file, each has a unique name, so I can unzip in the same directory

for a1 in MasterDict:    
    for b1 in MasterDict[a1]:
        MasterDict[a1][b1]['top_MutSig']='NA'
        for j1 in range (1,31):
            MasterDict[a1][b1]['MutSig'+str(j1)]='NA'
        
MutSig={}        
for P1 in MasterDict:   
    os.chdir(home+'/Delete_myself/'+P1+'_MutSig')    
    MutSig[P1]={}
    for filename in os.listdir('.'):
        if filename.endswith('.zip'):
            zip_file = zipfile.ZipFile(filename, 'r')
            # print os.getcwd()
            # print filename
            zip_file.extractall(path='.')                       
            zip_file.close() 
        os.chdir(home+'/Delete_myself/'+P1+'_MutSig')
        for filename in os.listdir('.'):
            if filename.endswith('weights.txt'):
                N1=filename.strip().split('_')
                name1=str(N1[1])
                MutSig[P1][name1]={}
                infile=open(filename, 'rU')
                MutSig[P1][name1]['top_MutSig']='NA'
                for m1 in range(1, 31):
                    MutSig[P1][name1]['Signature'+str(m1)]='NA'           
                c=0
                for L1 in infile:
                    c+=1
                    if c==2:
                        M2=L1.strip().split()
                        for s1 in range (1, 31):
                            MutSig[P1][name1]['Signature'+str(s1)]=M2[s1]
                            if M2[s1]==max(M2[1:]):
                                MutSig[P1][name1]['top_MutSig']='Signature'+str(s1)               
                infile.close()
            
                                            
for P1 in MasterDict:
    for name1 in MutSig[P1]:
        if MasterDict[P1].has_key(name1):
            MasterDict[P1][name1]['top_MutSig']=MutSig[P1][name1]['top_MutSig']
            for s1 in range(1,31):
                MasterDict[P1][name1]['MutSig'+str(s1)]=MutSig[P1][name1]['Signature'+str(s1)]
        else:        
            for t1 in  MasterDict[P1]:               
                if name1==MasterDict[P1][t1]['normal']:
                    MasterDict[P1][t1]['top_MutSig']=MutSig[P1][name1]['top_MutSig']
                    for s1 in range(1,31):
                        MasterDict[P1][t1]['MutSig'+str(s1)]=MutSig[P1][name1]['Signature'+str(s1)]


##################### End of MutSig ########################################################

##################### Adding HLA  #########################################################

os.chdir(home+'/HLA/tumor')
HLA={}        

for filename in os.listdir('.'):
    if filename.endswith('.result.tsv'):
        for k1 in MasterDict:
            for k2 in MasterDict[k1]:
                MasterDict[k1][k2]['tumor_A1']='NA'
                MasterDict[k1][k2]['tumor_A2']='NA'
                MasterDict[k1][k2]['tumor_B1']='NA'
                MasterDict[k1][k2]['tumor_B2']='NA'
                MasterDict[k1][k2]['tumor_C1']='NA'
                MasterDict[k1][k2]['tumor_C2']='NA'
                
                MasterDict[k1][k2]['normal_A1']='NA'
                MasterDict[k1][k2]['normal_A2']='NA'
                MasterDict[k1][k2]['normal_B1']='NA'
                MasterDict[k1][k2]['normal_B2']='NA'
                MasterDict[k1][k2]['normal_C1']='NA'
                MasterDict[k1][k2]['normal_C2']='NA'
                
                if MasterDict[k1][k2]['PN'] in filename:
                    PN=MasterDict[k1][k2]['PN']
                    HLA[PN]={}
                    infile=open(filename, 'rU')
                    HLA[PN]['tumor_A1']='NA'
                    HLA[PN]['tumor_A2']='NA'
                    HLA[PN]['tumor_B1']='NA'
                    HLA[PN]['tumor_B2']='NA'
                    HLA[PN]['tumor_C1']='NA'
                    HLA[PN]['tumor_C2']='NA'
                    
                    HLA[PN]['normal_A1']='NA'
                    HLA[PN]['normal_A2']='NA'
                    HLA[PN]['normal_B1']='NA'
                    HLA[PN]['normal_B2']='NA'
                    HLA[PN]['normal_C1']='NA'
                    HLA[PN]['normal_C2']='NA'
                          
                    c=0
                    for L1 in infile:
                        c+=1
                        if c==2:
                            M2=L1.strip().split()
                            HLA[PN]['tumor_A1']=M2[1]
                            HLA[PN]['tumor_A2']=M2[2]
                            HLA[PN]['tumor_B1']=M2[3]
                            HLA[PN]['tumor_B2']=M2[4]
                            HLA[PN]['tumor_C1']=M2[5]
                            HLA[PN]['tumor_C2']=M2[6]
                                      
                    infile.close()
          
                                            
for k1 in MasterDict:
    for k2 in MasterDict[k1]:
        for PN in HLA:
            if PN==MasterDict[k1][k2]['PN']:
                MasterDict[k1][k2]['tumor_A1']=HLA[PN]['tumor_A1']
                MasterDict[k1][k2]['tumor_A2']=HLA[PN]['tumor_A2']
                MasterDict[k1][k2]['tumor_B1']=HLA[PN]['tumor_B1']
                MasterDict[k1][k2]['tumor_B2']=HLA[PN]['tumor_B2']
                MasterDict[k1][k2]['tumor_C1']=HLA[PN]['tumor_C1']
                MasterDict[k1][k2]['tumor_C2']=HLA[PN]['tumor_C2']

os.chdir(home+'/HLA/normal')
HLA={}        

for filename in os.listdir('.'):
    if filename.endswith('.result.tsv'):
        for k1 in MasterDict:
            for k2 in MasterDict[k1]:
                if MasterDict[k1][k2]['PN'] in filename:
                    PN=MasterDict[k1][k2]['PN']
                    HLA[PN]={}
                    infile=open(filename, 'rU')
                    c=0
                    for L1 in infile:
                        c+=1
                        if c==2:
                            M2=L1.strip().split()
                            HLA[PN]['normal_A1']=M2[1]
                            HLA[PN]['normal_A2']=M2[2]
                            HLA[PN]['normal_B1']=M2[3]
                            HLA[PN]['normal_B2']=M2[4]
                            HLA[PN]['normal_C1']=M2[5]
                            HLA[PN]['normal_C2']=M2[6]
                                      
                    infile.close()
          
                                            
for k1 in MasterDict:
    for k2 in MasterDict[k1]:
        for PN in HLA:
            if PN==MasterDict[k1][k2]['PN']:
                MasterDict[k1][k2]['normal_A1']=HLA[PN]['normal_A1']
                MasterDict[k1][k2]['normal_A2']=HLA[PN]['normal_A2']
                MasterDict[k1][k2]['normal_B1']=HLA[PN]['normal_B1']
                MasterDict[k1][k2]['normal_B2']=HLA[PN]['normal_B2']
                MasterDict[k1][k2]['normal_C1']=HLA[PN]['normal_C1']
                MasterDict[k1][k2]['normal_C2']=HLA[PN]['normal_C2']


###################### End of HLA  ########################################################

#################### Metrics info  ##################################################
os.chdir(home+'/Metrics')
Metrics={}
metrics_header=[]
for folder in os.listdir('.'):
    if folder.startswith('EA'):
        os.chdir(home+'/Metrics/'+folder)
        for filename in os.listdir('.'):
            M1=filename.strip().split('.')
            name1=M1[0]
            Metrics[name1]={} 
            infile=open(filename,'rU')
            c=0
            for L1 in infile:
                c+=1
                if c==7:
                    metrics_header=L1.strip('\n').split()
                if c==8:
                    Metrics[name1]=L1.strip('\n')
            infile.close()     
        os.chdir(home+'/Metrics')
for k1 in MasterDict:
    for k2 in MasterDict[k1]:
        MasterDict[k1][k2]['tumor_Metrics']='NA\t'*55+'NA'
        MasterDict[k1][k2]['normal_Metrics']='NA\t'*55+'NA'
        for name1 in Metrics:
            if MasterDict[k1].has_key(name1):
                MasterDict[k1][name1]['tumor_Metrics']=Metrics[name1]
            if name1==MasterDict[k1][k2]['normal']:
                MasterDict[k1][k2]['normal_Metrics']=Metrics[name1]          
                        
###################### End of Metrics data  ###################################

###################### Start of Neoiantigens ##################################
threshold=0.80 # cutoff for the neoantigens combined_prediction_score value
os.chdir(home+'/Delete_myself')
neoantigen={}
for folder in os.listdir('.'):
    if folder.endswith('Neoantigen'):
        os.chdir(home+'/Delete_myself/'+folder)
        for filename in os.listdir('.'):
            name1=filename.strip().split('.')
            neoantigen[name1[0]]={'total_number_of_neoantigens':0,
                                  'number_of_neoantigens_threshold_'+str(threshold):0}               
            infile=open(filename, 'rU')
            c=0
            for L1 in infile:
                c+=1
                if c>1:
                    neoantigen[name1[0]]['total_number_of_neoantigens']+=1
                    M1=L1.strip().split()
                    if float(M1[5])>= float(threshold):
                        neoantigen[name1[0]]['number_of_neoantigens_threshold_'+str(threshold)]+=1 
            infile.close()
for k1 in MasterDict:
    for k2 in MasterDict[k1]:
        MasterDict[k1][k2]['total_number_of_neoantigens']='NA'
        MasterDict[k1][k2]['number_of_neoantigens_threshold_'+str(threshold)]='NA'
        for k3 in neoantigen:
            if k3==MasterDict[k1][k2]['PN']:
                MasterDict[k1][k2]['total_number_of_neoantigens']=neoantigen[k3]['total_number_of_neoantigens']
                MasterDict[k1][k2]['number_of_neoantigens_threshold_'+str(threshold)]=neoantigen[k3]['number_of_neoantigens_threshold_'+str(threshold)]
                    
###################### End of Neoantigens  #####################################

################# writing the output! ##########################################

os.chdir(home) 
for p1 in MasterDict:
    outfile=open(p1+'_summary.txt','w')
outfile.write('BMSPROJECTID'+'\t')
outfile.write('trial_ID\t'+
              'PN\t'+
              'tumor_ID\t'+
              'normal_ID\t'+
              'gender\t'+
              'type\t'+
              
              'TP-strelka\t'+
              'TP-Tnsnv\t'+
              'TP-VarDict\t'+
              
              'TH_strelka_cluster_count\t'+
              'TH_Tnsnv_cluster_count\t'+
              'TH_VarDict_cluster_count')
for i1 in cluster_list:
    outfile.write('\t'+'TH_strelka_cluster'+i1+'_percent')
    outfile.write('\t'+'TH_strelka_cluster'+i1+'_number_of_mut')
for i1 in cluster_list:    
    outfile.write('\t'+'TH_Tnsnv_cluster'+i1+'_percent')
    outfile.write('\t'+'TH_Tnsnv_cluster'+i1+'_number_of_mut')
for i1 in cluster_list:    
    outfile.write('\t'+'TH_VarDict_cluster'+i1+'_percent')
    outfile.write('\t'+'TH_VarDict_cluster'+i1+'_number_of_mut')
    
outfile.write('\t'+'Top_MutationSignature')    
for h1 in range(1,31):
    outfile.write('\t'+'MutationSignature'+str(h1))   
    
outfile.write('\t'+ 'MSI_Total_Number_of_Sites'+
              '\t'+  'MSI_Number_of_Somatic_Sites'+
              '\t'+  'MSI_%')

outfile.write('\t'+'tumor_HLA-A1'+
              '\t'+'tumor_HLA-A2'+
              '\t'+'tumor_HLA-B1'+
              '\t'+'tumor_HLA-B2'+
              '\t'+'tumor_HLA-C1'+
              '\t'+'tumor_HLA-C2'+
              '\t'+'normal_HLA-A1'+
              '\t'+'normal_HLA-A2'+
              '\t'+'normal_HLA-B1'+
              '\t'+'normal_HLA-B2'+
              '\t'+'normal_HLA-C1'+
              '\t'+'normal_HLA-C2') 

outfile.write( '\t'+'Neoantigens_total_number'+
               '\t'+'Neoantigens_with_threshold_'+str(threshold))                                                                                  
               
outfile.write(
              '\t'+'manifest_ACCNUM'+
              '\t'+'manifest_BARCODE'+
              '\t'+'manifest_BATCHID'+
              '\t'+'manifest_WELLID')
              
for j1 in range(0, len(metrics_header)):
    outfile.write('\t'+'tumor_Metrics_'+metrics_header[j1])
for j1 in range(0, len(metrics_header)):
    outfile.write('\t'+'normal_Metrics_'+metrics_header[j1])                                         
             
outfile.write('\n')
              
              
for k1 in MasterDict:
    for k2 in MasterDict[k1]:
        outfile.write(str(MasterDict[k1][k2]['BMSPROJECTID'])+'\t')              
        outfile.write(str(k1)+'\t'+
                      str(MasterDict[k1][k2]['PN'])+'\t'+
                      str(k2)+'\t'+
                      str(MasterDict[k1][k2]['normal'])+'\t'+
                      str(MasterDict[k1][k2]['gender'])+'\t'+
                      str(MasterDict[k1][k2]['type'])+'\t'+
                      
                      str(MasterDict[k1][k2]['TumorPurity-Strelka'])+'\t'+
                      str(MasterDict[k1][k2]['TumorPurity-Tnsnv'])+'\t'+
                      str(MasterDict[k1][k2]['TumorPurity-VarDict'])+'\t'+
                      
                      str(MasterDict[k1][k2]['TH_strelka_cluster_count'])+'\t'+
                      str(MasterDict[k1][k2]['TH_Tnsnv_cluster_count'])+'\t'+
                      str(MasterDict[k1][k2]['TH_VarDict_cluster_count']))
        for i1 in cluster_list:
            outfile.write('\t%s'% str(MasterDict[k1][k2]['TH_strelka_cluster'+i1+'_percent']))
            outfile.write('\t%s'% str(MasterDict[k1][k2]['TH_strelka_cluster'+i1+'_number_of_mut']))
        for i1 in cluster_list:    
            outfile.write('\t%s'% str(MasterDict[k1][k2]['TH_Tnsnv_cluster'+i1+'_percent']))
            outfile.write('\t%s'% str(MasterDict[k1][k2]['TH_Tnsnv_cluster'+i1+'_number_of_mut']))
        for i1 in cluster_list:    
            outfile.write('\t%s'% str(MasterDict[k1][k2]['TH_VarDict_cluster'+i1+'_percent']))
            outfile.write('\t%s'% str(MasterDict[k1][k2]['TH_VarDict_cluster'+i1+'_number_of_mut']))
            
        outfile.write('\t'+str(MasterDict[k1][k2]['top_MutSig']))
        for h1 in range(1,31):
            outfile.write('\t%s'% str(MasterDict[k1][k2]['MutSig'+str(h1)]))            
            
        outfile.write('\t'+str(MasterDict[k1][k2]['MSI_Total_Number_of_Sites'])+
                      '\t'+str(MasterDict[k1][k2]['MSI_Number_of_Somatic_Sites'])+
                      '\t'+str(MasterDict[k1][k2]['MSI_%']))
                      
        outfile.write('\t'+str(MasterDict[k1][k2]['tumor_A1'])+
                      '\t'+str(MasterDict[k1][k2]['tumor_A2'])+
                      '\t'+str(MasterDict[k1][k2]['tumor_B1'])+
                      '\t'+str(MasterDict[k1][k2]['tumor_B2'])+
                      '\t'+str(MasterDict[k1][k2]['tumor_C1'])+
                      '\t'+str(MasterDict[k1][k2]['tumor_C2'])+
                      
                      '\t'+str(MasterDict[k1][k2]['normal_A1'])+
                      '\t'+str(MasterDict[k1][k2]['normal_A2'])+
                      '\t'+str(MasterDict[k1][k2]['normal_B1'])+
                      '\t'+str(MasterDict[k1][k2]['normal_B2'])+
                      '\t'+str(MasterDict[k1][k2]['normal_C1'])+
                      '\t'+str(MasterDict[k1][k2]['normal_C2']))
                      
        outfile.write('\t'+str(MasterDict[k1][k2]['total_number_of_neoantigens'])+
                      '\t'+str(MasterDict[k1][k2]['number_of_neoantigens_threshold_'+str(threshold)]))                                          
        
        outfile.write(
                      '\t'+str(MasterDict[k1][k2]['ACCNUM'])+
                      '\t'+str(MasterDict[k1][k2]['BARCODE'])+
                      '\t'+str(MasterDict[k1][k2]['BATCHID'])+
                      '\t'+str(MasterDict[k1][k2]['WELLID']))
        outfile.write('\t'+MasterDict[k1][k2]['tumor_Metrics']+
                      '\t'+MasterDict[k1][k2]['normal_Metrics'])             
                                           
        outfile.write('\n')   
                      
                             
                          
outfile.close()                      

             
############################# end of writing the output file, it is outside the Delet_muyself folder  ##############
    
########## Cleaning up!!! This is to delete the 'Delete_myself' folder and all the files in it!  ####################
rerun=raw_input("Check the output file! Is it good? Yes or No?")

os.chdir(home)

import shutil
if rerun=="Yes":
    shutil.rmtree(home+'/Delete_myself')
    print "Done!"
    
        
                          
                                                                              
############# This is the end of the programm !!! ###################                        
