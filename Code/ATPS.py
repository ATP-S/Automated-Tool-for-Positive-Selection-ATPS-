import importt
from ATPS_functions import *
import os
import os.path
import shutil
import pandas as pd
import glob
import gc
import mne
import numpy as np
#python3 runn.py -G TP53,RB1,BRCA1 -S Mus_musculus,Homo_sapiens,Rattus_norvegicus,canis_lupus,panthera_tigris,ovis_aries,orycteropus_afer,cervus_canadensis,prionailurus_bengalensis -I Rattus_norvegicus
############# handle interest error

n = len(sys.argv)
genes = ''
Species = ''
interest = ''
inp_path = ''
tax = ''
save = ''
model_csv = ''
csv_instance = ''
models = []
for i in range(1,n):
    if sys.argv[i] == "-G": ## required
        print(sys.argv[i])
        genes =  sys.argv[i+1].split(',')
    if sys.argv[i] == "-S": ## optional
        Species = sys.argv[i+1].split(',')
    if sys.argv[i] == "-O": ## optional required comment: create a folder and get the path of the folder using pwd command
        save = sys.argv[i+1]
    if sys.argv[i] == "-T": ## optional
        tax = sys.argv[i+1] 
    if sys.argv[i] == "-I": ## optional
        interest = sys.argv[i+1]
    if sys.argv[i] == "-IF": ##optional
        inp_path = sys.argv[i+1]
    if sys.argv[i] == "-A": ##optional
        align_type = sys.argv[i+1]

try:
    del_codeml_dir()
except:
    pass
path = path_dir()
deletion_files()
try:
    os.remove("Gene_Output.csv")
except:
    print("")
print("=========================================================")
print("downloading and setup Gblocks and Jmodeltest ........ pleas wait")
print("=========================================================")

download()

try:
    os.system("codeml")
except ImportError:
    os.system("sudo apt-get install -y paml")   

#try:
#    os.system("codeml")
#except ImportError:
#    print("Sorry you didn't install codeml on your device so you can install it by following the instruction in this link : http://abacus.gene.ucl.ac.uk/software/paml.html")
#    q = input("to quit press 1 :")
#    if q == "1" :
#         sys.exit(0)

try:
    os.remove("codeml.ctl")
except:
    print("")
try:
    os.remove("rst")
except:
    print("")
try:
    os.remove("rst1")
except:
    print("")
try:
    os.remove("rub")
except:
    print("")
try:
    os.remove("lnf")
except:
    print("")
try:
    os.remove("2NG.dN")
except:
    print("")
try:
    os.remove("2NG.dS")
except:
    print("")
try:
    os.remove("2NG.t")
except:
    print("")
try:
    os.remove("4fold.nuc")
except:
    print("")
try:
    os.remove("mlc")
except:
    print("")

counts = 0
inp_file = ''
fetch = 1
if inp_path:
    fetch = 0
    inp_files = os.listdir(inp_path)
    print(inp_files)
    genes = [gene.split('.')[0] for gene in inp_files]
#model = "Name,Model0, Model7,Model8,Model8a,Model2a,Model2,LRT(M7_vs_M8),LRT(M8a_vs_M8),LRT(M2a_vs_M2),p-v7vs8,p-v8avs8,p-v2avs2"
LF = glob.glob("*")
if "Study_Output.csv" not in LF:
    os.system("echo 'Name,Model0, Model7,Model8,Model8a,Model2a,Model2,LRT(M7_vs_M8),LRT(M8a_vs_M8),LRT(M2a_vs_M2),p-v7vs8,p-v8avs8,p-v2avs2' >> Study_Output.csv")

codeml_creating_file()

for g in genes:
    path = os.getcwd()
    #fetchingbyspecies("TP53",["Mus musculus","Homo sapiens","Rattus norvegicus"])
    deletion_files()
    print("=========================================================")
    print("sequence fetching........ pleas wait")
    print("=========================================================")
    if fetch == 0:
        inp_file = inp_files[counts]

    list_empty, interestt = fetchingbyspecies(g,Species, interest, fetch, inp_path, inp_file)
    if list_empty == False:
        continue
    interest = interestt
    #key_list = fetching(tax, g)
    print("=========================================================")
    print("sequence MSA has been started ........ pleas wait")
    print("=========================================================")
    file_path = "ProteinSequences.fasta"
    diff_aligners(file_path , align_type)
    os.system("python3 exit.py")
    if interest != '':
        reversedd(interest)
    else: reversedd(interest)
    print("=========================================================")
    print("sequence filtration........ pleas wait")
    print("=========================================================")
    Gblocks()
    convert_fst_phy()
    rem_spaces()
    jmodel()
    try:
        partition, freq, pinvar = parsing_jmodeltest()
    except:
        os.system("sudo apt install default-jdk")
        convert_fst_phy()
        rem_spaces()
        print("=========================================================")
        print("model building........ pleas wait")
        print("=========================================================")
        jmodel()
        partition, freq, pinvar = parsing_jmodeltest()
    if partition == '' or freq == '':
        partition, freq, pinvar = spare_parse()
    if pinvar == '':
        pinvar = 'e'
    print("=========================================================")
    print("bulding tree ........ pleas wait")
    print("=========================================================")
    phyml(partition, freq, pinvar)
    try:
        convert_to_newickTree()
    except:
        os.system("sudo apt-get install -y phyml")
        phyml(partition, freq, pinvar)
        convert_to_newickTree()
    parsing_treefile()
    try:
        del_codeml_dir()
    except:
        pass
    print("=========================================================")
    print("phylofit running and wingscore will excecute ........ pleas wait")
    print("=========================================================")
    #try:
    #    phast(g)
    #except:
    #    pass
    creat_codeml_dir()
    os.system("python3 exit.py")
    print("=========================================================")
    print("codeml has been started model 078 ........ pleas wait it may take some time")
    print("=========================================================")
    model078()
    try:
        shutil.move("rub", 'codeml078')
    except:
        os.system("sudo apt-get install -y paml")
        model078()
        shutil.move("rub", 'codeml078')
    shutil.move("rst1", 'codeml078')
    shutil.move("rst", 'codeml078')
    shutil.move("lnf", 'codeml078')
    shutil.move("codeml078_mlc.txt", 'codeml078')
    shutil.move("2NG.t", 'codeml078')
    shutil.move("2NG.dS", 'codeml078')
    shutil.move("2NG.dN", 'codeml078')
    shutil.move("4fold.nuc", 'codeml078')
    os.system("python3 exit.py")
    o = open("genes_without_BEB.txt", "w")
    try:
        BEB_list = BEB(path)
    except:
        o.writelines(g)
        o.writelines("/n")
    try:
        Positive_selection_sites(BEB_list ,interest,path)
    except:
        pass
    print("=========================================================")
    print("codeml has been started model 8a ........ pleas wait it may take some time")
    print("=========================================================")
    model8a()
    print(interest)
    shutil.move("rub", 'codeml8a')
    shutil.move("rst1", 'codeml8a')
    shutil.move("rst", 'codeml8a')
    shutil.move("lnf", 'codeml8a')
    shutil.move("codeml8a_mlc.txt", 'codeml8a')
    shutil.move("2NG.t", 'codeml8a')
    shutil.move("2NG.dS", 'codeml8a')
    shutil.move("2NG.dN", 'codeml8a')
    shutil.move("4fold.nuc", 'codeml8a')
    os.system("python3 exit.py")
    if interest != '':
        print(interest)
        hashing(interest)
        print("=========================================================")
        print("codeml has been started model 078 ........ pleas wait it may take some time")
        print("=========================================================")
        model2a()
        shutil.move("rub", 'codeml2a')
        shutil.move("rst1", 'codeml2a')
        shutil.move("rst", 'codeml2a')
        shutil.move("lnf", 'codeml2a')
        shutil.move("codeml2a_mlc.txt", 'codeml2a')
        shutil.move("2NG.t", 'codeml2a')
        shutil.move("2NG.dS", 'codeml2a')
        shutil.move("2NG.dN", 'codeml2a')
        shutil.move("4fold.nuc", 'codeml2a')
        os.system("python3 exit.py")
        print("=========================================================")
        print("codeml has been started model 078 ........ pleas wait it may take some time")
        print("=========================================================")
        model2()
        shutil.move("rub", 'codeml2')
        shutil.move("rst1", 'codeml2')
        shutil.move("rst", 'codeml2')
        shutil.move("lnf", 'codeml2')
        shutil.move("codeml2_mlc.txt", 'codeml2')
        shutil.move("2NG.t", 'codeml2')
        shutil.move("2NG.dS", 'codeml2')
        shutil.move("2NG.dN", 'codeml2')
        shutil.move("4fold.nuc", 'codeml2')
        os.system("python3 exit.py")
    gc.collect()
    csv_instance = ''
    model_csv = ''
    if interest != '':
        if counts > 0:
            model_instance, model_csv = codeml_output(1, g)
            for m in range(len(model_instance)):
                models[m].extend(model_instance[m])
        else:
            models, model_csv = codeml_output(1,g)
        os.system(f"echo {model_csv} >> Study_Output.csv")
    else:
        codeml_output(0,g)
    print(models)
    print("=========================================================")
    print("file saving ........")
    print("=========================================================")
    saving_(g)
    del_codeml_dir()
    counts += 1
    try:
        os.remove("Gene_Output.csv")
    except:
        pass
    try:
        del_paths = glob.glob(os.path.join(path + '/' +  g + '.gene','*.gene'))
        for del_path in del_paths:
            shutil.rmtree(del_path)
    except:
        pass
##
##
df = pd.read_csv(r'Study_Output.csv')
p78 = np.asfarray(df['p-v7vs8'])
a, x = mne.stats.bonferroni_correction(p78, alpha=0.05)
p88a = np.asfarray(df['p-v8avs8'])
b, y = mne.stats.bonferroni_correction(p88a, alpha=0.05)
p22a = np.asfarray(df['p-v2avs2'])
c, z = mne.stats.bonferroni_correction(p22a, alpha=0.05)
df['adj78'] = x
df['adj88a'] = y
df['adj22a'] = z
df.to_csv(r"Study_Output_adjPvlaue.csv")
number_of_fetched_species()
