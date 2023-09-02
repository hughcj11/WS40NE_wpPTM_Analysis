#Author: Chelsea Hughes
import csv
import re
base_path = "/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Proteomics Library"
#Below creates a dictionary within the code of unimod IDs and their biological relevance for later use.
with open("/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Testing code.csv") as csvfile:
    cellreader = csv.reader(csvfile, delimiter=',')
    dictionary={}
    for row in cellreader:
        key= row[0] + row[1]
        value=[row[2], row[3]]
        dictionary[key]=value

#The below code is opening the csv data from Skyline and assigning variables to each row

#Fix file name below


#with open(base_path+"EmbryohPTMs_Unimod.csv") as csvfile:
    #For a new csv file, change the above parenthesis to reflect your base path and specific document
    cellreader = csv.reader(csvfile, delimiter=',')
    abundance={}
    aminomap={}
    replicatename=[]
    count=0
    #For all rows, we collect "answers" in a list. These are the answers that we are getting as we go row by row
    answers=[]
    for row in cellreader:
        count=count+1
        if count==1:
            replicatename=row[8:]
        if count<=2:
            continue
        ua=row[3]
        #For each row, we will get an answer
        answer=[row[1],row[2],row[5],row[6],ua]
        #Below defines which row we need to begin reading to get several "answers", including the unimod and its position/residue
        pep_seq=row[3]
        #The below code removes any mass shifts that are not unimods, which appear in the Skyline data in brackets [].
        pattern = '\[.*\]'
        result = re.search(pattern, pep_seq)
        if result!=None:
            pep_seq=pep_seq.replace(result.group(),"")
        last_amino=""
        #unimod refers to the type of modifcation 
        #mod_amino refers to the residue being modified 
        #mod_pos refers to the position in that seqeunce being modified
        mod_amino=""
        unimod=""
        mod_pos=-1
        should_count=True
        replicatesclean=[]
        replicates=row[8:]
        for rep in replicates:
            if rep=="#N/A":
                rep=0
                replicatesclean.append(rep)
            else:
                replicatesclean.append(rep)
        replicates=replicatesclean
        for amino in pep_seq: 
            if amino==")":
                pep_mod_pos=mod_pos+int(row[5])
                answer.append(pep_mod_pos)
                answer.append(mod_amino)
                answer.append(unimod)
                key=row[1]+str(pep_mod_pos)+mod_amino+unimod
                
                if key in abundance.keys():
                    repsum=[] 
                    ziplist=zip(abundance[key],map(float,replicates))
                    for groupedreps in ziplist:
                        repsum.append(sum(groupedreps))  
                    abundance[key]=repsum
                else:
                    abundance[key]=list(map(float,replicates))
                #After the answer is recorded, the below text "resets" the code to count again
                mod_amino=""
                unimod=""
                should_count=True
                mod_pos = mod_pos -1
            else:
                if unimod!="":
                    unimod=unimod+amino
                if amino=="(":
                    should_count=False
                    mod_amino=last_amino
                if last_amino==":":
                    unimod=amino
            last_amino=amino
            if should_count:
                mod_pos=mod_pos+1
                #The below code is making a dictionary for each histon+amino acid+position that will have an associated replicates value. 
                #/this dictionary and its replicates will be referenced later to determine total possible abundance of hPTMs.       
                if amino!=")":
                    mapmap=row[1]+str(mod_pos+int(row[5]))+amino
                    if mapmap in aminomap.keys():
                        mapvalue=[] 
                        ziplist=zip(aminomap[mapmap],map(float,replicates))
                        for groupedreps in ziplist:
                            mapvalue.append(sum(groupedreps))  
                        aminomap[mapmap]=mapvalue
                    else:
                        aminomap[mapmap]=list(map(float,replicates))       
        answers.append(answer)
#The code below is transferring the information above into the first intermediate csv file (IntermediatePTMSheet). This intermediate sheet pulls out all the modifications on each peptide 
with open(base_path+'IntermediatePTMSheet.csv', 'w') as csvfile:
    cellwriter = csv.writer(csvfile, delimiter=',')
    cellwriter.writerow(["","","", "","", "Modification 1", "", "", "Modification 2","", "", "Modification 3"])
    cellwriter.writerow(["Protein","Protein Description","Begin Pos", "End Pos","Unimod Acession", "Position", "Residue", "Unimod", "Position", "Residue", "Unimod", "Position", "Residue", "Unimod"])
    for row_answer in answers:
        cellwriter.writerow(row_answer)


##Fix code below to not refer to histones

#The code below is transferring the information above into the second csv file (IntermediatePTMSheet2) that will ultimately produce hPTM IDs
with open(base_path+'IntermediatePTMSheet2.csv', 'w') as csvfile:
    cellwriter = csv.writer(csvfile, delimiter=',')
    cellwriter.writerow(["Protein","Protein Description","Position","Residue","Unimod", "hPTM_ID"])
    library=[]
    libraryspace=[]
    uniquehistone=[]
    for row_answer in answers:
        found_histone=False
        histone_result=""
        full_histone = row_answer[1] #Ex. PREDICTED: histone H1 1/2 [Austrofundulus limnaeus]
        parts = full_histone.split(" ")
        for histone in parts:
            if histone.startswith("["):
                found_histone= False
            if histone=="histone":
                found_histone=True
            elif found_histone==True:
                if histone_result!="":
                    histone_result=histone_result+"_"+histone
                elif histone_result=="":
                    histone_result=histone
        uniquehistone.append(histone_result)
        all_mods=row_answer[5:]
        for i in range(len(all_mods)):
            if i%3==0:
                position=all_mods[i]
            if i%3==1:
                residue=all_mods[i]
            if i%3==2:
                unimod=all_mods[i]
               #uniqueunimod.append(unimod)
                concat=row_answer[0]+"."+histone_result+"."+str(position)+"."+residue+"."+unimod
                concatspace=row_answer[0]+"$"+histone_result+"$"+str(position)+"$"+residue+"$"+unimod
                library.append(concat)
                libraryspace.append(concatspace)
                cellwriter.writerow([row_answer[0],histone_result,position,residue,unimod,concat])




#The below document lists only unique hPTM IDs and determines if they are biologically relevant
with open(base_path+'ProteinPTMLibrary.csv', 'w') as csvfile:
    cellwriter = csv.writer(csvfile, delimiter=',')
    cellwriter.writerow(["hPTM_ID","Protein Accession","Protein Description","Position","Amino Acid","Amino Acid + Position","Unimod","PTM Description","Biological Relevance"])
    pa=[]
    pd=[]
    pos=[] #"Position"
    aa=[] #"Amino Acid"
    aap=[]
    um=[]
    uniqueunimod=[]
    biorellist=[]
    library=list(set(library))
    libraryspace=list(set(libraryspace))
    for parts in libraryspace:
        allparts=parts.split("$")
        pa.append(allparts[0])
        pd.append(allparts[1])
        pos.append(allparts[2])
        aa.append(allparts[3])
        aap.append(allparts[2]+allparts[3])
        um.append(allparts[4])
        uniqueunimod.append(allparts[4]+"+"+allparts[3]) 
    for countindex in range(len(libraryspace)):
        if aa[countindex]+um[countindex] in dictionary.keys():
            unimodname=dictionary[aa[countindex]+um[countindex]][0]
            biorelevance=dictionary[aa[countindex]+um[countindex]][1]
            if biorelevance=="Yes" or biorelevance=="Multiple":
                biorellist.append([libraryspace[countindex].replace("$","."),pa[countindex],pd[countindex],pos[countindex],aa[countindex],aap[countindex],um[countindex],unimodname, biorelevance]) 
        else:
            unimodname=""
            biorelevance=""    
        cellwriter.writerow([libraryspace[countindex].replace("$","."),pa[countindex],pd[countindex],pos[countindex],aa[countindex],aap[countindex],um[countindex],unimodname, biorelevance]) 
     



#The below document is a list of the unique proteins identified.
with open(base_path+'UniqueProteinLibrary.csv', 'w') as csvfile:
    cellwriter = csv.writer(csvfile, delimiter=',')
    pa=set(pa)
    for uniqueanswer in pa:
        cellwriter.writerow([uniqueanswer])  


#The below document is a list of the unique unimods+residue identified. This can be used to update the master list of hPTMS and their biological relevance
with open(base_path+'UniqueUnimodLibrary.csv', 'w') as csvfile:
    cellwriter = csv.writer(csvfile, delimiter=',')
    cellwriter.writerow(["unimod+residue"])
    uniqueunimod=set(uniqueunimod)
    for uniquemod in uniqueunimod:
        cellwriter.writerow([uniquemod])

#The below document is a list of the unique residues that were modified. 
with open(base_path+'UniqueResidueLibrary.csv', 'w') as csvfile:
    cellwriter = csv.writer(csvfile, delimiter=',')
    cellwriter.writerow(["Amino Acid"])
    aa=set(aa)
    for uniqueaa in aa:
        cellwriter.writerow([uniqueaa])  

#The below document is a version of the ProteinPTMLibrary with only the biologically relevant modifications included
with open(base_path+'BioRelevantProteinPTMLibrary.csv', 'w') as csvfile:
    cellwriter = csv.writer(csvfile, delimiter=',')
    cellwriter.writerow(["hPTM_ID","Protein Accession","Protein Description","Position","Amino Acid","Amino Acid + Position","Unimod","PTM Description","Biological Relevance"])
    for rowanswer in biorellist:
        cellwriter.writerow(rowanswer)  

#The below document is a version of the biologically relevant ProteinPTMLibrary with the relative abundance calculated for each PTM type
with open(base_path+'Replicate calculations.csv', 'w') as csvfile:
    cellwriter = csv.writer(csvfile, delimiter=',')
    cellwriter.writerow(["hPTM_ID","Protein Accession","Protein Description","Position","Amino Acid","Amino Acid + Position","Unimod","PTM Description","Biological Relevance"]+replicatename+[""]+replicatename)
    for rowanswer in biorellist:
        key=rowanswer[1]+rowanswer[3]+rowanswer[4]+rowanswer[6]
        mapmap=rowanswer[1]+rowanswer[3]+rowanswer[4]
        if key in abundance.keys():
            rowanswer.extend(abundance[key])
        else:
            for something in replicatename:
                rowanswer.append(0)
        rowanswer.append("")
        if mapmap in aminomap.keys():
            rowanswer.extend(aminomap[mapmap])
        cellwriter.writerow(rowanswer)  
#The below document is a version of the biologically relevant ProteinPTMLibrary with the relative abundance calculated for each PTM type across the total possible residues (global)
with open(base_path+'Replicate_calculations_Total_PTMs.csv', 'w') as csvfile:
    cellwriter = csv.writer(csvfile, delimiter=',')
    cellwriter.writerow(["Unimod"]+replicatename+[""]+replicatename)
    PTMtotalpos={}
    PTMtotal={}
    for rowanswer in biorellist:
        key=rowanswer[1]+rowanswer[3]+rowanswer[4]+rowanswer[6]
        mapmap=rowanswer[1]+rowanswer[3]+rowanswer[4]
        PTMkey=rowanswer[7]
        PTMsums=[]
        PTMtotalsums=[]
        if key in abundance.keys():
            if PTMkey in PTMtotal.keys():
                ziplist=zip(PTMtotal[PTMkey],abundance[key])
                for eachrow in ziplist:
                    PTMsums.append(sum(eachrow))
                PTMtotal[PTMkey]=PTMsums
            else: 
                PTMtotal[PTMkey]=abundance[key]
        if mapmap in aminomap.keys():
            if PTMkey in PTMtotalpos.keys():
                ziplist=zip(PTMtotalpos[PTMkey],aminomap[mapmap])
                for eachrow in ziplist:
                    PTMtotalsums.append(sum(eachrow))
                PTMtotalpos[PTMkey]=PTMtotalsums
            else: 
                PTMtotalpos[PTMkey]=aminomap[mapmap]
    for PTManswers in PTMtotal.keys():
        cellwriter.writerow([PTManswers]+PTMtotal[PTManswers]+[""]+PTMtotalpos[PTManswers])



       