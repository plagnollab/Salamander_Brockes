#fadd columns to annovar output. takes directory containing only summary and extra files
import os
import sys
import csv
import cPickle as pickle

############

refseq_annotation = pickle.load(open("data/blast_dict.p"))

########


###### This is the directory containing the RAW csv tables. - an example file is in there called "Contig_normalised_counts.csv" - contig IDs must be the first column
directory = "data/diff_expression"

#'''
blast_master_dict = {}


##### all the blast annotations calculated are stored in this directory
blst_directory = "blast_output/"
#"""
for file_ in os.listdir(blst_directory):

    print file_
    
    idz = file_.split("_")[1]
    
    print idz
    
    blast_master_dict[idz] = {}
    read_csv = csv.reader(open(os.path.join(blst_directory,file_)), delimiter="\t")
    
    blas_deets = {}
    
    for row in read_csv:
        
        if row[0].split("_")[-1] == "seq1":
            
            #print row[0]
            #print 
            
            #print row[1].split("|")[3],row[4],row[10]
            #raw_input()
            
            if idz == "human":
                blast_master_dict[idz][row[0].split("_s")[0]] = [row[1].split("|")[3],row[4],row[10]]
            else:
                blast_master_dict[idz][row[0].split("_s")[0]] = [row[1],row[4],row[10]]
    
    
#'''


annotation_master_dict = {}


################### all annotations from biomart are stored this this directory for human, zebrafish, xenopus
annotation_directory = "annotations_biomart"

for file_ in os.listdir(annotation_directory):
    
    print file_
        
    idz = file_.split("_")[0]
    
    print idz   
    
    #if idz == "human":
    #    continue

    read_csv = csv.reader(open(os.path.join(annotation_directory,file_)), delimiter="\t")
    
    annotation_master_dict[idz] = {}
    for row in read_csv:
       
        #print row
       
        #print row
        
        
        
        if row[0] == "" and row[1] != "":
            
            annotation_master_dict[idz][row[1]] = row[2:]
            
            #print row[1]
            
            
        elif row[0] != "":
            
            
            #print row[0]
            annotation_master_dict[idz][row[0]] = row[1:]




################# these lines are for debugging, you can directly access the hash table saving time, this is only necessary when testing though.

#pickle.dump(blast_master_dict,open("blast_master_dict","w"))

#pickle.dump(annotation_master_dict,open("annotation_master_dict","w"))
#"""

#blast_master_dict = pickle.load(open("blast_master_dict"))
#annotation_master_dict = pickle.load(open("annotation_master_dict"))
#raw_input("asdf")

################



for spec_ in blast_master_dict:  
    
    print len(blast_master_dict[spec_])
    
    
# For each raw count file in the input directory:
### run through each annotation set and include the annotation if present by appending to a single array and then writing the array.
for file_ in os.listdir(directory):
    
    if file_.startswith("."):
        continue
    
    if file_.endswith(".csv") and file_.find("annot") == -1:

        write_csv = csv.writer(open(os.path.join(directory,file_+"_annotated.csv"),"w"), delimiter=',',quoting=csv.QUOTE_ALL)
        
        if file_.find("counts") > -1:
            write_csv.writerow(["Contig ID",'s_310D0_D0', 's_312D0_D0', 's_32D0_D0', 's_33D0_D0', 's_34D0_D0', 's_310D15_D15', 's_312D15_D15', 's_32D15_D15', 's_33D15_D15', 's_34D15_D15', 's_310D150_D150', 's_312D150_D150', 's_32D150_D150', 's_33D150_D150', 's_34D150_D150']+["Description Refseq","gi ID",'RefSeq Protein ID',"Percent match","Evalue",'Associated Gene Name', 'Description', 'Ensembl Protein ID']+['Xenopus Associated Gene Name', 'Xenopus Description', 'Xenopus Ensembl Gene ID', 'Xenopus Percent match', 'Xenopus Evalue']+['Zebrafish Associated Gene Name', 'Zebrafish Description', 'Zebrafish Ensembl Protein ID', 'Zebrafish Percent match', 'Zebrafish Evalue'])
        else:        
            write_csv.writerow(["id","Contig name","baseMean","baseMeanA","baseMeanB","foldChange","log2FoldChange","pval","padj","Description Refseq","gi ID",'RefSeq Protein ID',"Percent match","Evalue",'Associated Gene Name', 'Description', 'Ensembl Protein ID']+['Xenopus Associated Gene Name', 'Xenopus Description', 'Xenopus Ensembl Gene ID', 'Xenopus Percent match', 'Xenopus Evalue']+['Zebrafish Associated Gene Name', 'Zebrafish Description', 'Zebrafish Ensembl Protein ID', 'Zebrafish Percent match', 'Zebrafish Evalue'])
        
        
        
        read_csv = csv.reader(open(os.path.join(directory,file_),"rU"), delimiter=",") 
        
        print file_
    
        for row in read_csv:
            
            print row
            
            #print row[0]
            
            if row[6] == "":
                continue
            
            id_ = row[6]
            
            print id_
            
            tmp_annotation = []
            #### iterate through each specie
            
            for spec_ in ["human","xenopus","zebrafish"]:  
                
                
                if spec_ == "human":
                    
                    blas_deets = blast_master_dict[spec_]
                    
                    biomart_dict = annotation_master_dict[spec_]
                    
                    
                    #print blas_deets.keys()

                    if id_ in blas_deets:
                        #print blas_deets[id_]
                        
                        refseq_prot_id = blas_deets[id_][0].split(".")[0]
                        
                        bio_row = ["","",""]
                        if refseq_prot_id in biomart_dict:
                            #print biomart_dict[refseq_prot_id]
                            bio_row = biomart_dict[refseq_prot_id][1:]
                            
                            #print len(bio_row)
                        #else:
                            #print refseq_prot_id
                            
                            #print "Not included"
                            #raw_input()
                            
                        if len(bio_row) != 3:
                            
                            bio_row = ["NA"]+bio_row
                            #raw_input("halt")
                            
                        #print refseq_annotation[refseq_prot_id]
                        # print refseq_prot_id
                    
                        #print len(refseq_annotation[refseq_prot_id]+blas_deets[id_]+bio_row)
                        #print refseq_annotation[refseq_prot_id]+blas_deets[id_]+bio_row
                    
                        #print bio_row
                        #raw_input()
                        tmp_annotation += refseq_annotation[refseq_prot_id]+blas_deets[id_]+bio_row
                        
                        #print len(refseq_annotation[refseq_prot_id]+blas_deets[id_]+bio_row)
                        #raw_input()

                    
                
                
                
                    else:
                        tmp_annotation += [""]*8
                        
                else:
                    
                    blas_deets3 = blast_master_dict[spec_]
                    
                    biomart_dict = annotation_master_dict[spec_]

                    if id_ in blas_deets3:
                        #print blas_deets3[id_]
                        #print id_
                        #raw_input("idee>")
                        
                        refseq_prot_id_2 = blas_deets3[id_][0]
                        
                        bio_row = ["","",""]
                        if refseq_prot_id_2 in biomart_dict:
                            #print biomart_dict[refseq_prot_id]
                            bio_row = biomart_dict[refseq_prot_id_2]
                            
                            #print len(bio_row)
                        else:
                            print "refseq_prot_id_2",refseq_prot_id_2
                            
                            #print "Not included"
                            #raw_input()
                            
                        #print refseq_annotation[refseq_prot_id]
                        # print refseq_prot_id
                        tmp_annotation += bio_row+blas_deets3[id_]
                        
                        #print len(bio_row+blas_deets3[id_])
                        #print "non spec"
                        #raw_input()
                    else:
                        
                        #print "No id in here, need to add empty space"
                        
                        tmp_annotation += [""]*5
                        
                        
                
            #print row+tmp_annotation
            #raw_input()
            write_csv.writerow(row+tmp_annotation)
            #raw_input()
    
    
    
    #["id",,"baseMean","baseMeanA","baseMeanB","foldChange","log2FoldChange","pval","padj",'RefSeq Protein ID [e.g. NP_001005353]', 'Associated Gene Name', 'Description', 'Ensembl Gene ID']
