# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:18:27 2020

@author: rmagni
"""
import itertools 
from itertools import repeat
import multiprocessing
import concurrent.futures
from pathlib import Path
import pandas as pd
import time
import argparse
import glob
import os

### START: RAFA MODIFICATION ###
import numpy as np
### END: RAFA MODIFICATION ###


def Enzyme(aa_cut,not_cut,isoleu):
    c=()
    

    for i in aa_cut:
    
        for a in  not_cut:
        
            c1=(i,i+",")
            nc1=(i+","+a,i+a)
            c=c+(c1,)+(nc1,)
   
    if isoleu !="NOCHANGE" :
       
        
       ci=tuple([(i,isoleu) for i in ["I","L","J"] if i!=isoleu])
       c=c+ci
        
    return c     



def find_offsets(haystack, needle):

    offs = -1
    while True:
        offs = haystack.find(needle, offs+1)
        if offs == -1:
            break
        else:
            yield offs



def prev_next_aa(i,a):
    
    if a>0 and a<len(i)-1:
        return (i[a-1][-1],i[a+1][0])
    
    elif a==0:
        return ("NtermP",i[a+1][0])
    
    elif a==len(i)-1: 
        
         return (i[a-1][-1],"CtermP")




def Digest(seq,nm,minpeptidelength,maxpeptidelength,aa_cut,not_cut,isoleu):
    
    c=Enzyme(aa_cut,not_cut,isoleu)     
    plen=len(seq)
    
    for i in c:
        seq=seq.replace(*i)
        if seq.endswith(','):
            seq = seq[:-1]
    
    if "," in seq[:-1]:

        pos=[offs-n for n,offs in enumerate(find_offsets(seq, ","))]
        pos.insert(0,0)
        pos.insert(plen,plen)
        seq=seq.split(",")
        while "" in seq: seq.remove("")  
        npaa=[prev_next_aa(seq,a) for a,i in enumerate(seq)]
        
        
        sq = [[] for i in range(nm+1)]
        ps= [[] for i in range(nm+1)]
        npe=[[] for i in range(nm+1)]
    
        for i,j in enumerate (sq):
        
            sq1=[seq[i:] for i in range(i+1)]
            ps1=[pos[i:] for i in range(i+1)]
            np1=[npaa[i:] for i in range(i+1)]
            sq[i]=["".join(i) for  i in list(zip(*sq1))]
            ps[i]=[str(i[0]+1)+"_"+str(i[1]) for i in zip(ps1[0][:len(ps1[-1][1:])],ps1[-1][1:])]
            npe[i]=[str(i[0][0])+"_"+str(i[1][1]) for i in zip(np1[0][:len(ps1[-1])],np1[-1])]


        sq=[y for y in sq if y]
        sq=[list(itertools.zip_longest(ps[i],npe[i],sq[i],[i],fillvalue=i)) for i,j in enumerate(sq)]
        sq=[item for sublist in sq for item in sublist if ((len(item[2])>=minpeptidelength) and (len(item[2])<=maxpeptidelength))]
        ps="_".join([str(len([x  for x in sq if x[3]==i])) for i,j in enumerate(ps) ])
        
        if not sq:
            
            sq=[["N/A_N/A","N/A_N/A","No valid peptides","N/A"]]
        
        return sq,str(plen),ps,str(len(sq))
    
    else:
        
        if ((len(seq)>=minpeptidelength) and (len(seq)<=maxpeptidelength)):
            
            sq=[["1"+str(len(seq)),"NtermP_CtermP",seq,0]]
        
            return sq,str(len(seq)),str(1),str(1)
    
        else:
        
            return [["N/A_N/A","N/A_N/A","No valid peptides","N/A"]],str(0),str(0),str(0)


def expand(sh):
    
    f=[[tup[0]+list(tup[1]) for tup in (itertools.zip_longest([[a,b,c,d,e]],f,fillvalue=[a,b,c,d,e]))] for a,b,c,d,e,f in sh] 
    f=[item for sublist in f for item in sublist]
    return f


def Fastareader(Fastapath,nm,decoy_tag,minpeptidelength,maxpeptidelength,aa_cut,not_cut,isoleu):
    
    fasta=open(Fastapath,"r+")
    header=list()
    headerabr=list()
    sequence=list()
    proteinlen=list()
    nomisspeptides=list()
    misspeptides=list()
    seq=""
    qlen=""
    nmp=""
    mp=""
    f=fasta.readlines()
    f.append(">")
    
    for i in f:
    	
        if i[0] == '>':
            if seq!="":

                seq,qlen,nmp,mp=Digest(seq,nm,minpeptidelength,maxpeptidelength,aa_cut,not_cut,isoleu)
                sequence.append(seq)
                proteinlen.append(qlen)
                nomisspeptides.append(nmp)
                misspeptides.append(mp)
                seq = '';
            
            else:
                
                 if len(header)>1:
                     
                    del header[-1]
                    del headerabr[-1]
        
            header.append(i.rstrip())
            
            if "|" in i:
                
                headerabr.append(i.rstrip().split("|")[1]) 
            else:
                
                headerabr.append("NoUniprotID") 
                
        else:
            seq+=i.rstrip()
    
    fasta.close()
    
    return (i for i in zip(headerabr,header,proteinlen,nomisspeptides,misspeptides,sequence) if decoy_tag not in i[1])


def joiner(Filepath,dffasta,PDO,PDS,PHR,pepcolumn):    
    dffile=pd.read_csv(Filepath,sep="\t",header=PHR)
    dffile=dffile.merge(dffasta, how='left', on=pepcolumn) 

    ### START: RAFA MODIFATION ###
    #start = time.time()
    dffile = getMostProbableProtein(dffile, pepcolumn)
    #print(f"Calculated most probable protein in {divmod(time.time()-start,60)[0]}m and {round(divmod(time.time()-start,60)[1],4)}s")
    ### END: RAFA MODIFICATION ###

    dffile.to_csv(os.path.splitext(Filepath)[0]+PDS+os.path.splitext(Filepath)[1],sep="\t",index=False)

    

### START: RAFA MODIFICATION ###

def getMostProbableProtein(dffile, pepcolumn):
    '''
    
    '''
    # Columns containing information from multiple protein separated by ";". Name of columns got from "col" variable in main
    semicolon_col_list = ["ProteinsLength","NpMiss","NpTotal","PeptidePos","PrevNextaa"]
    semicolon_col_protein_list = ["Proteins", "UniProtIDs"]

    semicolon_col_work_list = [i for i in semicolon_col_protein_list+semicolon_col_list if i in dffile.columns]
    semicolon_col_protein_work_list = [i for i in semicolon_col_protein_list if i in semicolon_col_work_list]

    if len(semicolon_col_protein_work_list) == 0:
        print(f"Error: Most probable protein cannot be calculated. At least one column is required: {semicolon_col_protein_list}")
        return dffile


    # Name of the column with protein information (it can be UniProtIDs or Proteins)
    sc_prot_col = semicolon_col_protein_work_list[0]

    # Get indexes of dffile
    df_index = dffile.index.to_list()

    # This array will contain the position of the most probable protein in each row
    df_index_result_pos = np.ones_like(df_index)*(-1)

    # Get boolean with non decoy PSMs (has no real protein assigned)
    non_decoy_bool = ~dffile[sc_prot_col].isnull()

    # Extract sc_prot_col from dffile. Generate List of Lists [ [p1, p2...], [p1, p2..] ... ]
    sc_prot_list = [tuple(i.split(" // ")) for i in dffile.loc[non_decoy_bool, sc_prot_col].to_list()]

    # Extract peptide sequence of each scan: pepcolumn of dffile
    p_seq_list = dffile.loc[non_decoy_bool, pepcolumn].to_list()

    # Get set of pairs (peptide sequence, [protein list]). Do not repeat peptide sequence!
    pseq_prot_pair_list = list(set(list(zip(p_seq_list, sc_prot_list))))
    
    # Get flat list with all proteins contributed by each peptide sequence to get number of peptides per protein
    protein_from_pseq = sorted([j for i in pseq_prot_pair_list for j in i[1]])
    protein2npep = {k : len(list(g)) for k, g in itertools.groupby(protein_from_pseq)}

    # Get flat list with all proteins contributed by each scan to get number of scans per protein
    protein_from_scan = sorted([j for i in sc_prot_list for j in i])
    protein2nscan = {k : len(list(g)) for k, g in itertools.groupby(protein_from_scan)}

    # Extract elements of sc_prot_list with more than one protein
    sc_prot_gt1_bool_arr = np.array([len(i) for i in sc_prot_list]) > 1
    sc_prot_gt1_list = [i for i,j in zip(sc_prot_list, sc_prot_gt1_bool_arr.tolist()) if j]

    
    # Resolve elements with one protein only: From full index, get non-decoy. And from them, get those with one protein
    # aux_arr is analogous to non-decoy elements. We first index those with one protein only.
    # aux_arr is subset of df_index_result_pos --> aux_arr = df_index_result_pos[non_decoy_bool]
    
    aux_arr = np.ones_like(np.arange(0, len(sc_prot_list)))*(-1)
    aux_arr[~sc_prot_gt1_bool_arr] = 0 # protein in position 0 for elements with one protein only

    # aux_arr2 is analogous to elements with more than one protein. We first index those with only one maximum
    # aux_arr2 is a subset of aux_arr --> aux_arr2 = aux_arr[sc_prot_gt1_bool_arr]
    aux_arr2 = np.ones_like(np.arange(0, len(sc_prot_gt1_list)))*(-1)

    
    # Replace protein by its number of peptides
    sc_prot_npep_gt1_list = [[protein2npep[j] for j in i] for i in sc_prot_gt1_list]

    # Get List of pairs (index of elem, position of maximum protein)
    sc_prot_npep_gt1_1max = [[n, np.argmax(i)] for n, i in enumerate(sc_prot_npep_gt1_list) if np.sum(np.max(i) == np.array(i)) == 1]
    sc_prot_npep_gt1_1max_index = [i for i,j in sc_prot_npep_gt1_1max]
    sc_prot_npep_gt1_1max_position = [j for i,j in sc_prot_npep_gt1_1max]
    
    # Add solved position to aux_arr2
    aux_arr2[sc_prot_npep_gt1_1max_index] = sc_prot_npep_gt1_1max_position

    # Resolve element with more than one maximum. We now use number of scans instead of number of peptides
    sc_prot_npep_gt1_gt1max_index = np.arange(0, len(sc_prot_gt1_list))
    sc_prot_npep_gt1_gt1max_index[sc_prot_npep_gt1_1max_index] = -1
    sc_prot_npep_gt1_gt1max_index = sc_prot_npep_gt1_gt1max_index[sc_prot_npep_gt1_gt1max_index != -1].tolist()

    sc_prot_nscan_gt1_list = \
        [[protein2nscan[j] for j in sc_prot_gt1_list[i]] for i in sc_prot_npep_gt1_gt1max_index]

    sc_prot_nscan_gt1_1max = [np.argmax(i) for i in sc_prot_nscan_gt1_list] #if np.sum(np.max(i) == np.array(i)) == 1]

    # Add to aux_arr2, the position of elements resolved using PSMs (or arbitrary position)
    aux_arr2[sc_prot_npep_gt1_gt1max_index] = sc_prot_nscan_gt1_1max

    # Incorporate to aux_arr its subset aux_arr2
    aux_arr[sc_prot_gt1_bool_arr] = aux_arr2

    # Incorporate to df_index_result_pos its subset aux_arr
    df_index_result_pos[non_decoy_bool] = aux_arr

    # Generate new columns
    mppSuffix = "_MPP"

    new_columns_list = [[j.split(" // ")[k] if k!=-1 else "" for j,k in zip(dffile[i].to_list(), df_index_result_pos)] \
        for i in semicolon_col_work_list]

    new_columns_df = pd.DataFrame({i+mppSuffix: j for i,j in zip(semicolon_col_work_list, new_columns_list)})

    # Remove > from the first character in protein column
    if semicolon_col_protein_list[1]+mppSuffix in new_columns_df.columns:
        new_columns_df[semicolon_col_protein_list[1]+mppSuffix] = new_columns_df[semicolon_col_protein_list[1]+mppSuffix].str.replace('^>', '', regex=True)
    
    
    # Generate final dataframe
    dffile_MPP = pd.concat([dffile.reset_index(drop=True), new_columns_df.reset_index(drop=True)], axis=1)
    
    return dffile_MPP

### END: RAFA MODIFICATION ###

    
def main(Fastapath,nm,decoy_tag,minpeptidelength,maxpeptidelength,aa_cut,not_cut,isoleu,enzymesuffix,Filepath,pepcolumn,collogic,FDO,FDS,PDO,PDS,PHR,num_threads):
        
    sh=expand(Fastareader(Fastapath,nm,decoy_tag,minpeptidelength,maxpeptidelength,aa_cut,not_cut,isoleu))
    sh.sort(key=lambda x: x[-2])
    sh=[ k+[" // ".join(i) for i in list(zip(*g))[:7]] for k, g in itertools.groupby(sh, lambda x: x[-2:])]
    [i.insert(2,len(i[4].split(" // "))) for i in sh]
    [i.insert(3,"Proteo") if i[2]==1 else i.insert(3,"NoProteo") for i in sh]
    col=["MissCleavage","ProteinsNumber","Proteotypic","UniProtIDs","Proteins","ProteinsLength","NpMiss","NpTotal","PeptidePos","PrevNextaa"]
    dfcol=[pepcolumn]+col
    dffasta=pd.DataFrame(sh,columns=dfcol)
    collogic=[pepcolumn]+[i for (i, v) in zip(col, collogic) if v=="1"]
    
    if  PDO=="1":

        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            
            executor.map(joiner,Filepath,repeat(dffasta[collogic]),repeat(PDO),repeat(PDS),repeat(PHR),repeat(pepcolumn)) 
            
        
    if  FDO=="1":
        
        
        dffasta.to_csv(os.path.splitext(Fastapath)[0]+enzymesuffix+FDS+".tsv",sep="\t",index=False)
    
    
   

    
if __name__ == '__main__':
    
    multiprocessing.freeze_support()  
    start = time.time()
    
    parser = argparse.ArgumentParser()    
    parser.add_argument("-p", "--params", help="params file",type=Path)
    args = parser.parse_args()
    paramsfile= str(args.params)
    
    #########get params###############
    
   
    
    with open(paramsfile ,"r") as f:   
        x = f.read().splitlines()
    
    defaults=[(" ",""),("\t",""),("add_column_",""),("decoy_tag=#","decoy_tag=NODECOY#"),("Isoleucine_Leucine_aa_Symbol=#","Isoleucine_Leucine_aa_Symbol=NOCHANGE#"),
       ("enzyme_name=#","enzyme_name=x#"),("Fasta_digested_suffix=#","Fasta_digested_suffix=digested#"),("PeptideFile_suffix=#","PeptideFile_suffix=over#")
       ,("Column_names_row=#","Column_names_row=1#")]
    for i in defaults:
        x=[a.replace(*i) for a in x]

    x=dict(list(filter(None,[list(filter(None,s.split("#")[0].split("="))) for s in x])))   
    num_threads=int(x["num_threads"])
    Fastapath=x["Fasta_path"]
    FDO=x["Fasta_digested_outputfile"]
    FDS=x["Fasta_digested_suffix"]
    FDS=("_"+FDS)
    Filepath=glob.glob(x["PeptideFile_path"])
    PDO=x["PeptideFile_outputfile"]
    PDS=x["PeptideFile_suffix"]
    PDS=("_"+PDS).replace("_over","")
    PHR=int(x["Column_names_row"])-1
    nm=int(x["allowed_missed_cleavage"])
    decoy_tag=x["decoy_tag"]
    minpeptidelength=int(x["digest_min_length"])
    maxpeptidelength=int(x["digest_max_length"])						 		
    enzymesuffix=x["enzyme_name"]
    enzymesuffix=("_"+enzymesuffix).replace("_x","")
    aa_cut=list(x["enzyme_cutafter"])
    not_cut=list(x["enzyme_butnotafter"])
    pepcolumn=x["Peptide_column_name"]
    isoleu=x["Isoleucine_Leucine_aa_Symbol"]
    
    collogic=[x["MissCleavage"],x["ProteinsNumber"],x["Proteotypic"],x["UniProtIDs"],x["Proteins"],x["ProteinsLength"],x["NpMiss"]
              ,x["NpTotal"],x["PeptidePos"],x["PrevNextaa"]]

    main(Fastapath,nm,decoy_tag,minpeptidelength,maxpeptidelength,aa_cut,not_cut,isoleu,enzymesuffix,Filepath,pepcolumn,collogic,FDO,FDS,PDO,PDS,PHR,num_threads)
    
    
    end = time.time()
    timer=divmod(end-start,60)  
    etxt1="FINISHED FASTA DIGESTION IN "+str(timer[0])+" MINUTES AND "+ str(round(timer[1],4))+" SECONDS"
    print(etxt1)
