import pymysql
import os
import sys
import re

def read_fna(file):
    seq = ""
    
    with open(file,"r") as f:
        for string in f:
            n = string.strip('\n')
            seq += n
    return seq       



dict_dna = {'A': 135, 'T' : 126, 'C' : 111, 'G' : 151}
dict_2chain = {'A': 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}
dict_rna = {'A': 135, 'U' : 112, 'C' : 111, 'G' : 151}
class origin:
    dicti = dict_dna
    def __init__(self, seq):
        self.chain = seq
        self.len = len(seq)
        
    def weight(self):
        weight = 0
        for residue in self.chain:
            weight += self.dicti[residue]
        return weight  
   
            
class Class_dna(origin):
    dicti = dict_dna

        
    def transcription(self):
        rna = ''
        for residue in self.chain:
            if residue == 'T':
                rna += 'U'
            else:
                rna += residue
        return Class_rna(rna)    
    def second_chain(self):
        chain2 = ''
        for residue in self.chain:
            chain2 += dict_2chain[residue]
        return Class_dna(chain2)
class Class_rna(origin):
    dicti = dict_rna
         
    def translation(self):
        translation = { 
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T', 
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',                  
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P', 
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R', 
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A', 
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G', 
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S', 
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L', 
        'UAC':'Y', 'UAU':'Y', 'UAA':'', 'UAG':'', 
        'UGC':'C', 'UGU':'C', 'UGA':'', 'UGG':'W', 
            } 
        protein = ''
        rna = self.chain
        i = 0
        for i in range(0,len(rna)-3,3):
            protein += translation[rna[i:i + 3]]
        return Class_protein(protein) 

class Class_protein(origin):
    dicti = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
           'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
           'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
           'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 , '_' : 0}
    
    
def find_orf(seq):
    end = ['TAG', 'TAA', 'TGA']
    i = 0
    dna = ''
    orf = ''
    j = 0
    k = 0
    result = []
    for k in range(3):
        orf = ''
        for i in range(k, len(seq)-3):
            if seq[i:i + 3] == 'ATG':
                j = i + 3
                l = 0
                while j < len(seq) - 3:
                    orf += seq[j:j+3]
                    if seq[j:j + 3]  in end:
                        l = 1
                        break
                    j += 3    
                    
                if len(orf) > 100 and l ==1:
                    result.append(orf)
                
    
    return result    
    
#file = sys.argv[1]
#seq = read_fna(file)
seq = read_fna("sample.fasta")
chromas = re.compile('>chr[^A-Z]{1,3}')
chromasomes = chromas.findall(seq)
conn = pymysql.connect('localhost', 'root','0000', 'mariadb' )
while True:
    print("Select a mode: 1 - create orf base and add orf from file to a base, 2 - extract part of dna or protein with weight")
    mode =  int(input())
    if mode == 1:
        break
    else:
        print("Which chromasome? ", chromasomes)
        chr_name = input()
        print("Enter start")
        start = int(input())
        print("Enter end")
        end = int(input())
        break

if mode == 1:
    with conn.cursor() as cur:
        chr_len = 0
        chr_name_coord = []
        cur.execute('''CREATE TABLE IF NOT EXISTS ORF_base (id int NOT NULL AUTO_INCREMENT, chromasome VARCHAR(200), DNA_orf LONGTEXT, DNA_weight FLOAT,RNA LONGTEXT, RNA_weight FLOAT, Protein LONGTEXT,Protein_weight FLOAT, PRIMARY KEY (id))''')
        for chrom_name in chromasomes:
            inf = []
            chromasome = re.compile(chrom_name)
            chr_inf = chromasome.search(seq)
            inf.append(chrom_name)
            inf.append(chr_inf.start())
            inf.append(chr_inf.end())
            chr_name_coord.append(inf)
        #chr_name_coord.append(['ewf', 12, 15])  
        #chr_name_coord.append(['ewfd', 111, 115])
        chr_coord = []
        for i in range(len(chr_name_coord) - 1):
            inf = []
            inf.append(chr_name_coord[i][0])
            inf.append(chr_name_coord[i][2] + 1)
            inf.append(chr_name_coord[i+1][1])
            chr_coord.append(inf)
        inf_end = []    
        inf_end.append(chr_name_coord[len(chr_name_coord) - 1][0])
        inf_end.append(chr_name_coord[len(chr_name_coord) - 1][2])
        inf_end.append(len(seq))
        chr_coord.append(inf_end)    
        
        for chrom in chr_coord:
            chr_seq = seq[chrom[1]:chrom[2]]
            chain1 = Class_dna(chr_seq)
            chain2 = chain1.second_chain()
            orf1 = find_orf(chain1.chain)
            orf2 = find_orf(chain2.chain)
            for el in orf1:
                dna = Class_dna(el)
                rna = dna.transcription()
                protein = rna.translation()
                dna_chain = '{}'.format(dna.chain)
                rna_chain = '{}'.format(rna.chain)
                protein_chain = '{}'.format(protein.chain)
                
                insert = 'INSERT INTO ORF_base (chromasome, DNA_orf, DNA_weight, RNA, RNA_weight, Protein, Protein_weight) VALUES (%s, %s, %s, %s,  %s, %s, %s)'
                
                lisst = [chrom[0], dna_chain, dna.weight(), rna_chain, rna.weight(), protein_chain, protein.weight()]
                cur.execute(insert, lisst) #, chrom[0], dna_chain, dna.weight(), rna_chain, rna.weight(), protein_chain, protein.weight())
                
          
                conn.commit()
                
            for el in orf2: 
                dna = Class_dna(el)
                rna = dna.transcription()
                protein = rna.translation()
                dna_chain = '{}'.format(dna.chain)
                rna_chain = '{}'.format(rna.chain)
                protein_chain = '{}'.format(protein.chain)
                
                insert = 'INSERT INTO ORF_base (chromasome, DNA_orf, DNA_weight, RNA, RNA_weight, Protein, Protein_weight) VALUES (%s, %s, %s, %s,  %s, %s, %s)'
                
                lisst = [chrom[0], dna_chain, dna.weight(), rna_chain, rna.weight(), protein_chain, protein.weight()]
                cur.execute(insert, lisst)
                dna = Class_dna(el)
                rna = dna.transcription()
                protein = rna.translation()
                dna_chain = '{}'.format(dna.chain)
                rna_chain = '{}'.format(rna.chain)
                protein_chain = '{}'.format(protein.chain)
                
                conn.commit()
else:
    with conn.cursor() as cur:
        find = 'SELECT DNA from DNA_base WHERE chromasome=%s'
        
        dna_find = cur.execute(find, chr_name)
        f = cur.fetchall()
        for el in f:
            seq = Class_dna(el[0][start:end])
            rna = seq.transcription()
            protein = rna.translation()
            print("DNA {0}, weight: {1}".format(seq.chain, seq.weight()))
            print("RNA {0}, weight: {1}".format(rna.chain, rna.weight()))
            print("DNA {0}, weight: {1}".format(protein.chain, protein.weight()))
            
