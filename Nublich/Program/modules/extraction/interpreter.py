from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 
import gzip
import re 
import time
import numpy as np

def get_file_format(filename):
    format_search = re.search(r'\.(\w+)$',filename)
    file_format = format_search.group(1)
    return (file_format)
    
def report():
    print("ok") 

def get_num_sequence(filename, max=10):
    start_time = time.time()
    file_input = filename
    format = get_file_format(filename)
    reads = 0
    bases = 0
    if format == "gz":
        file_input = gzip.open(filename,"rt")
    records = parse(file_input, "fastq") 
    for idx, record in enumerate(records,1):
        reads +=1
        bases += len(record)
        if idx == max:
            break
        if (idx+1) % 10000 == 0:
            print(f"{idx+1} reads have been scanned successfully ", "Total runtime : {time_used:.4f} seconds".format(time_used = (time.time() - start_time)))
            start_time = time.time()
    return [reads,bases]

def get_file_info(filename,max=10,first=True):
    barcode_set = set()
    length = []
    phred_score = []
    start_time = time.time()
    file_input = filename
    format = get_file_format(filename)
    column_set = set()
    reads = 0
    bases = 0
    if format == "gz":
        file_input = gzip.open(filename,"rt")
    records = parse(file_input, "fastq") 
    print("Scanning through files...")
    for idx, record in enumerate(records,1):
        for m in re.finditer(r" (\w+)=(\w+)", record.description):
            column_set.add(m.group(1))
            if m.group(1) == "barcode":
                barcode_set.add(m.group(2))
        length.append(len(record))
        phred_score_list = record.letter_annotations["phred_quality"]
        phred_avg = round(sum(phred_score_list)/len(phred_score_list),3)
        phred_score.append(phred_avg)
        reads+=1
        bases+=len(record)
        if idx == max:
            break
        if (idx+1) % 10000 == 0 and first:
            print(f"{idx+1} reads have been scanned successfully ", "Total runtime : {time_used:.4f} seconds".format(time_used = (time.time() - start_time)))
            start_time = time.time()
    column_set.update(["mean_phred_score","no.","read_length"])
    phred_stat = [np.amin(phred_score),np.percentile(phred_score,10),np.percentile(phred_score,25),np.percentile(phred_score,50),np.percentile(phred_score,75), np.percentile(phred_score,90), np.max(phred_score), np.average(phred_score)]
    length_stat = [np.amin(length),np.percentile(length,10),np.percentile(length,25),np.percentile(length,50),np.percentile(length,75), np.percentile(length,90), np.max(length), np.average(length)]    
    info = [sorted(column_set),phred_score,length,barcode_set,[reads,bases],phred_stat,length_stat]
    # for i, a in enumerate(info,0):
    #     print (i,a)
    return info

# def add_column(file_name):
def get_phred_score(filename,max=10):
    start_time = time.time()
    phred_score = []
    file_input = filename
    format = get_file_format(filename)
    if format == "gz":
        file_input = gzip.open(filename,"rt")
    records = parse(file_input, "fastq") 
    print("Scanning through files..")
    for idx, record in enumerate(records,1):
        phred_score_list = record.letter_annotations["phred_quality"]
        phred_avg = float(f"{sum(phred_score_list)/len(phred_score_list):.03f}")
        phred_score.append(phred_avg)
        if idx == max:
            break
        if idx % 10000 == 0:
            print(f"{idx} reads have been scanned successfully ", "Total runtime : {time_used:.4f} seconds".format(time_used = (time.time() - start_time)))
            start_time = time.time()
    return phred_score

def get_length(filename,max=10):
    start_time = time.time()
    length = []
    file_input = filename
    format = get_file_format(filename)
    if format == "gz":
        file_input = gzip.open(filename,"rt")
    records = parse(file_input, "fastq") 
    print("Scanning through files..")
    for idx, record in enumerate(records,1):
        length.append(len(record))
        if idx == max:
            break
        if idx % 10000 == 0:
            print(f"{idx} reads have been scanned successfully ", "Total runtime : {time_used:.4f} seconds".format(time_used = (time.time() - start_time)))
            start_time = time.time()
    return length

def get_barcode(filename,max=10):
    start_time = time.time()
    barcode_set = set()
    file_input = filename
    format = get_file_format(filename)
    if format == "gz":
        file_input = gzip.open(filename,"rt")
    records = parse(file_input, "fastq") 
    print("Scanning through files..")
    for idx, record in enumerate(records,1):
        m = re.search(f"barcode=(\w+)", record.description)
        barcode_set.add(m.group(1))
        if idx == max:
            break
        if idx % 10000 == 0:
            print(f"{idx} reads have been scanned successfully ", "Total runtime : {time_used:.4f} seconds".format(time_used = (time.time() - start_time)))
            start_time = time.time()
    return list(barcode_set)

def extract_info(filename, chocol, max = 10):
    start_time = time.time()
    file_input = filename
    format = get_file_format(filename)
    column_dict = {}
    for i in chocol:
        column_dict[i] = []
    if format == "gz":
        file_input = gzip.open(filename,"rt")
    records = parse(file_input, "fastq") 
    for idx, record in enumerate(records, 1): 
        for k in set(chocol).difference(set(["mean_phred_score","no.","read_length"])):
            s = re.search(f"{k}=(\w+)", record.description)
            if bool(s):
                column_dict[k].append(s.group(1))
            else:
                column_dict[k].append(" ")
        if "no." in chocol:
            column_dict["no."].append(idx)
        if "mean_phred_score" in chocol:
            phred_score_list = record.letter_annotations["phred_quality"]
            phred_avg = sum(phred_score_list)/len(phred_score_list)
            column_dict["mean_phred_score"].append(f"{phred_avg:.03f}")
        if "read_length" in chocol:
            column_dict["read_length"].append(len(record))
        if idx == max:
            break
        if idx % 10000 == 0:
            print(f"{idx} reads have been scanned successfully ", "Total runtime : {time_used:.4f} seconds".format(time_used = (time.time() - start_time)))
            start_time = time.time()
    print ("Extraction Complete")
    return column_dict

if __name__ == "__main__":
# def report(dict_to_print,column_name):
#     tem = ""
#     for i in dict_to_print.keys():
#         for j in column_name:
#         dict_to_print[i]
#         tem = ""
#         print(tem)
    
    # print(f"{Name:}")
    # for i in dict_to_print.keys():
# print(get_file_format("ont.exp2.fastq.gz"))
# print(extract_info("ont.exp2.fastq.gz",100))
# extract_info("ont.exp2.fastq.gz",30000)
# print(extract_info("example.fastq",100))
    print(get_file_info("ont.exp2.fastq.gz",100000) )
# print(get_file_col_des("example.fastq", 5))
# print(get_num_sequence("example.fastq"))
# print(get_num_sequence("ont.exp2.fastq.gz"))
# print("ok")





