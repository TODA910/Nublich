from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 
from modules.extraction.interpreter import get_file_format
from modules.serializer.storer import phred_read_store, length_read_store, barcode_read_store, info_read_store
import os
import numpy as np
import gzip
import re 
import time

def seq_filter(filename, seq, max = 1000000):
    start_time = time.time()
    count = 0
    stem = os.path.splitext(filename)[0]
    file_input = filename
    format = get_file_format(filename)
    if format == "gz":
        file_input = gzip.open(filename,"rt")
        ext = ""
    else:
        ext = ".fastq"
    records = parse(file_input, "fastq") 
    if len(seq) >6:
        Name = f"{seq[:3]}...{seq[-3:]}"
    else:
        Name = f"{seq}"
    with gzip.open(f"{Name}_{stem}{ext}.gz", "wt") as file_output:
        for idx, record in enumerate(records,1):
            if bool(re.search(seq, str(record.seq))):
                count +=1
                file_output.write(record.format("fastq"))
            if idx == max:
                break
            if idx % 10000 == 0:
                print(f"{idx} reads have been scanned successfully ", "--- %s seconds ---" % (time.time() - start_time))
                start_time = time.time()
        print(f"{count} sequences matched the given sequence")
    file_output.close()
    print("File created successfully")
    

def barcode_filter(filename, max = 1000000):
    start_time = time.time()
    stem = os.path.splitext(filename)[0]
    file_input = filename
    format = get_file_format(filename)
    if format == "gz":
        file_input = gzip.open(filename,"rt")
        ext = ""
    else:
        ext = ".fastq"
    Barcode_list = info_read_store(filename, max)[3]
    if len(Barcode_list) == 0 or len(Barcode_list) == 1:
        print(f"There is already a single barcode in this file")
        if format == "gz":
            file_input.close()
    else: 
        records = parse(file_input, "fastq") 
        print("Creating filtered files...")
        for idx, record in enumerate(records,1):
            m = re.search(f"barcode=(\w+)", record.description)
            with gzip.open(f"{m.group(1)}_{stem}{ext}.gz", "at") as file_output:
                file_output.write(record.format("fastq"))
                file_output.close() 
            if idx == max:
                break
            if (idx+1) % 10000 == 0:
                print(f"{idx+1} reads have been scanned successfully ", "--- %s seconds ---" % (time.time() - start_time)) 
                start_time = time.time()
        if format == "gz":
            file_input.close() 
        
# def barcode_filter(filename, max = 1000000):
#     stem = os.path.splitext(filename)[0]
#     file_input = filename
#     format = get_file_format(filename)
#     if format == "gz":
#         file_input = gzip.open(filename,"rt")
#         ext = ""
#     else:
#         ext = ".fastq"
#     Barcode_list = phred_read_store(filename, 1000000)
#     if len(Barcode_list) == 0 or 1:
#         print(f"There is already a single barcode in this file")
#     else: 
        
#         for idx, record in enumerate(records,1):
#              m = re.finditer(f"barcode=(\w+)", record.description)

#         for i in Barcode_list:
#             Name = i
#             with open(f"{filename}", 'rt') as file_input:
#                 records = parse(file_input, "fastq") 
#                 with gzip.open(f"{Name}_{stem}{ext}", "wt") as file_output:
#                     for idx, record in enumerate(records,1):
#                         m = re.finditer(f"barcode=(\w+)", record.description)
#                         if m.group(1) == i:
#                             file_output.write(record)
#                         elif:
#                             break
#                         if idx == max:
#                             break
#                         if idx % 10000 == 0:
#                                 print(f"{idx} reads have been scanned successfully")
#             with open(f"{Name}_{stem}{ext}", 'rb') as f_in:
#                 with gzip.open(f"{Name}_{stem}{ext}.gz", 'wb') as f_out:
#                     shutil.copyfileobj(f_in, f_out)


def phred_filter(filename, max = 1000000, phred_interval = [10, 90]):
    ##Limit interval
    start_time = time.time()
    stem = os.path.splitext(filename)[0]
    file_input = filename
    format = get_file_format(filename)
    if format == "gz":
        file_input = gzip.open(filename,"rt")
        ext = ""
    else:
        ext = ".fastq"
    print(phred_interval)
    Name = f"phred_[{str(phred_interval[0])}-{str(phred_interval[1])}]"
    phred_list = info_read_store(filename, max)[1]
    permin = np.percentile(phred_list, phred_interval[0])
    permax = np.percentile(phred_list, phred_interval[1])
    records = parse(file_input, "fastq") 
    with gzip.open(f"{Name}_{stem}{ext}.gz", "wt") as file_output:
        print("Creating filtered files...")
        for idx, record in enumerate(records,1):
            phred_score_list = record.letter_annotations["phred_quality"]
            phred_avg = sum(phred_score_list)/len(phred_score_list)
            if phred_avg >= permin and phred_avg <= permax:
                file_output.write(record.format("fastq"))
            if idx == max:
                break
            if idx % 10000 == 0:
                print(f"{idx} reads have been scanned successfully ", "--- %s seconds ---" % (time.time() - start_time))
                start_time = time.time()
    file_output.close()
    print("File created successfully")

def phred_filter_val(filename, max = 1000000, phred_interval = [0, 40]):
    ##Limit interval
    start_time = time.time()
    stem = os.path.splitext(filename)[0]
    file_input = filename
    format = get_file_format(filename)
    if format == "gz":
        file_input = gzip.open(filename,"rt")
        ext = ""
    else:
        ext = ".fastq"
    records = parse(file_input, "fastq") 
    Name = f"phred_[{str(phred_interval[0])}-{str(phred_interval[1])}]_val"
    permin = phred_interval[0]
    permax = phred_interval[1] 
    with gzip.open(f"{Name}_{stem}{ext}.gz", "wt") as file_output:
        print("Creating filtered files...")
        for idx, record in enumerate(records,1):
            phred_score_list = record.letter_annotations["phred_quality"]
            phred_avg = sum(phred_score_list)/len(phred_score_list)
            if phred_avg >= permin and phred_avg <= permax:
                file_output.write(record.format("fastq"))
            if idx == max:
                break
            if idx % 10000 == 0:
                print(f"{idx} reads have been scanned successfully ", "--- %s seconds ---" % (time.time() - start_time))
                start_time = time.time()
    file_output.close()
            
def length_filter(filename, max = 1000000, Length_interval = [10, 90]):
    ##Limit interval
    start_time = time.time()
    stem = os.path.splitext(filename)[0]
    file_input = filename
    format = get_file_format(filename)
    if format == "gz":
        file_input = gzip.open(filename,"rt")
        ext = ""
    else:
        ext = ".fastq"
    records = parse(file_input, "fastq") 
    Name = f"Length_[{str(Length_interval[0])}-{str(Length_interval[1])}]"
    length_list = info_read_store(filename, max)[2]
    permin = np.percentile(length_list, Length_interval[0])
    permax = np.percentile(length_list, Length_interval[1])
    with gzip.open(f"{Name}_{stem}{ext}.gz", "wt") as file_output:
        print("Creating filtered files...")
        for idx, record in enumerate(records,1):
            if len(record) >= permin and len(record) <= permax:
                file_output.write(record.format("fastq"))
            if idx == max:
                break
            if idx % 10000 == 0:
                print(f"{idx} reads have been scanned successfully ", "--- %s seconds ---" % (time.time() - start_time))
                start_time = time.time()
    file_output.close()
    print("File created successfully")

def length_filter_val(filename, max = 1000000, phred_interval = [0, 100000]):
    ##Limit interval
    start_time = time.time()
    stem = os.path.splitext(filename)[0]
    file_input = filename
    format = get_file_format(filename)
    if format == "gz":
        file_input = gzip.open(filename,"rt")
        ext = ""
    else:
        ext = ".fastq"
    records = parse(file_input, "fastq") 
    Name = f"Length_[{str(phred_interval[0])}-{str(phred_interval[1])}]_val"
    permin = phred_interval[0]
    permax = phred_interval[1]
    with gzip.open(f"{Name}_{stem}{ext}.gz", "wt") as file_output:
        print("Creating filtered files...")
        for idx, record in enumerate(records,1):
            if len(record) >= permin and len(record) <= permax:
                file_output.write(record.format("fastq"))
            if idx == max:
                break
            if idx % 10000 == 0:
                print(f"{idx} reads have been scanned successfully ", "--- %s seconds ---" % (time.time() - start_time))
                start_time = time.time()
    file_output.close()
    print("File created successfully")

def phred_info(filename, max = 1000000):
    # print(info_read_store(filename,max))
    phred_list = info_read_store(filename, max)[5]
    # print(phred_list)
    amin = phred_list[0]
    per10 = phred_list[1]
    per25 = phred_list[2]
    per50 = phred_list[3]
    per75 = phred_list[4]
    per90 = phred_list[5]
    amax = phred_list[6]
    avg = phred_list[7]

    print("\n---- Basic Phred Information ----")
    print (f"Lowest value: {amin:.03f}")
    print (f"10th percentile: {per10:.03f}")
    print (f"25th percentile: {per25:.03f}")
    print (f"50th percentile: {per50:.03f}")
    print (f"75th percentile: {per75:.03f}")
    print (f"90th percentile: {per90:.03f}")
    print (f"Highest value: {amax:.03f}")
    print (f"Average value: {avg:.03f}")
    # return [per10,per25,per50,per75,per90]

def length_info(filename, max = 1000000):
    length_list = info_read_store(filename, max)[6]
    # print(length_list)
    amin = length_list[0]
    per10 = length_list[1]
    per25 = length_list[2]
    per50 = length_list[3]
    per75 = length_list[4]
    per90 = length_list[5]
    amax = length_list[6]
    avg = length_list[7]
    print (f"\n---- Basic Length Information ----")
    print (f"Lowest value: {amin:.01f}")
    print (f"10th percentile: {per10:.01f}")
    print (f"25th percentile: {per25:.01f}")
    print (f"50th percentile: {per50:.01f}")
    print (f"75th percentile: {per75:.01f}")
    print (f"90th percentile: {per90:.01f}")
    print (f"Highest value: {amax:.01f}")
    print (f"Average value: {avg:.01f}")
    # return [per10,per25,per50,per75,per90]



    
