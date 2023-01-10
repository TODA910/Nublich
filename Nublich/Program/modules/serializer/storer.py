import os
import pickle
from modules.extraction.interpreter import get_file_info, get_phred_score, get_barcode, get_length, get_file_info

def info_read_store(filename, reads):    
    subject = "info"
    ##Load last file
    if os.path.isfile(f"{filename}_{subject}.pickle"):
        # print(f"Reading {filename}_{subject}.pickle")
        with open(f"{filename}_{subject}.pickle", 'rb') as input:
            previous_info = pickle.load(input)
            num_reads = previous_info[0]
            info = previous_info[1]
            # print(previous_info)
        input.close()
    else:
        # print("ERROR")
        with open(f'{filename}_{subject}.pickle', 'wb') as outputfile:
            print(f"Creating {filename}_{subject}.pickle")
            var = list(get_file_info(filename, reads))
            pickle.dump([reads, var], outputfile)
        outputfile.close()
        with open(f"{filename}_{subject}.pickle", 'rb') as input:
            previous_info = pickle.load(input)
            num_reads = previous_info[0]
            info = previous_info[1]
        input.close()
    if reads != num_reads:
        # print("ERROR")
        with open(f'{filename}_{subject}.pickle', 'wb') as outputfile:
            var = list(get_file_info(filename, reads))
            pickle.dump([reads, var], outputfile)
        outputfile.close()
        with open(f"{filename}_{subject}.pickle", 'rb') as input:
            previous_info = pickle.load(input)
            num_reads = previous_info[0]
            info = previous_info[1]
        input.close()
    return(info)

# def col_read_store(filename, reads):    
#     subject = "column"
#     ##Load last file
#     if os.path.isfile(f"{filename}_{subject}.pickle"):
#         with open(f"{filename}_{subject}.pickle", 'rb') as input:
#             previous_info = pickle.load(input)
#             num_reads = previous_info[0]
#             info = previous_info[1]
#             # print(previous_info)
#         input.close()
#     else:
#         with open(f'{filename}_{subject}.pickle', 'wb') as outputfile:
#             var = list(get_file_col_des(filename, reads))
#             pickle.dump([reads, var], outputfile)
#         outputfile.close()
#         with open(f"{filename}_{subject}.pickle", 'rb') as input:
#             previous_info = pickle.load(input)
#             num_reads = previous_info[0]
#             info = previous_info[1]
#         input.close()
#     if reads > num_reads:
#         with open(f'{filename}_{subject}.pickle', 'wb') as outputfile:
#             var = list(get_file_col_des(filename, reads))
#             pickle.dump([reads, var], outputfile)
#         outputfile.close()
#     return(info)

def phred_read_store(filename, reads):    
    subject = "phred"
    ##Load last file
    if os.path.isfile(f"{filename}_{subject}.pickle"):
        with open(f"{filename}_{subject}.pickle", 'rb') as input:
            previous_info = pickle.load(input)
            num_reads = previous_info[0]
            info = previous_info[1]
            # print(previous_info)
        input.close()
    else:
        with open(f'{filename}_{subject}.pickle', 'wb') as outputfile:
            var = list(get_phred_score(filename, reads))
            pickle.dump([reads, var], outputfile)
        outputfile.close()
        with open(f"{filename}_{subject}.pickle", 'rb') as input:
            previous_info = pickle.load(input)
            num_reads = previous_info[0]
            info = previous_info[1]
        input.close()
    if reads != num_reads:
        with open(f'{filename}_{subject}.pickle', 'wb') as outputfile:
            var = list(get_phred_score(filename, reads))
            pickle.dump([reads, var], outputfile)
        outputfile.close()
    return(info)

def length_read_store(filename, reads):    
    subject = "length"
    ##Load last file
    if os.path.isfile(f"{filename}_{subject}.pickle"):
        with open(f"{filename}_{subject}.pickle", 'rb') as input:
            previous_info = pickle.load(input)
            num_reads = previous_info[0]
            info = previous_info[1]
            # print(previous_info)
        input.close()
    else:
        with open(f'{filename}_{subject}.pickle', 'wb') as outputfile:
            var = list(get_length(filename, reads))
            pickle.dump([reads, var], outputfile)
        outputfile.close()
        with open(f"{filename}_{subject}.pickle", 'rb') as input:
            previous_info = pickle.load(input)
            num_reads = previous_info[0]
            info = previous_info[1]
        input.close()
    if reads != num_reads:
        with open(f'{filename}_{subject}.pickle', 'wb') as outputfile:
            var = list(get_length(filename, reads))
            pickle.dump([reads, var], outputfile)
        outputfile.close()
    return(info)

def barcode_read_store(filename, reads):
    subject = "barcode"
    ##Load last file
    if os.path.isfile(f"{filename}_{subject}.pickle"):
        with open(f"{filename}_{subject}.pickle", 'rb') as input:
            previous_info = pickle.load(input)
            num_reads = previous_info[0]
            info = previous_info[1]
            print(previous_info)
        input.close()
    else:
        with open(f'{filename}_{subject}.pickle', 'wb') as outputfile:
            var = list(get_barcode(filename, reads))
            pickle.dump([reads, var], outputfile)
        outputfile.close()
        with open(f"{filename}_{subject}.pickle", 'rb') as input:
            previous_info = pickle.load(input)
            num_reads = previous_info[0]
            info = previous_info[1]
        input.close()
    if reads > num_reads:
        with open(f'{filename}_{subject}.pickle', 'wb') as outputfile:
            var = list(get_barcode(filename, reads))
            pickle.dump([reads, var], outputfile)
        outputfile.close()
    return(info)


# def store_dict(filename, data):    
#     ##Load last file
#     if os.path.isfile("datadict.pickle"):
#         with open("datadict.pickle", 'rb') as input:
#             previous_info = pickle.load(input)
#             last_file = previous_info[0]
#             data_dict = previous_info[1]
#             # print(previous_info)
#         input.close()
#         if filename == last_file:
#             data = data_dict

#     else:
#         with open('datadict.pickle', 'wb') as outputfile:
#             pickle.dump([filename, data], outputfile)
#         outputfile.close()
#         with open("datadict.pickle", 'rb') as input:
#             previous_info = pickle.load(input)
#             last_file = previous_info[0]
#             data_dict = previous_info[1]
#         input.close()

#     if filename == last_file:
#         data = data_dict

#     else:
#         with open('datadict.pickle', 'wb') as outputfile:
#             col = list(get_file_col_des(filename, reads))
#             pickle.dump([filename, col], outputfile)
#         outputfile.close()
#     return(data)
# # read_store(filename, reads)  