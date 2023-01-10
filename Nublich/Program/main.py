from modules.extraction.interpreter import get_num_sequence, extract_info
from modules.serializer.storer import *
from modules.export.report_gen import *
from modules.extraction.filter import *
from modules.export.visualizer import *
import os
import pickle
import pandas as pd
import time
import shutil
import subprocess
import platform
import csv

def argparserLocal():
    from argparse import ArgumentParser
    parser = ArgumentParser(prog = "nublich", description= "Program to analyze your fastq/gz files")

    subparsers = parser.add_subparsers(
        title = "Available commands", dest = "command"
    )
    subparsers.reqiured = True

    crb_command = subparsers.add_parser("count", help = "Count total reads and bases")
    crb_command.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    crb_command.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")

    get_phred = subparsers.add_parser("getphredinfo", help = "Get phred score basic statistics")
    get_phred.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    get_phred.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")

    get_length = subparsers.add_parser("getlengthinfo", help = "Get base length basic statistics")
    get_length.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    get_length.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")

    get_command = subparsers.add_parser("getallinfo", help = "Get all basic information")
    get_command.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    get_command.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")

    seq_command = subparsers.add_parser("filter_seq", help = "Filter by base sequences")
    seq_command.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    seq_command.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")
    seq_command.add_argument("-s", "--seq", type=str, default=None, dest="sequence", help = "Insert specific base sequence to scan through")

    barcode_command = subparsers.add_parser("filter_barcode", help = "Filter by barcode")
    barcode_command.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    barcode_command.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")

    score_command = subparsers.add_parser("filter_score", help = "Filter specified range of phred score")
    score_command.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    score_command.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")
    score_command.add_argument("-v", "--var", action = "store_true", dest="value", help = "Percentile to absolute value")
    score_command.add_argument("--min", type=int, default=0, dest= "interval_min", help = "Insert minimum")
    score_command.add_argument("--max", type=int, default=None, dest= "interval_max", help = "Insert maximum")
    
    length_command = subparsers.add_parser("filter_length", help = "Filter specified range of base length")
    length_command.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    length_command.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")
    length_command.add_argument("-v", "--var", action = "store_true", dest="value", help = "Percentile to absolute value")
    length_command.add_argument("--min", type=int, default=0, dest= "interval_min", help = "Insert minimum")
    length_command.add_argument("--max", type=int, default=None, dest= "interval_max", help = "Insert maximum")

    report_command = subparsers.add_parser("overview", help = "See an overview table of a file")
    report_command.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    report_command.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")
    report_command.add_argument("-t", "--title", type=str, nargs = "*", default=None, dest="columns", help = "Insert column titles")

    get_report_command = subparsers.add_parser("get_report", help = "Generate a HTML report")
    get_report_command.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    get_report_command.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")
   
    convert_command = subparsers.add_parser("zip", help = "Zip fastq files into gz files")
    convert_command.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")

    csv_command = subparsers.add_parser("csv", help= "Convert a fastq/gz files into a csv files")
    csv_command.add_argument("-f", "--file", type=str, default=None, dest="filename", help = "Insert file")
    csv_command.add_argument("-r", "--reads", type=int, default=1000000, dest="max_reads", help = "Insert number of reads to scan")
    csv_command.add_argument("-t", "--title", type=str, nargs = "*", default=None, dest="columns", help = "Insert column titles")

    # subparsers.add_parser("get_report", help = "Get a HTML report of the last scanned file")
    check_command = subparsers.add_parser("check")
    check_command.add_argument("-f", "--file", dest="filename")

    subparsers.add_parser("clear", help ="Clear pickle files")

# report()
    return parser
def test():
    parser = argparserLocal()
    args = parser.parse_args(["report","-f", "example.fastq"])
    print(get_num_sequence(args.filename))

def main():
    start_time = time.time()
    parser = argparserLocal()
    args = parser.parse_args()
    # print(args)
    print(f"Running...")

    if args.command == "check":
        with open(args.filename, 'rb') as file_input:
            previous_info = pickle.load(file_input)
            last_file = previous_info[0]
            column = previous_info[1]
        file_input.close()
        print(f"Reads: {last_file} \nAvailable info: {column}")

    if args.command == "clear":
        my_dir = "."
        deleted_files, all_files = 0, 0
        for i in os.listdir(my_dir):
            if i.endswith("pickle"):
                print(i)
                all_files +=1
                a = True
                while a:
                    x = input("Do you want to delete this file? [y/n]: ")
                    if x == "y":
                        os.remove(os.path.join(my_dir, i))
                        deleted_files += 1
                        a = False
                    elif x == "n":
                        a = False
                        break
        if all_files == 0:
            print(f"There are currently no pickle files.")
        else: 
            print(f"{deleted_files} out of {all_files} pickle files have been deleted.")

    if args.command == "zip":
        if args.filename == None:
            exit(parser.parse_args(['count','-h']))
        stem = os.path.splitext(args.filename)[0]
        file_input = args.filename
        format = get_file_format(args.filename)
        if format == "gz":
            exit(print("The input file is already a gz file!"))
        else:
            ext = ".fastq"
            with open(f"{stem}{ext}", 'rb') as f_in:
                with gzip.open(f"{stem}{ext}.gz", 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
    
    if args.command == "csv":
        if args.filename == None:
            exit(parser.parse_args(['csv','-h']))
        col = info_read_store(args.filename,args.max_reads)[0]
        if args.columns == None:
            # exit(print(f"Please enter titles from the following: {col}"))
            chosen_column = col
        else:
            chosen_column = args.columns
            for i in chosen_column:
                if i not in col:
                    exit(print(f"{i} cannot be a title \nPlease choose from the following: {col}"))
        extracted_info = extract_info(args.filename, chosen_column, args.max_reads)
        num = len(list(extracted_info.values())[0])
        csv_file = open("data.csv", "w")
        writer = csv.writer(csv_file)
        writer.writerow([i for i in extracted_info.keys()])
        for k in range(0,num):
            value_list = []
            for i in extracted_info.values():
                value_list.append(i[k])
            writer.writerow(value_list)
        

    if args.command == "count":
        if args.filename == None:
            exit(parser.parse_args(['count','-h']))
        reads, bases = info_read_store(args.filename, args.max_reads)[4]
        print(f"Number of reads: {reads}, Number of bases: {bases}")

    if args.command == "get_report":
        if args.filename == None:
            exit(parser.parse_args(['get_report','-h']))
        col = info_read_store(args.filename,args.max_reads)[0]
        extracted_info = extract_info(args.filename, col, args.max_reads)
        df = pd.DataFrame(extracted_info)
        sum_df = get_sum_df(df)
        barcodesum_df = get_barcodesum_df(df)
        con1 = getTable(sum_df, 'General Summary')
        con2 = getTable(barcodesum_df, 'Barcode Summary')
        con3 = graph_pie(df, 'Barcode Percentage')
        con4 = graph_length_density(df, 'Length Distribution')
        con5 = graph_phred_density(df, 'Mean PHRED Distribution', info_read_store(args.filename,args.max_reads)[5])
        con6 = graph_length_phred(df, 'Length vs PHRED quality')
        createReport(con1,con2,con3,con4,con5,con6,filename = args.filename)
        path = f"{args.filename}_report.html"
        current_platform = platform.system()
        if current_platform == "Linux":
            subprocess.call(["xdg-open", path])
        elif current_platform == "Windows":
            os.system("start "+path)
        elif current_platform == "Darwin":
            subprocess.call(["open", path])

    if args.command == "getallinfo":
        if args.filename == None:
            exit(parser.parse_args(['getallinfo','-h']))
        info = info_read_store(args.filename,args.max_reads)
        col = info[0]
        if len(info[3]) < 2:
            barcode = 1
        elif len(info[3]) >= 2:
            barcode = len(info[3])
        else:
            barcode = "ERROR"
        reads, bases = info[4]
        print(f"\n---- Available information ----")
        print(f"Number of reads: {reads}, Number of bases: {bases}")
        print(f"Available data description: {col}")
        print(f"Number of barcodes: {barcode}")
        phred_info(args.filename,args.max_reads)
        length_info(args.filename,args.max_reads)

    if args.command == "getphredinfo":
        if args.filename == None:
            exit(parser.parse_args(['getphredinfo','-h']))
        phred_info(args.filename,args.max_reads)

    if args.command == "getlengthinfo":
        if args.filename == None:
            exit(parser.parse_args(['getlengthinfo','-h']))
        length_info(args.filename,args.max_reads)

    if args.command == "overview":
        if args.filename == None:
            exit(parser.parse_args(['overview','-h']))
        col = info_read_store(args.filename,args.max_reads)[0]
        if args.columns == None:
            # exit(print(f"Please enter titles from the following: {col}"))
            chosen_column = col
        else:
            chosen_column = args.columns
            for i in chosen_column:
                if i not in col:
                    exit(print(f"{i} cannot be a title \nPlease choose from the following: {col}"))
        extracted_info = extract_info(args.filename, chosen_column, args.max_reads)
        # print(extracted_info)
        df = pd.DataFrame(extracted_info)
        print(df)

    if args.command == "filter_seq":
        if args.filename == None:
            exit(parser.parse_args(["filter_seq",'-h']))
        seq_filter(args.filename, args.sequence, args.max_reads)
    
    if args.command == "filter_barcode":
        if args.filename == None:
            exit(parser.parse_args(['filter_barcode','-h']))
        barcode_filter(args.filename, args.max_reads)
    
    if args.command == "filter_score":
        if args.filename == None:
            exit(parser.parse_args(['report','-h']))
        if args.value:
            if args.interval_max == None:
                interval_max = 100000
                if args.interval_min == 0:
                    exit(print("Please insert a percentile interval"))
                if args.interval_min < 0:
                    exit(print("Invalid percentile interval"))
            else:
                if args.interval_min > args.interval_max:
                    exit(print("Error! Min > Max"))
                interval_max = args.interval_max
            interval = [args.interval_min, interval_max]
            phred_filter_val(args.filename, args.max_reads, interval)
        else:
            if args.interval_max == None:
                interval_max = 100
                if args.interval_min == 0:
                    exit(print("Please insert a percentile interval"))
                if args.interval_min < 0:
                    exit(print("Invalid percentile interval"))
            else:
                if args.interval_max > 100:
                    exit(print("Invalid percentile interval"))
                if args.interval_min == 0 and args.interval_max == 100:
                    exit(print("Please insert a percentile interval"))
                if args.interval_min > args.interval_max:
                    exit(print("Error! Min > Max"))
                interval_max = args.interval_max
            interval = [args.interval_min, args.interval_max]
            phred_filter(args.filename, args.max_reads, interval)

    if args.command == "filter_length":
        if args.filename == None:
            exit(parser.parse_args(['report','-h']))
        if args.value:
            if args.interval_max == None:
                interval_max = 100000
                if args.interval_min == 0:
                    exit(print("Please insert a percentile interval"))
                if args.interval_min < 0:
                    exit(print("Invalid percentile interval"))
            else:
                if args.interval_min > args.interval_max:
                    exit(print("Error! Min > Max"))
                interval_max = args.interval_max
            interval = [args.interval_min, interval_max]
            length_filter_val(args.filename, args.max_reads, interval)
        else:
            if args.interval_max == None:
                interval_max = 100
                if args.interval_min == 0:
                    exit(print("Please insert a percentile interval"))
                if args.interval_min < 0:
                    exit(print("Invalid percentile interval"))
            else:
                if args.interval_max > 100:
                    exit(print("Invalid percentile interval"))
                if args.interval_min == 0 and args.interval_max == 100:
                    exit(print("Please insert a percentile interval"))
                if args.interval_min > args.interval_max:
                    exit(print("Error! Min > Max"))
                interval_max = args.interval_max
            interval = [args.interval_min, interval_max]
            length_filter(args.filename, args.max_reads, interval)

    print("\nTotal runtime : {time_used:.4f} seconds".format(time_used = (time.time() - start_time)))


if __name__ == "__main__":
    test()