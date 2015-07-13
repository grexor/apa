import apa
import pybio
import glob
import os
import sys
import gzip
import glob
import struct
import time

# parameters
min_read_length = 10
tail_window = 20
A_percentage = 0.85

def process_lib_ok(library_id, force=False):
    
    start_time = time.time()
    
    # call check_unique_id to see if we can use the original file id
    # assert(check_unique_id(library_id)==True)
    
    status_folder = os.path.join(apa.path.data_folder, library_id, "status")
    
    if os.path.exists(os.path.join(status_folder, "extract.status")) and force==False:
        print library_id, ": extract already processed or processing at the moment"
        return False
    
    if not os.path.exists(status_folder):
        os.makedirs(status_folder)
    f_log = open(os.path.join(status_folder, "extract.log"), "wt")
    f_status = open(os.path.join(status_folder, "extract.status"), "wt")

    lib = apa.annotation.libs[library_id]
    demulti_codes = {}
    fastq_files = {}
    for exp_id, exp_data in lib.experiments.items():
        assert(demulti_codes.get(exp_data["dcode"], None)==None) # only allow unique demulti codes
        demulti_codes[exp_data["dcode"]] = exp_id
        experiment_folder = os.path.join(apa.path.data_folder, library_id, "e%s" % exp_id)
        if not os.path.exists(experiment_folder):
            os.makedirs(experiment_folder)
        fastq_files[exp_id] = gzip.open(os.path.join(apa.path.data_folder, library_id, "e%s" % exp_id, "%s_e%s.fastq.gz" % (library_id, exp_id)), "w")
    aremoved_file = open(os.path.join(apa.path.data_folder, library_id, "%s.aremoved.bin" % (library_id)), "wb")    
    rnd_file = open(os.path.join(apa.path.data_folder, library_id, "%s.rnd.bin" % (library_id)), "wb")
    
    # take care of unmatched reads
    experiment_folder = os.path.join(apa.path.data_folder, library_id, "unmatched")
    if not os.path.exists(experiment_folder):
        os.makedirs(experiment_folder)
    fastq_files["unmatched"] = gzip.open(os.path.join(apa.path.data_folder, library_id, "unmatched", "%s_unmatched.fastq.gz" % (library_id)), "w")
    
    db_random = {}
    read_number = 1 # number reads, starting with 1
    rnd_code_number = 1 # number random barcodes, starting with 1
    unmatched_stats = {}
    lib_stats = {"all_reads" : 0}
    
    for fastq_filename in lib.fastq_files: # iterate over all fastq files in library
        print fastq_filename
        f = pybio.data.Fastq(os.path.join(apa.path.data_folder, library_id, fastq_filename))
        while f.read():
            lib_stats["all_reads"] += 1
            read_id = f.id.rstrip("/2")
            sequence = f.sequence.rstrip("A")
            cutat = len(sequence)
            new_cutat = remove_tail(f.sequence)
            if new_cutat < cutat:
                sequence = sequence[:new_cutat].rstrip("A")
            # final read processing
            sequence = sequence.rstrip("A")
            quality = f.quality[:len(sequence)]

            if len(sequence)<min_read_length: # do not process too short reads
                lib_stats["too_short"] = lib_stats.get("too_short", 0) + 1
                continue
            
            dcode = read_id.split("#")[-1][:lib.dcode_len] # experiment demulti code
            rnd_code = read_id.split("#")[-1][lib.dcode_len:] + sequence[:6]  # now: rnd_code #N from first_5 definition and first 6 nt of read
            
            # for each rnd_code of the read, check if you already have it in the database
            # if not, assign a new integer id
            rnd_code_id, rnd_code_freq = db_random.get(rnd_code, (None, None))
            if rnd_code_id==None:
                db_random[rnd_code] = (rnd_code_number, 1)
                rnd_code_id = rnd_code_number
                rnd_code_number += 1
            else:
                rnd_code_freq += 1
                db_random[rnd_code] = (rnd_code_id, rnd_code_freq)
            rnd_file.seek(read_number*4)
            rnd_file.write(struct.pack("I", rnd_code_id))
                
            exp_id = demulti_codes.get(dcode, "unmatched")
            lib_stats[exp_id] = lib_stats.get(exp_id, 0) + 1
            if exp_id=="unmatched":
                unmatched_stats[dcode] = unmatched_stats.get(dcode, 0) + 1
                          
            tail_len = len(f.sequence) - len(sequence)        
            if tail_len>0:
                aremoved_file.seek(read_number)
                aremoved_file.write(struct.pack("B", tail_len))
            fastq_files[exp_id].write("@%s\n%s\n+\n%s\n" % (read_number, sequence, quality))
    
            if f.count%100000==0:
                print library_id, ": EXTRACT : %.1fM reads" % (f.count/float(1000000)), ", unique random barcodes = %.1fM" % (rnd_code_number/float(1000000))
                f_log.write(library_id + ": EXTRACT : %.1fM reads" % (f.count/float(1000000)) + ", unique random barcodes = %.1fM\n" % (rnd_code_number/float(1000000)))
            
            #if read_number>300000:
            #    break
                
            read_number += 1 # next read
       
    # close files
    aremoved_file.close()
    rnd_file.close()

    for f in fastq_files.values():
        f.close()
    
    # write down random barcodes and their ids (just for readability)
    rnd_txt_file = open(os.path.join(apa.path.data_folder, library_id, "%s.rnd.txt" % (library_id)), "wt")
    rnd_txt_file.write("freq\trandom_code\trandom_code_id\n")
    L = [(rnd_code_freq, rnd_code, rnd_code_id) for rnd_code, (rnd_code_id, rnd_code_freq) in db_random.items()]
    L.sort(reverse=True)
    for (rnd_code_freq, rnd_code, rnd_code_id) in L:
        rnd_txt_file.write("%s\t%s\t%s\n" % (rnd_code_freq, rnd_code, rnd_code_id))
    rnd_txt_file.close()    

    # print out unmatched code stats
    L = [(val, k) for k, val in unmatched_stats.items()]
    L.sort(reverse=True)
    f = open(os.path.join(apa.path.data_folder, library_id, "unmatched", "demulti_codes.txt"), "wt")
    f.write("# codes not matching any annotated experiment code\n")
    f.write("code\tfreq\n")
    for (val, k) in L:
        f.write("%s\t%s\n" % (k, val))
    f.close()

    # stats (library)
    L = [(number, exp_id) for exp_id, number in lib_stats.items()]
    L.sort(reverse=True)
    f = open(os.path.join(apa.path.data_folder, library_id, "%s.stats.txt" % library_id), "wt")
    f.write("exp_id\tnum_reads\n")
    for (val, k) in L:
        f.write("%s\t%s\n" % (k, val))
    f.close()

    f_log.close()
    stop_time = time.time()
    f_status.write(library_id + " : time=%sm\n" % (int(stop_time-start_time)/60))
    f_status.close()


def remove_tail(sequence):
    # 2. search for A rich region
    for i in range(10, len(sequence)):
        tail = sequence[i:i+tail_window]
        tail_len = len(tail)
        num_A = tail.count("A")
        if tail.startswith("AAAAAAAAAAAAAAA"): # 15A? trim here
            return i
        if not tail.startswith("A"):
            continue
        if len(tail)>1 and not tail.startswith("AA"):
            continue
        if len(tail)>2 and not tail.startswith("AAA"):
            continue
        percentage_A = num_A / float(tail_len)
        if percentage_A>=A_percentage:
            return i
    return N

def process_lib(library_id, force=False):
    
    start_time = time.time()
    
    # call check_unique_id to see if we can use the original file id
    # assert(check_unique_id(library_id)==True)
    
    status_folder = os.path.join(apa.path.data_folder, library_id, "status")
    
    if os.path.exists(os.path.join(status_folder, "extract.status")) and force==False:
        print library_id, ": extract already processed or processing at the moment"
        return False
    
    if not os.path.exists(status_folder):
        os.makedirs(status_folder)
    f_log = open(os.path.join(status_folder, "extract.log"), "wt")
    f_status = open(os.path.join(status_folder, "extract.status"), "wt")

    lib = apa.annotation.libs[library_id]
    demulti_codes = {}
    fastq_files = {}
    for exp_id, exp_data in lib.experiments.items():
        assert(demulti_codes.get(exp_data["dcode"], None)==None) # only allow unique demulti codes
        demulti_codes[exp_data["dcode"]] = exp_id
        experiment_folder = os.path.join(apa.path.data_folder, library_id, "e%s" % exp_id)
        if not os.path.exists(experiment_folder):
            os.makedirs(experiment_folder)
        fastq_files[exp_id] = gzip.open(os.path.join(apa.path.data_folder, library_id, "e%s" % exp_id, "%s_e%s.fastq.gz" % (library_id, exp_id)), "w")
    rnd_file = open(os.path.join(apa.path.data_folder, library_id, "%s.rnd.bin" % (library_id)), "wb")
    
    # take care of unmatched reads
    experiment_folder = os.path.join(apa.path.data_folder, library_id, "unmatched")
    if not os.path.exists(experiment_folder):
        os.makedirs(experiment_folder)
    fastq_files["unmatched"] = gzip.open(os.path.join(apa.path.data_folder, library_id, "unmatched", "%s_unmatched.fastq.gz" % (library_id)), "w")
    
    db_random = {}
    read_number = 0
    rnd_code_number = 1 # number random barcodes, starting with 1
    unmatched_stats = {}
    lib_stats = {"all_reads" : 0}
    
    for fastq_filename in lib.fastq_files: # iterate over all fastq files in library
        print fastq_filename
        f = pybio.data.Fastq(os.path.join(apa.path.data_folder, library_id, fastq_filename))
        while f.read():
            read_number += 1 # id of current read (starts with 1)
            lib_stats["all_reads"] += 1
            read_id = f.id.rstrip("/2")
            sequence = f.sequence
            quality = f.quality

            if len(sequence)<min_read_length: # do not process too short reads
                lib_stats["too_short"] = lib_stats.get("too_short", 0) + 1
                continue
            
            dcode = read_id.split("#")[-1][:lib.dcode_len] # experiment demulti code
            rnd_code = read_id.split("#")[-1][lib.dcode_len:] + sequence[:6]  # now: rnd_code #N from first_5 definition and first 6 nt of read
            
            # for each rnd_code of the read, check if you already have it in the database
            # if not, assign a new integer id
            rnd_code_id, rnd_code_freq = db_random.get(rnd_code, (None, None))
            if rnd_code_id==None:
                db_random[rnd_code] = (rnd_code_number, 1)
                rnd_code_id = rnd_code_number
                rnd_code_number += 1
            else:
                rnd_code_freq += 1
                db_random[rnd_code] = (rnd_code_id, rnd_code_freq)
            rnd_file.seek(read_number*4)
            rnd_file.write(struct.pack("I", rnd_code_id))
                
            exp_id = demulti_codes.get(dcode, "unmatched")
            lib_stats[exp_id] = lib_stats.get(exp_id, 0) + 1
            if exp_id=="unmatched":
                unmatched_stats[dcode] = unmatched_stats.get(dcode, 0) + 1
            fastq_files[exp_id].write("@%s\n%s\n+\n%s\n" % (read_number, sequence, quality))
    
            if f.count%100000==0:
                print library_id, ": EXTRACT : %.1fM reads" % (f.count/float(1000000)), ", unique random barcodes = %.1fM" % (rnd_code_number/float(1000000))
                f_log.write(library_id + ": EXTRACT : %.1fM reads" % (f.count/float(1000000)) + ", unique random barcodes = %.1fM\n" % (rnd_code_number/float(1000000)))       
    # close files
    rnd_file.close()

    for f in fastq_files.values():
        f.close()
    
    # write down random barcodes and their ids (just for readability)
    rnd_txt_file = open(os.path.join(apa.path.data_folder, library_id, "%s.rnd.txt" % (library_id)), "wt")
    rnd_txt_file.write("freq\trandom_code\trandom_code_id\n")
    L = [(rnd_code_freq, rnd_code, rnd_code_id) for rnd_code, (rnd_code_id, rnd_code_freq) in db_random.items()]
    L.sort(reverse=True)
    for (rnd_code_freq, rnd_code, rnd_code_id) in L:
        rnd_txt_file.write("%s\t%s\t%s\n" % (rnd_code_freq, rnd_code, rnd_code_id))
    rnd_txt_file.close()    

    # print out unmatched code stats
    L = [(val, k) for k, val in unmatched_stats.items()]
    L.sort(reverse=True)
    f = open(os.path.join(apa.path.data_folder, library_id, "unmatched", "demulti_codes.txt"), "wt")
    f.write("# codes not matching any annotated experiment code\n")
    f.write("code\tfreq\n")
    for (val, k) in L:
        f.write("%s\t%s\n" % (k, val))
    f.close()

    # stats (library)
    L = [(number, exp_id) for exp_id, number in lib_stats.items()]
    L.sort(reverse=True)
    f = open(os.path.join(apa.path.data_folder, library_id, "%s.stats.txt" % library_id), "wt")
    f.write("exp_id\tnum_reads\n")
    for (val, k) in L:
        f.write("%s\t%s\n" % (k, val))
    f.close()

    f_log.close()
    stop_time = time.time()
    f_status.write(library_id + " : time=%sm\n" % (int(stop_time-start_time)/60))
    f_status.close()