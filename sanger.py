import os
import zipfile
import math
from Bio import SeqIO

def extractZip(zip_file: str) -> str:
    extract_to = os.path.join(os.path.dirname(zip_file), f"_Extracted_{os.path.basename(zip_file)[:-4]}")
    os.makedirs(extract_to)
    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
        zip_ref.extractall(extract_to)
    print(f"Extracted files to {extract_to}\n")

    return extract_to

def findAb1Files(directory: str) -> list:
    ab1_files = []
    for filename in os.listdir(directory):
        file = os.path.join(directory, filename)
        if file.endswith(".ab1") and os.path.isfile(file):
            ab1_files.append(file)
    
    return ab1_files

def convertAb1Fastq(ab1_files: list, output_folder: tuple[str, str]) -> str:
    if not ab1_files:
        print(f"No .ab1 files found")
        return "empty"
    
    output = os.path.join(output_folder[0], "_" + output_folder[1] + "_FASTQ") 
    if not os.path.exists(output):
        os.makedirs(output)
        print(f"Created output folder: {output}\n")

    for ab1_file in ab1_files:
        fastq_file = os.path.join(output, 
                                f"{os.path.splitext(os.path.basename(ab1_file))[0]}.fastq")
        
        print(f"Converting {ab1_file} to {fastq_file}\n")
        
        try:
            with open(ab1_file, "rb") as input_handle, open(fastq_file, "w") as output_handle:
                records = SeqIO.parse(input_handle, "abi")
                SeqIO.write(records, output_handle, "fastq")
            print(f"Successfully converted {ab1_file} to {fastq_file}\n")
        except Exception as e:
            print(f"Failed to convert {ab1_file}: {e}")
    
    return output

def asciiToPhred(ascii: chr) -> int:
    return ord(ascii) - 33

def avgQuality(quality_string: str):
    count = 0
    probs = []
    for q in quality_string:
        probs.append(10 ** (asciiToPhred(q) / -10))
        if asciiToPhred(q) < 20:
            count += 1
    avg_probs = sum(probs) / len(probs)
    return (-10 * math.log10(avg_probs), count) # Returns average phred score

def calcFastqQuality(fastq: str):
    with open(fastq, "r") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break  # Exit loop when no more lines are available
            sequence = f.readline().strip()
            f.readline()  # Skip the '+' line
            try:
                qual = f.readline().strip()[30:801]  # Adjusted for quality scores range
            except Exception as e:
                qual = f.readline().strip()  # Fallback if range is incorrect
            avg_qual, count = avgQuality(qual)

    header = os.path.basename(fastq)

    return (header, sequence, qual, avg_qual, count)

