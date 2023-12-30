'''
Description: calculation of Theta watterson, Theta Pi, and Tajima's D

Requirements:
- autosomes assumed as diploid by default
- chrY assumed as haploid by default
- whether to collapse the region into a relatively long sequences -> there are extremely short coding sequences in the genome, which is evolutionary constraint, but can provide evolutionary info.

reference:
- Watterson GA. On the number of segregating sites in genetical models without recombination. Theor Popul Biol. 1975;7(2):256–76. doi: 10.1016/0040-5809(75)90020-9.
'''
import time
import pandas as pd
import numpy as np
import math
from functools import reduce
#
import gzip
import multiprocessing
import argparse


def parseWindow(window_shift):
    """
    sliding windows parameters
    """
    windowsize = int(window_shift.split('@')[0])
    stepsize = int(window_shift.split('@')[1])
    overlapsize = windowsize - stepsize
    return windowsize, stepsize, overlapsize


def loadRegion(regionFile):
    """
    for only load regino file once
    """
    region = pd.read_csv(regionFile, 
                         sep='\s+', 
                         header=None, 
                         usecols=[0,1,2,3], 
                         names=['regionID','chr','start','end'],
                         low_memory=False)
    # skip expanding
    region['chr'] = region['chr'].astype(str)
    region.sort_values(by=['chr','start','end'],ascending=True,inplace=True)
    region.reset_index(inplace=True,drop=True)
    
    return region


def selectSampleIndex(inputFile, sampleFile):
    """
    parse the sample list for sample indexing and sfs estimation
    Notes: 
    - only open the start line of sfs, and generate the sample index
    - sample file only contains one column 
    """
    sampleList = []
    with gzip.open(inputFile, "rt") as input:
        for line in input:
            if line.startswith("#CHROM"):
                # headerLine = line.strip().split("\t")
                rawSample = line.strip().split("\t")[9:]
                sampleList = [x.split(".")[0] for x in rawSample]
                break
    # print(sampleList)
    
    # load sampleFile
    inputSample = []
    with open(sampleFile, "r") as input:
        for line in input:
            inputSample.append(line.strip())
    
    indexes = [sampleList.index(element) for element in inputSample if element in sampleList]
    return indexes


def parseGenotype(hapList, selected_samples):
    """
    split each line of genotype and convert to haplotype counts
    Notes: 
    - using phased vcf by default
    """
    new_hapList = [hapList[i] for i in selected_samples]
    separate_hapList = [allele for haps in new_hapList for allele in haps.split("|")]  # suitable for haploid calculation
    hapCounts = [separate_hapList.count("0"), separate_hapList.count("1")]
    # print(separate_hapList)
    return hapCounts


def ThetaEst(hapCounts, nSegregatingSites, nSeq):
    """
    calculation of nucleotide diversity based on 1) the number of segregating sites, 2) pairwise nucleotide difference

    - hapCountList: 
    - nSegregatingSites: the number of segregating sites
    - nSeq: total number of sequences, for example, for n diploids, there will be 2n seq, and n for n haploids.
    """
    pi = np.sum((hapCounts[:, 0] * hapCounts[:, 1])) * 1.0/(nSeq*(nSeq-1.0)/2.0)
    k = nSegregatingSites * 1.0 / reduce(lambda x,y: x+1.0/y, range(1,nSeq))

    return pi, k


def splitRegion(regionID, chromID, start, end, windowsize, stepsize, overlapsize):
    """
    split the region into multiple window for the theta calculation
    
    return 1) regionID; 2) chromID; 3) window start; 4) window end
    """
    # start_time = time.time()

    # 
    length = end - start + 1  # calculate region length
    bin_num = max(int(math.ceil((length - overlapsize)*1.0 / stepsize)),1)  # calculate the bin number
    ex_len = bin_num * stepsize + overlapsize  # if using overlapping windows, the original length -> the extended length
    ex_start = int(max(start-(ex_len-length)/2.0, 1.0)) # revise the start position to calculate the theta

    # construct split region array
    starts = np.arange(ex_start, ex_start + stepsize * bin_num, stepsize)
    ends = starts + windowsize - 1

    # stack the same elements
    splitRegion = np.column_stack(([regionID] * bin_num, [chromID] * bin_num, starts, ends))
    # print(f"stacked split region:\n{splitRegion}")
    # print(f"flatten split region:\n{splitRegion.flatten()}")
    
    # end_time = time.time()
    # elapsed_time = end_time - start_time
    # print(f"Elapsed time: {elapsed_time} seconds")
    
    return splitRegion


def flatRegions(regionRecords, windowsize, stepsize, overlapsize, targetRegion=False):
    """
    construct the region for the down-stream calculation
    """
    # sliding window parameters
    ## if window parameter is set to `target_region` -> using 5kb as window size
    if targetRegion:  # if the `window_shift` parameter is set to "target_region", calculate the theta for the corresponding region
        
        return regionRecords # should be converted to array
    else:
        # apply `splitRegion` to create region object for the down-stream calculation
        ## Notes: implement the dataframe column to save dataset may be memory-consuming
        
        totalLength = np.sum(regionRecords['end'] - regionRecords['start'])
        numRecords = max(int(math.ceil((totalLength - overlapsize)*1.0 / stepsize)), 1) 
        # print(f"numRecords: {numRecords}")

        # list with index as window index
        preList = np.empty((numRecords, 4), dtype=object)  # np is homogeneous -> only contains one type of data type
        # print(f"preList: {len(preList)}")
        # print(f"preList: {preList}")

        # 
        start_time = time.time()
        currentIndex = 0
        for _, row in regionRecords.iterrows():
            # retrieve split-up region
            tempRecords = splitRegion(row['regionID'], row['chr'], row['start'], row['end'], windowsize, stepsize, overlapsize)
            # print(tempRecords.shape)
            num_rows = tempRecords.shape[0]


            if currentIndex + num_rows > len(preList):
                # resize preList to accommodate more rows
                additional_rows_needed = currentIndex + num_rows - len(preList)
                preList = np.resize(preList, (len(preList) + additional_rows_needed, 4))

            
            preList[currentIndex:currentIndex + num_rows, :] = tempRecords
            # preList[currentIndex] = tempRecords
            # print(pd.DataFrame(tempValues, columns=['regionID', 'chr', 'start', 'end']))
            currentIndex += num_rows
            # print(currentIndex)
        
        # first test whether the input list is instance or not -> nested list can be directly convert to df, with single line
        # frameList = [pd.DataFrame(valueList, columns=['regionID', 'chr', 'start', 'end']) for valueList in preList]
        # split_region = np.concatenate(preList, axis=0)
        # print(split_region)

        # split_region = pd.DataFrame(split_region, columns=['regionID', 'chr', 'start', 'end'])
        # split_region.sort_values(by=['chr','start','end'], ascending=True, inplace=True)

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time} seconds")
        
        return preList


def process_variant_line(line, regions, selected_samples):
    """
    process one vcf line from the VCF file to find the region of the variant.
    
    :param line: each variants line
    :param regions: split_region
    """
    # extract variatn position, genotype information for the theta calculation
    split_line = line.strip("\n").split("\t")
    _, pos, genoTable = split_line[0], split_line[1], split_line[9:]  # ignore chromID by default

    for index, region in enumerate(regions):
        _, _, start, end = region
        if int(start) <= int(pos) <= int(end):
            # 
            # return index, genoTable

            # convert vcf to haplotype
            hapCounts = parseGenotype(genoTable, selected_samples)
            return index, hapCounts
    return None, None


def parallel_process_vcf(vcf_file, regions, samples, num_processes):
    """
    aggregate vcf lines parallelly
    """
    pool = multiprocessing.Pool(processes=num_processes)

    # initialize list to store results
    results = []
    with gzip.open(vcf_file, 'rt') as file:
        for line in file:
            # print(f"The dtype of parsed line is: type(line)")
            if line.startswith('#'):
                continue
            result = pool.apply_async(process_variant_line, (line, regions, samples)) # async mode multiprocessing
            results.append(result)
    pool.close()
    pool.join()

    # aggregate hapCounts
    regions_hapCounts_List = [[] for _ in regions]
    for result in results:
        region_index, hapCounts = result.get()  # retrieve region idnex and genotype data
        if region_index is not None:
            regions_hapCounts_List[region_index].append(np.array(hapCounts))
    
    return regions_hapCounts_List


def parallel_Theta(regions_hapCounts_List, nSeq):
    """
    calculate theta for each region
    """
    # loop version
    for index, hapCounts in enumerate(regions_hapCounts_List):
        new_hapCoutns = np.array(hapCounts)
        if new_hapCoutns.size > 0:
            # print(hapCounts)
            nSegregatingSites = np.count_nonzero(new_hapCoutns[:, 1])
            pi, k = ThetaEst(new_hapCoutns, nSegregatingSites, nSeq)

            return pi, k
            # print(f"{index}\t{split_region[index][0]}\t{split_region[index][2]}\t{split_region[index][3]}\t{pi}\t{k}")
        else:
            # print(f"{index}\t{split_region[index][0]}\t{split_region[index][2]}\t{split_region[index][3]}\t{0}\t{0}")
            return 0, 0



def main():
    
    parser = argparse.ArgumentParser(description='Theta_D_H.Est, https://github.com/Shuhua-Group/Theta_D_H.Est for more details')
    parser.add_argument("--gzvcf", type=str, required = True, \
                        help="phased.vcf.gz, format:GT (i.e., 0|1). able to deal with diploids and haploids, seperately.")
    parser.add_argument("--region", type=str, required = True, \
                        help="region.bed, variants in regions to be used, 4 columns: <region ID> <chrom ID> <start pos> <end pos>, no header line, tab or space sperated")
    parser.add_argument("--window_shift", type=str, required = False, default='target_region', \
                        help="windowsize@increment, for example, 50000@10000.")
    parser.add_argument("--out", type=str, required = False, default='out.txt', \
                        help="output file name. default: out.txt, will be automatically gzipped.")
    parser.add_argument("--samples", type=str, required = False, \
                        help="sample list")
    args = parser.parse_args()


    ### nSeq
    with gzip.open(args.gzvcf, "rt") as input:
        for line in input:
            if line.startswith("#"):
                continue
            else:
                split_line = line.strip("\n").split("\t")
                genoTable = split_line[9:]
                nSeq = len( [allele for geno in genoTable for allele in geno.split("|")] )
                break
    
    ### load original region
    windowsize, stepsize, overlapsize = parseWindow(args.window_shift)
    region = loadRegion(args.region)
    
    ### make samples 
    selected_samples = selectSampleIndex(args.gzvcf, args.samples)

    ### make regions
    split_region = flatRegions(region, windowsize, stepsize, overlapsize)
    print(split_region)

    ### separate variants to regions
    regions_hapCounts_List = parallel_process_vcf("/home/chenhongpu/opt/biosoft/biopipeline/population_genetics/theta/test/CHB.chr22.coding.vcf.gz", split_region, selected_samples, 8)

    
    output = open(args.out, "w")
    output.write(f"regionID\tchrom\twindow_start\twindwo_end\tThetaPi\tThetaW\n")
    ### calculate theta
    # loop version
    for index, hapCounts in enumerate(regions_hapCounts_List):
        new_hapCounts = np.array(hapCounts)
        # print(hapCounts)
        if new_hapCounts.size == 0:
            # print(hapCounts)
            regionID = split_region[index][0]
            chrom = split_region[index][1]
            region_start = split_region[index][2]
            region_end = split_region[index][3]
            pi = 0
            k = 0

            output.write(f"{regionID}\t{chrom}\t{region_start}\t{region_end}\t{str(pi)}\t{str(k)}\n")
        else:
            regionID = split_region[index][0]
            chrom = split_region[index][1]
            region_start = split_region[index][2]
            region_end = split_region[index][3]

            nSegregatingSites = np.count_nonzero(new_hapCounts[:, 1])
            pi, k = ThetaEst(new_hapCounts, nSegregatingSites, nSeq)

            output.write(f"{regionID}\t{chrom}\t{region_start}\t{region_end}\t{str(pi)}\t{str(k)}\n")

    output.close()






if __name__ == "__main__":
    # run
    main()

    # nSeq = 206 # -> calculate once
    # window_shift = "50000@50000"
    # windowsize, stepsize, overlapsize = parseWindow(window_shift)

    # # load region

    # # # Testing
    # # tmpRegion = region.iloc[1, ]
    # # print(tmpRegion)
    # # regionRecords = splitRegion(tmpRegion[0], tmpRegion[1], tmpRegion[2], tmpRegion[3], windowsize, stepsize, overlapsize)
    # # print(regionRecords)

    # split_region = flatRegions(region, windowsize, stepsize, overlapsize)
    # print(split_region)
    # print(len(split_region))

    # regions_hapCounts_List = parallel_process_vcf("/home/chenhongpu/opt/biosoft/biopipeline/population_genetics/theta/test/CHB.chr22.coding.vcf.gz", split_region, 8)
    # print(regions_hapCounts_List)

    # # parallel_Theta(regions_hapCounts_List)
    

    # # main("./test/CHB.chr22.coding.vcf.gz")