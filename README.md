# LARGE-PD_Phase2_Paper

`GP2 ❤️ Open Science 😍`

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15864760.svg)](https://doi.org/10.5281/zenodo.15864760)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Last Updated:** July 2025

This repository contains a detailed description of all analyses performed as part of the LARGE-PD Phase 2 manuscript. The Latin American Research Consortium on the Genetics of Parkinson’s Disease (LARGE-PD) is an international collaborative effort aimed at understanding the genetic architecture of Parkinson's disease in Latin American populations, with a particular focus on admixed and underrepresented groups.

The analyses documented here include:

🧬 Genotyping Quality Control (QC)

🌍 Global and Local Ancestry Inference

📊 Genome-Wide Association Studies (GWAS) and Admixture Mapping

🌐 Trans-ethnic Meta-analysis

# Genotyping calling
## Prerequisites 
This steps requires:
- Plink2 installed (https://www.cog-genomics.org/plink/2.0/)
- IAAP installed (https://emea.support.illumina.com/downloads/iaap-genotyping-cli.html)
- Python3
- GP2 cluster file (.egt). 
    - The egt used in our work was provided by Mike Nalls (mike@datatecnica.com)
- Illumina Sample Sheet 
    - The sample sheet is a CSV file with two header lines and one entry per line.

Example: sampleSheet.csv
```
[Data]
Sample_ID,SentrixBarcode_A,SentrixPosition_A,Path
ExternalID1,123456789,R10C11,/home/input/idat/123456789/
ExternalID2,123456789,R11C11,/home/input/idat/123456789/
ExternalID3,123456789,R12C11,/home/input/idat/123456789/
ExternalID4,123456710,R12C11,/home/input/idat/123456710/
```

## Execution
### Genotyping calling usign IAAP CLI

We downloaded the Illumina Array Analysis Platform (IAAP) Genotyping Command Line Interface (CLI) and used it to generate a PLINK PED file for each entry in the sample sheet using the following command:
```
iaap-cli gencall /home/input/calling/NeuroBooster_20042459_A4.bpm \
/home/input/calling/recluster_09092022.egt /home/output/genotyping/ \
-s /home/input/calling/sampleSheet.csv -p -t 8 -c 0.15 
```

We do not recommend using more than eight threads.
We set the GenCall score cutoff to 0.15 because some poorly performing variants were retained with the default cutoff of 0.2.

### Merging multiples PED files

After running the IAAP CLI, we obtained a single .map file (same prefix name as the .bpm) and a .ped file for each sample. The .ped filenames follow the pattern \<SentrixBarcode\>_\<SentrixPosition\>.ped.
We merged all .ped files using the script below:

```python=
fileOut = open("/home/output/genotyping/NeuroBooster_20042459_A4.ped", "w")
fileIn = open("/home/input/calling/sampleSampleSheet.csv")

header = True
for line in fileIn:
    if header:
        if line.startswith("Sample_ID"):
            header = False
    else:
        newID, barcode, position, path = line.strip().split(",")
        for i in range(1,15):
            ID = newID.replace(f"_B{i}", "") # Removing the batch from ID
        else:
            ped = open(f"{barcode}_{position}.ped")
            print(f"Adding {barcode}_{position}.ped")
            for linePed in ped:
                split = linePed.strip().split()
                fileOut.write(f"{ID} {ID}")
                for j in range(2,len(split)):
                    fileOut.write(f" {split[j]}")
                fileOut.write(f"\n")
            ped.close()
fileOut.close()
fileIn.close()
```
Note: LARGE-PD sample IDs contain batch information (e.g., _B1, _B2, etc.) in the Sample_ID column of the sample sheet. This script removes the batch suffix during processing.
    
## Converting to BED/BIM/FAM
    
After merging all samples, we convert to BED/BIM/FAM using plink2
    
```bash=
plink2 --make-pgen --out /home/output/genotyping/NB --pedmap /home/output/genotyping/NeuroBooster_20042459_A4 --sort-vars
```

## Removing the NeuroBooster Array problems
    
After extensive research, members of the Mata Lab dry lab team, led by Dr. Thiago Peixoto Leal, identified three types of issues in the array data by comparing multiple Illumina annotation tables: (i) LiftOver errors (variants with different chromosomes on hg37 and hg38), (ii) non-biallelic variants genotyped on the minus which one genotyped allele is the reference allele on + (causing error while merging/preparing to imputation) (“–”) strand, and
(iii) variants where the reference allele was not genotyped.

We compiled a list of all such variants (exclude.txt) and removed them from the dataset
Removing the NBA problems 

```
plink2 --pfile /home/output/genotyping/NB --chr 1-27 --exclude /home/input/calling/genotyping/exclude.txt --make-pgen --out /home/output/genotyping/NB_NoLift 
```

# Genotyping Quality Control
## Prerequisites 
- Plink 1.90
- Plink 2
- NAToRA (https://github.com/ldgh/NAToRA_Public)
- king
- Python3
    -NetworkX installed 

## Execution
At the Mata Lab, we do not use GenoTools for performing genetic quality control (QC) because some of its steps can lead to unnecessary sample loss, particularly due to the high levels of admixture present in Latin American populations.
    
The Mata Lab Genotyping QC can be find on Mata Lab Github (https://github.com/MataLabCCF/GWASQC)

```bash=
python3.8 main.py --plink1 /home/programs/plink --plink2  /home/programs/plink2 \
    -O /home/output/QCed/ -o LARGE -i /home/output/genotyping/NB_NoLift \
    -P /home/input/Demographic/County.txt -S -I /home/input/Demographic/Covar.txt \
    --NAToRA /home/programs/NAToRA_Public.py
```
    
We use NAToRA (https://github.com/ldgh/NAToRA_Public) to perform relationship control because it outperforms PLINK and KING in relatedness removal, offering the option to select the optimal result that guarantees minimal sample loss. For more details about the LARGE-PD QC, please check the GitHub

## Imputation

To send to TOPMed imputation server we have to flip the variants and fix the ref/alt switches. We use bcftools +fix-ref to prepare the data

```
/home/programs/plink2 --export vcf-iid bgz --out /home/output/QCed/Test \
--output-chr chr26 --pfile /home/output/QCed/LARGE_HWECase

/home/programs/bcftools +fixref /home/output/QCed/Test.vcf.gz-Oz -o /home/output/QCed/FixRef.vcf.gz -- -f /home/references/hg38/Homo_sapiens_assembly38.fasta -m flip -d
```

```python=
import os
for i in range(1,23):
    os.system(f"/home/peixott/beegfs/Programs/plink2 "
              f"--vcf FixRef.vcf.gz --output-chr chr26 "
              f"--chr {i} --recode vcf-iid bgz "
              f"--out LARGE_chr{i}_ToTOPMED")

```
We uploaded the LARGE_chr\<chrom>_ToTOPMED.vcf.gz. We later downloaded the imputed data, extracted and indexed the data using 


```
/home/programs/bcftools index -c -f chrom<chrom>.dose.vcf.gz
```

# Admixture inference
## Prerequisites 
- Plink 1.90
- Plink 2
- ADMIXTURE
    
## Execution
We performed global ancestry inference using supervised ADMIXTURE. As a reference panel, we used the dataset described by Shriner et al. (DOI: 10.1016/j.xhgg.2023.100235). The sample IDs used can be found on the Mata Lab GitHub repository:
https://github.com/MataLabCCF/GnomixPretrained/tree/main/Parentals

## Extract the reference samples
We downloaded the 1000 Genomes 30x dataset
(https://www.internationalgenome.org/data-portal/data-collection/30x-grch38)
and extracted the reference samples listed in AllRef.txt using PLINK2:
    
```
plink2 --vcf /home/references/1KGP/PGEN/1KGP30x_chrom<chrom>.vcf.gz \
--chr <chrom> --keep /home/references/AllRef.txt \
--make-pgen --out /home/references/1KGP/PGEN/OneThousand_chr<chrom>
```
Note: AllRef.txt contains the list of reference individuals from Shriner et al.
   
## Merge all chromosomes
```
plink2 --pmerge-list /home/references/1KGP/PGEN/ListOne.txt --make-pgen --out /home/references/1KGP/PGEN/OneThousand_All  
```

<details>
    <summary>ListOne.txt</summary>
/home/references/1KGP/PGEN/OneThousand_chr10 <br>
/home/references/1KGP/PGEN/OneThousand_chr11 <br>
/home/references/1KGP/PGEN/OneThousand_chr12 <br>
/home/references/1KGP/PGEN/OneThousand_chr13 <br>
/home/references/1KGP/PGEN/OneThousand_chr14 <br>
/home/references/1KGP/PGEN/OneThousand_chr15 <br>
/home/references/1KGP/PGEN/OneThousand_chr16 <br>
/home/references/1KGP/PGEN/OneThousand_chr17 <br>
/home/references/1KGP/PGEN/OneThousand_chr18 <br>
/home/references/1KGP/PGEN/OneThousand_chr19 <br>
/home/references/1KGP/PGEN/OneThousand_chr1 <br>
/home/references/1KGP/PGEN/OneThousand_chr20 <br>
/home/references/1KGP/PGEN/OneThousand_chr21 <br>
/home/references/1KGP/PGEN/OneThousand_chr22 <br>
/home/references/1KGP/PGEN/OneThousand_chr2 <br>
/home/references/1KGP/PGEN/OneThousand_chr3 <br>
/home/references/1KGP/PGEN/OneThousand_chr4 <br>
/home/references/1KGP/PGEN/OneThousand_chr5 <br>
/home/references/1KGP/PGEN/OneThousand_chr6 <br>
/home/references/1KGP/PGEN/OneThousand_chr7 <br>
/home/references/1KGP/PGEN/OneThousand_chr8 <br>
/home/references/1KGP/PGEN/OneThousand_chr9
</details>

# Supervised Admixture

With the merged reference file, we used the script from the Mata Lab repository (https://github.com/MataLabCCF/AdmixtureWithReference)
to generate the input files required to run ADMIXTURE:
    
```
python3.8 prepareAdmixture.py \
-i /home/output/QCed/LARGE_HWECase \
-I /home/references/1KGP/OneThousand_All \
-o /home/output/Admixture -F /home/references/Parentals \
-P /home/programs/plink2 -p /home/programs/plink -c 
```

```
cd /home/output/Admixture 
/home/programs/admixture --cv LD_Bfile1And2.bed 5 -j40 --supervised | tee log5_1.out 
```

We also used the file Common_Bfile2 generated by this pipeline. It was split by chromosome, converted to VCF, and submitted to the TOPMed Imputation Server for phasing. We later downloaded, extracted on /home/input/Phased/
    
# PCA
## Prerequisites 
- Plink 1.90
- Plink 2
- GCTA
    
## Execution
Our work used both supervised and unsupervised principal component analysis (PCA).
Both PCAs were generated using the publicly available code from the LARGE-PD GitHub repository:
https://github.com/MataLabCCF/ProjectedPCAAndModelSelection

The first PCA we generated was a supervised PCA, using the same reference panel employed in the supervised ADMIXTURE analysis. This PCA was used to identify ancestry outliers for analyses that include ancestry as a covariate in the association model.
    
```bash=
# Supervised PCA
python covarProjectedPCA.py \
   -A /home/output/QCed/LARGE_HWECase  \
   -t /home/input/Demographic/Covar.txt \ 
   -R /home/references/1KGP/OneThousand_All \
   -n LARGE -f /home/output/PCA/Projected \
   --gcta /home/programs/gcta64 \
   --selectModel /home/output/PCA/selectModel.R \
   --plink1 /home/programs/plink --plink2 /home/programs/plink2
```

We also generated an unsupervised PCA using all QCed samples. This PCA was used as covariates in the SAIGE GWAS. SAIGE incorporates a Genetic Relatedness Matrix (GRM) into its statistical model, which allows the inclusion of related individuals without introducing bias into the association results.

```bash=
# Unsupervised PCA with all samples
python covarProjectedPCA.py \
   -A /home/output/QCed/LARGE_HWECase   \
   -t /home/input/Demographic/Covar.txt \
   -n LARGE -f /home/output/PCA/NonProjected \
  --gcta /home/programs/gcta64 \
  --selectModel /home/output/PCA/selectModel.R \
  --plink1 /home/programs/plink --plink2 /home/programs/plink2
```

We also generated a PCA excluding genetic outliers (identified using the projected PCA) to be used in analyses that include ancestry in the statistical model, such as ATT, TRACTOR, and admixture mapping.

To achieve this, we first created a list of genetic outliers using the script bellow.


<details>
    <summary>outliers.py</summary>
    
```python=
import numpy as np

fileIn = open("AutosomalPCA_TargetPC.proj.eigenvec")
dictPCA = {}
sumStat = {}

for line in fileIn:
    split = line.strip().split()
    for i in range(2, len(split)):
        PC = f"PC{i-1}"
        if PC not in dictPCA:
            dictPCA[PC] = []
        dictPCA[PC].append(float(split[i]))
fileIn.close()

for PC in dictPCA:
    sumStat[PC] = {}
    sumStat[PC]["mean"] = np.mean(dictPCA[PC])
    sumStat[PC]["sd"] = np.std(dictPCA[PC])

toRemoveListPC1_10 = []
toRemoveListPCStepwise = []
fileIn = open("AutosomalPCA_TargetPC.proj.eigenvec")
fileOutToRemove = open("toRemovePC1-10.txt", "w")
fileOutToRemoveList = open("toRemoveExplanation.txt", "w")
for line in fileIn:
    split = line.strip().split()
    for i in range(2, len(split)):
        PC = f"PC{i-1}"
        PCInd = float(split[i])

        minRange = sumStat[PC]["mean"] - (3*sumStat[PC]["sd"])
        maxRange = sumStat[PC]["mean"] + (3*sumStat[PC]["sd"])

        if PCInd \< minRange or PCInd > maxRange:
            if i \<= 7:
                if split[1] not in toRemoveListPC1_10:
                    fileOutToRemove.write(f"0\t{split[1]}\n")
                toRemoveListPC1_10.append(f"0\t{split[1]}\n")
                fileOutToRemoveList.write(f"{split[1]}\t{PC}\n")
fileOutToRemove.close()
```
                     
</details>

We performed PCA without reference with all samples    

```bash=
# Unsupervised PCA with all samples
python covarProjectedPCA.py \
   -A /home/output/QCed/LARGE_HWECase   \
   -t /home/input/Demographic/Covar.txt \
   -n LARGE -f /home/output/PCA/NonProjected \
  --gcta /home/programs/gcta64 \
  --selectModel /home/output/PCA/selectModel.R \
  --plink1 /home/programs/plink --plink2 /home/programs/plink2
```

We also performed PCA without reference without outliers

```bash=
#Unsupervised non outlier
python covarProjectedPCA.py \
   -A /home/output/QCed/LARGE_HWECase   \
   -t /home/input/Demographic/Covar.txt \
   -n LARGE -f /home/output/PCA/NonProjectedNonOutlier \
  --gcta /home/programs/gcta64 --remove toRemovePC1-10.txt \
  --selectModel /home/output/PCA/selectModel.R \
  --plink1 plink --plink2 plink2 
```

We merged /home/output/PCA/NonProjected/AutosomalPCA_PCs.eigenvec with /input/Demographic/Covar.txt and created a new file named /input/Demographic/CovarToSAIGE.txt. We merged /home/output/PCA/NonProjectedNonOutlier/AutosomalPCA_PCs.eigenvec with /input/Demographic/Covar.txt and created a new file named /input/Demographic/CovarToAK.txt

# Local Ancestry
## Prerequistes
- Imputed data to GWAS from TOPMed
- Phased reference data from TOPMed
- bcftools

## Execution
The first step is extracting the TYPED variants from the imputed dataset
                     
```python=
import os

for i in range(1,23):
    os.system(f"/home/programs/bcftools view -i \"INFO/TYPED=1\"" 
              f"/home/Imputed/chr{i}.dose.vcf.gz -Oz -o "
              f"/home/output/Typed/CHR{i}_OnlyTyped.vcf.gz")      
    os.system(f"/home/programs/bcftools index /home/output/Typed/CHR{i}_OnlyTyped.vcf.gz" 
```

After this, we build the reference panel by merging the LARGE-PD non admixed Native American samples with the phased reference dataset using the build refTarget.py

<details>
    <summary>refTarget.py</summary>
            
```python= 
import argparse
import gzip
import os


def execute(command):
    print(f"{command}")
    os.system(command)
    print(f"===========================================================================================")

def buildTargetFile(VCF, chrom, listFolder):
    VCF = gzip.open(VCF.replace("&", str(chrom)))
    listOut = open(f"{listFolder}/TargetList.txt", "w")

    for line in VCF:
        line = line.decode("utf-8")
        if line.startswith("#CHROM"):
            header = line.strip().split()

            for i in range(9, len(header)):
                listOut.write(f"{header[i]}\n")
            break
    VCF.close()
    listOut.close()

    return f"{listFolder}/TargetList.txt"

def buildReferenceFile(correspondence, listFolder):
    correspondenceFile = open(correspondence)
    listOut = open(f"{listFolder}/RefList.txt", "w")

    for line in correspondenceFile:
        listOut.write(f"{line.strip().split()[0]}\n")
    correspondenceFile.close()
    listOut.close()

    return f"{listFolder}/RefList.txt"

def mergeRefAndTarget(target, reference, chrom, mergedFolder, bcftools):
    refWithChr = reference.replace("&", str(chrom))
    targetWithChr = target.replace("&", str(chrom))

    #Change the variant ID
    execute(f"{bcftools} annotate -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' {targetWithChr} -Oz "
            f"-o {mergedFolder}/TARGET_Annotated_chr{chrom}.vcf.gz")

    execute(f"{bcftools} annotate -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' {refWithChr} -Oz "
            f"-o {mergedFolder}/REFERENCE_Annotated_chr{chrom}.vcf.gz")

    execute(f"{bcftools} index {mergedFolder}/TARGET_Annotated_chr{chrom}.vcf.gz")
    execute(f"{bcftools} index {mergedFolder}/REFERENCE_Annotated_chr{chrom}.vcf.gz")

    #Get list of variants in common
    execute(f"{bcftools} isec -n=2 {mergedFolder}/TARGET_Annotated_chr{chrom}.vcf.gz "
            f"{mergedFolder}/REFERENCE_Annotated_chr{chrom}.vcf.gz -o {mergedFolder}/CommonVariants_chr{chrom}")

    #Merge
    execute(f"{bcftools} merge -R {mergedFolder}/CommonVariants_chr{chrom} "
            f"{mergedFolder}/TARGET_Annotated_chr{chrom}.vcf.gz {mergedFolder}/REFERENCE_Annotated_chr{chrom}.vcf.gz "
            f"-Oz -o {mergedFolder}/Merged_chr{chrom}.vcf.gz")
    execute(f"{bcftools} index {mergedFolder}/Merged_chr{chrom}.vcf.gz")

    return f"{mergedFolder}/Merged_chr{chrom}.vcf.gz"


def filterVCF(VCF, fileToFilter, bcftools, outFolder, name):
    execute(f"{bcftools} view {VCF} -S {fileToFilter} -Oz -o {outFolder}/{name}_chr{chrom}.vcf.gz")
    execute(f"{bcftools} index {outFolder}/{name}_chr{chrom}.vcf.gz")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build panel based on TOPMed imputation server')

    requiredIn = parser.add_argument_group("Input files arguments")
    requiredIn.add_argument('-r', '--reference',
                                 help='Reference database phased with TOPMed with chromosome replaced by &', required=True)
    requiredIn.add_argument('-t', '--target',
                                 help='Target database from TOPMed with chromosome replaced by &', required=True)
    requiredIn.add_argument('-c', '--correspondence',
                                 help='Correspondence list. It may contain samples from reference and target',
                                 required=True)
    requiredIn.add_argument('-b', '--begin', help='First chromosome (default = 1)', default=1, type=int)
    requiredIn.add_argument('-e', '--end', help='Last chromosome  (default = 22)', default=22, type=int)

    requiredOut = parser.add_argument_group("Output arguments")
    requiredOut.add_argument('-O', '--outputFolder', help='Name of output folder', required=True)

    programs = parser.add_argument_group("Programs")
    programs.add_argument('--bcftools', help='bacftools path (default: bcftools)',
                             required=False, default = "bcftools")


    args = parser.parse_args()

    os.system(f"mkdir {args.outputFolder}")

    listFolder = f"{args.outputFolder}/Lists/"
    os.system(f"mkdir {args.outputFolder}/Lists/")

    os.system(f"mkdir {args.outputFolder}/VCF/")
    mergedFolder = f"{args.outputFolder}/VCF/Merged/"
    os.system(f"mkdir {args.outputFolder}/VCF/Merged/")

    refFolder = f"{args.outputFolder}/VCF/Ref/"
    os.system(f"mkdir {refFolder}")

    targetFolder = f"{args.outputFolder}/VCF/Target/"
    os.system(f"mkdir {targetFolder}")

    targetFile = buildTargetFile(args.target, args.begin, listFolder)
    refFile = buildReferenceFile(args.correspondence, listFolder)
    for chrom in range(args.begin, args.end+1):
        VCFMerged = mergeRefAndTarget(args.target, args.reference, chrom, mergedFolder, args.bcftools)
        refVCF = filterVCF(VCFMerged, refFile, args.bcftools, refFolder, "REFERENCE")
        targetVCF = filterVCF(VCFMerged, targetFile, args.bcftools, targetFolder, "TARGET")
```
    
</details>

Example of command line:
```python=
import os
for i in range(1,23):
    os.system(f"python buildRefTarget.py "
                f"-r /home/input/Phased/chr{i}.reheaderphased.vcf.gz "
                f"--target /home/output/LA/Typed/CHR{i}_OnlyTyped.vcf.gz "
                f"-b {i} -e {i} "
                f"--correspondence /home/output/LA/Correspondence.txt "
                f"-O /home/output/LA/Panels")
                                    
```
Correspondece is a two columns file, the first with Individual ID and the second is the ancestry.
                     
To run Gnomix after generating the panels (ref and target):
                     
```python=
for count in range(1,23):
    os.system(f"mkdir /home/output/LA/Inference/LPD_chrom{count}")
    os.system(f"python3.8 /home/programs/gnomix/gnomix.py "
        f"/home/output/LA/Panels/VCF/Target/TARGET_chr{count}.vcf.gz "
        f"/home/output/LA/Inference/LPD_chrom{count} {count} False "
        f"/home/input/LA/G_chr{count}.b38.gmap "
        f"/home/output/LA/Panels/VCF/Ref/REFERENCE_chr{count}.vcf.gz "
        f"/home/output/LA/Correspondence.txt  "
        f"/home/input/configBest.yaml")                 
```

You can find the configBest.yaml and the pre-trained models at https://github.com/MataLabCCF/GnomixPretrained

# Associations
## Prerequisites
- SAIGE GWAS installed (https://saigegit.github.io/SAIGE-doc/)
- GENESIS installed (https://rdrr.io/bioc/GENESIS/)
- king installed
- Admix-kit installed (https://github.com/KangchengHou/admix-kit)
- Python3

## SAIGE GWAS
### Fit model
#### Convert PFILE to BED/BIM/FAM
```
mkdir /home/output/Association/PLINK
/home/programs/plink2 --pfile /home/output/QCed/LARGE_HWECase --chr 1-22 \
    --output-chr 26 --make-bed --out /home/output/Association/PLINK/LPD_Phase2
```

#### Generate the GRM
```
Rscript createSparseGRM.R \
    --plinkFile /home/output/Association/SAIGE/PLINK/LPD_Phase2 \
    --outputPrefix=/home/output/Association/SAIGE/GRM \
    --relatednessCutoff=0.088 --nThreads=5 \
    --numRandomMarkerforSparseKin 500000
```

#### Fit NULL model
```
Rscript step1_fitNULLGLMM.R \
    --plinkFile=/home/output/Association/SAIGE/PLINK/LPD_Phase2\
    --phenoFile=/home/input/Demographic/covarToSAIGE.txt \
    --phenoCol=DISEASE \
    --covarColList=AGE,SEX,GCTA_PC1,GCTA_PC2,GCTA_PC3,GCTA_PC4,GCTA_PC5,GCTA_PC6,GCTA_PC7,GCTA_PC8,GCTA_PC9,GCTA_PC10 \
    --qCovarColList=SEX --sampleIDColinphenoFile=IID \
    --isCateVarianceRatio=TRUE --traitType=binary \
    --outputPrefix=/home/output/Association/SAIGE/Model_LPD_Phase2_PC1-10 \
    --nThreads=20 --IsOverwriteVarianceRatioFile=TRUE --useSparseGRMtoFitNULL=TRUE \
    --sparseGRMFile /home/output/Association/SAIGE/GRM_relatednessCutoff_0.088_500000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile /home/output/Association/SAIGE/GRM_relatednessCutoff_0.088_500000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
``` 

#### Association test

Here we present the line for chromosome 22.

```
mkdir /home/output/Association/SAIGE/Output/
Rscript step2_SPAtests.R \
--vcfFile=/home/input/Imputed/chr<chrom>.dose.vcf.gz \
--vcfFileIndex=/home/input/Imputed/chr<chrom>.dose.vcf.gz.csi \
--vcfField=DS --SAIGEOutputFile=/home/output/Association/SAIGE/Output/Result_LPD_Phase2_Chrom<chrom> \
--minMAF=0 --minMAC=20 \
--GMMATmodelFile=/home/output/Association/SAIGE/Model_LPD_Phase2_PC1-10.rda \
--varianceRatioFile=/home/output/Association/SAIGE/Model_LPD_Phase2_PC1-10.varianceRatio.txt \
--sparseGRMFile /home/output/Association/SAIGE/GRM_relatednessCutoff_0.088_500000_randomMarkersUsed.sparseGRM.mtx \
--sparseGRMSampleIDFile /home/output/Association/SAIGE/GRM_relatednessCutoff_0.088_500000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
--LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_output_moreDetails=TRUE \
--AlleleOrder=ref-first --chrom=chr22
```

## Admix-kit analysis

TRACTOR and ATT was performed using the Admix-kit

### Prepare the data

```python=
import os

#Setting folders
plink2 = f"/home/programs/plink2"
pgenFolder = f"/home/output/Association/AK/PGEN/"
infoFolder = f"/home/input/Demographic/covarToAK.txt "
outFolder = f"/home/output/Association/AK/Results/"
imputationFolder = f"/home/input/Imputed/"
mspFolder = f"/home/output/LA/Inference/LPD_chrom{chrom}/query_results.msp"
#removeToAK.txtis the NAToRA relationship file + outliers 
relationship = "/home/output/Association/AK/removeToAK.txt"

os.system(f"mkdir {outFolder}\n")
os.system(f"mkdir {pgenFolder}\n")
os.system(f"{plink2} --vcf {imputationFolder}/chr{i}.dose.vcf.gz "
            f"--remove {relationship} --make-pgen "
            f"--out {pgenFolder}/LPD_{i} "
            f" --maf 0.005 --max-alleles 2 --rm-dup exclude-all "
            F"--snps-only --set-missing-var-ids @:#:\\$r:\\$a\n")
os.system(f"python3.8 newMSP.py -p {pgenFolder}/LPG_{i}.psam "
            f"-m {mspFolder} -o {pgenFolder}/MSP_{i}.msp")
os.system(f"admix lanc-convert {pgenFolder}/LPG_{i} "
            f"--rfmix {pgenFolder}/MSP_{i}.msp --out "
            f"{pgenFolder}/LPG_{i}.lanc\n")
```


The newMSP is a script that build a msp with the samples that are listed on psam file.

    
<details>
    <summary>newMSP.py</summary>
            
```python=
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='newMSP')

    requiredGeneral = parser.add_argument_group("Required arguments")
    requiredGeneral.add_argument('-p', '--psam', help='PSAM file', required=True)
    requiredGeneral.add_argument('-m', '--msp', help='MSP file', required=True)
    requiredGeneral.add_argument('-o', '--output', help='Name of output folder', required=True)
    
    args = parser.parse_args()
    
    file = open(f"{args.psam}")
    header = True
    
    sampleList = []
    
    for line in file:
        if header:
            header = False
        else:
            split = line.strip().split()
            sampleList.append(split[0])
    
    file.close()
    fileIn = open(args.msp)
    fileOut = open(args.output, "w")
    
    dictIDs = {}
    
    header = True
    for line in fileIn:
        if header:
            if "#chm" in line:
                IDs = line.strip().split("\t")
                for i in range(6, len(IDs),2):
                    IDComponent = IDs[i].split(".")
                    IDNoHap = IDComponent[0]
                    for j in range(1, len(IDComponent)-1):
                        IDNoHap = IDNoHap+f".{IDComponent[j]}"
#                    print(f"{IDs[i]} - {IDNoHap}")
                    dictIDs[IDNoHap] = i
                    
                #Header of new MSP
                
                fileOut.write(f"{IDs[0]}")
                for i in range(1,6):
                    fileOut.write(f"\t{IDs[i]}")
                    
                for sample in sampleList:
                    fileOut.write(f"\t{sample}.0\t{sample}.1")
                fileOut.write("\n")
                
                header = False
            else:
                fileOut.write(line)
        else:
            data = line.strip().split("\t")
            fileOut.write(f"{data[0].replace('chr', '')}")
            for i in range(1,6):
                fileOut.write(f"\t{data[i]}")
            for sample in sampleList:
                index = dictIDs[sample]
                fileOut.write(f"\t{data[index]}\t{data[index+1]}")
            fileOut.write("\n")
            
    fileOut.close()

for ID in dictIDs:
	print(f"{ID}: {dictIDs[ID]}")

```
</details>
    
### Run TRACTOR and ATT

```python=
import os

for i in range(22,0, -1):
    for method in ["ADM", "ATT", "TRACTOR"]:
        plink2 = f"/home/programs/plink2"
        pgenFolder = f"/home/output/Association/AK/PGEN/"
        infoFolder = f"/home/input/Demographic/covarToAK.txt "
        outFolder = f"/home/output/Association/AK/Results/"
        os.system(f"admix assoc --pfile {pgenFolder}/LPD_{i} "
                    f"--family binary --pheno {infoFolder}/Covar.txt "
                    f"--method {method} --quantile-normalize True "
                    f"--out {outFolder}/LPD_{method}_{i}")
```

It is mandatory IID on the first column, and phenotype on the second collumn on Covar.txt. All other collumns will be included on the model as covariate.

## Admixture mapping

Mata lab members built a code to run Admxiture Mapping using GENESIS. You can find it on https://github.com/MataLabCCF/GenesisAM

```
python GENESIS_AM.py \
    --msp /home/output/LA/Inference/LPD_chrom\*/query_results.msp \
    -b 1 -e 22 -p /output/SAIGE/PLINK/LPD_Phase2 \
    --plink1 /home/programs/plink -R Rscript -K /home/programs/king \
    -o Phase2 -O /home/output/AM/Phase2 \
    --covar /home/input/Demographic/covarToAK.txt

```

Note: The author of the python scripts (Thiago Peixoto Leal) apologizes for not using pandas. Sometimes, pure Python gets the job done—just not elegantly. 😅

## Meta-analysis with GWAMA
### General framework

First, create the working directory for the meta-analysis 
```
mkdir LPD_meta-analysis/
```
For running the meta-analysis using GWAMA, the sumstats from each cohort intended to meta-analyze must have a given format that match the headers expected by GWAMA. Once the format is achieved, run the python scripts Meta-analysis_SAIGE.py and Meta-analysis_admixkit.py. 

The process to manipulate the output sumstats from GWAMA in order to plot them is outlined as well. 

When running the scripts: Meta-analysis_SAIGE.py and Meta-analysis_admixkit.py, please take into account the following flags and options in order to run your desired meta-analysis:

    Options 

    -l, --list, help='List of hybrid files. Format: <NAME> <HYBRID> <PGEN prefix file> (mandatory) 
    -f, --folder, specify the name of the output folder you want (mandatory) 
    -n, --name, to specify the name of the output file (mandatory)
    -G, --gwama, Path to gwama (mandatory)
    -P, --plink2, Path to PLINK2 (mandatory)
    -o, --odds, Use OR instead BETA (by default, BETA) (optional)
    -F, --firth, Signalize that we used firth on PLINK2 (optional)
    -s, --sex, Run gender-differentiated and gender-heterogeneity analysis (optional)
    -r, --random, Random effect correction from GWAMA (optional, by default a fixed-effect meta-analysis is done)
    -g, --genomic, Use genomic control for adjusting studies result files from GWAMA (optional)

### GWAMA using SAIGE sumstats

Genotype-Phenotype association softwares like SAIGE and the Admix-kit suite output summary statistics for each chromosome, you must merge them in a single file (here the file is called SAIGE_ALL_1-10.txt). 
GWAMA expects the sample size in an additional column inside the sumstats, this could be fixed using a a simple command-line argument in bash like this:



```
awk 'BEGIN {OFS="\t"} NR==1 {print $0, "N"; next} {print $0, 1485}' /home/output/Association/SAIGE/Output/SAIGE_P1_ALL_1-10.txt >  /home/output/Association/SAIGE/Output/SAIGE_P1_ALL_1-10_N.txt

awk 'BEGIN {OFS="\t"} NR==1 {print $0, "N"; next} {print $0, 4303}' /home/output/Association/SAIGE/Output/SAIGE_P2_ALL_1-10.txt >  /home/output/Association/SAIGE/Output/SAIGE_P2_ALL_1-10_N.txt
```
    
    
This will add the number of samples of that study in the last column. You must do this for the total number of sumstats files intended to meta-analyze (this same process is aplicable to ATT phase 1 and 2, TRACTOR_NAT phase 1 and 2, TRACTOR_AFR phase 1 and 2, TRACTOR_EUR phase 1 and phase 2 sumstats as well).

Example of the expected header to run Meta-analysis_SAIGE.py:
    
    ID	Chrom	Pos	Allele1	Allele2	P	BETA	SE	OR	U95	L95	N
    rs200188737	1	730869	C	T	4.529689E-01	0.351674	0.468602	1.4214450570890105	3.56132946366326970.5673460068601568	1485
    rs61769351	1	758443	G	C	5.893046E-01	-0.0909793	0.168529	0.9130366105614922	1.27040727749360770.6561957468240437	1485
    rs745495619	1	762107	A	C	7.755619E-01	-0.0739648	0.259428	0.9287043836602479	1.544212512869177	0.5585318245007834	1485
    rs1221434219	1	763097	C	T	5.617923E-01	0.283	0.487779	1.3271051618171572	3.4523212814705655	0.5101518563682264	1485
    rs535793062	1	767578	T	C	3.673955E-01	0.432908	0.480281	1.541734374613348	3.95214614811429450.6014314230253184	1485
    rs142559957	1	769257	G	A	4.220278E-01	-0.126379	0.157401	0.8812807780665058	1.19976650656652820.6473391326885151	1485
    rs992276675	1	769283	G	A	5.726310E-01	0.390115	0.691471	1.477150656440849	5.728181754280549	0.3809191389908478	1485


After this, within the working directory, create a text file (input_list.txt) with the path for the input sumstats and the plink files containing the genetic information (in order to annotate the allele frequencies of each SNP to run GWAMA and downstream tools, this is made by the python script)(this same process is aplicable to admixkit meta-analysis as well):
    
```
> input_list_SAIGE_Meta_LPDs.txt
    
P1 SAIGE_P1_ALL_1-10_N.txt /home/output/Association/SAIGE/PLINK/LPD_Phase1.bed/bim.fam
P2 SAIGE_P2_ALL_1-10_N.txt /home/output/Association/SAIGE/PLINK/LPD_Phase2.bed/bim.fam
```
   


Now you are ready to run the Meta-analysis_SAIGE.py script with the adequate flags depending on the intended meta-analysis. This is the example for running a random-effects meta-analysis with SAIGE results, without genomic control and indicating that the Firth correction was implemented (SAIGE includes this correction in the logistic regression itself):
    
```
python Meta-analysis_SAIGE_results.py -l input_list_SAIGE_Meta_LPDs.txt -f GWAMA_SAIGE -n GWAMA_SAIGE_results -o -F -r -G /home/programs/GWAMA -P /home/programs/plink2
```

This script will reformat the sumstats in the way GWAMA expects them, it generates two files called inputP1.in and inputP2.in, for phase 1 and phase 2 in each meta-analysis. This files will prove to be useful when running the Meta-analysis with meta-regression with SAS PD GWASs and LARGE-PD cohorts (see further section)

    
<details>
    <summary>Meta-analysis_SAIGE.py</summary>
            
```python=
import numpy as np
import argparse
import os

def openListFile(fileName, sex):
    dictFiles = {}

    listFile = open(fileName)

    for fileLine in listFile:
        split = fileLine.strip().split()
        if not sex:
            if len(split) == 3:
                dictFiles[split[0]] = {}
                dictFiles[split[0]]["hybrid"] = split[1]
                dictFiles[split[0]]["pgen"] = split[2]
            else:
                print(f"Ignoring line {fileLine} because it has only ({len(split)} columns)")
        else:
            if len(split) == 4:
                dictFiles[split[0]] = {}
                dictFiles[split[0]]["hybrid"] = split[1]
                dictFiles[split[0]]["pgen"] = split[2]
                dictFiles[split[0]]["sex"] = split[3]
            else:
                print(f"Ignoring line {fileLine} because it has only ({len(split)} columns)")
    print (dictFiles)
    return dictFiles
    

def execute(command):
    print(command)
    os.system(command)

def calculateFrq(fileDict, plink2, folder):
    freqDict = {}
    for pop in fileDict:
        pgenFile = fileDict[pop]["pgen"]

        command = f"{plink2} --pfile {pgenFile} --freq --out {folder}/{pop}_freq"
        execute(command)

        # open the .afreq file once, read its header, and build a name→index map
        path = f"{folder}/{pop}_freq.afreq"
        with open(path) as file:
            header = file.readline().strip().split()
            idx    = { name: i for i, name in enumerate(header) }

            # now process each data line exactly as before…
            for line in file:
                split = line.strip().split()

                # SNP ID is still at split[1], but you *could* also do:
                #   snp = split[idx["ID"]]
                snp = split[1]

                if snp not in freqDict:
                    freqDict[snp] = {}
                if pop not in freqDict[snp]:
                    freqDict[snp][pop] = {}

                    # replace hard-coded 3 and 4 with the header names
                    freqDict[snp][pop]["ALT"]      = split[idx["ALT"]]
                    freqDict[snp][pop]["ALT_FREQ"] = split[idx["ALT_FREQS"]]

    return freqDict

def getStatus(fileDict):
    statusDict = {}
    for pop in fileDict:
        pgenFile = fileDict[pop]["pgen"]
        pvarFile = open(f"{pgenFile}.pvar")

        header = True
        for line in pvarFile:
            if header:
                if "#CHROM" in line:
                    header = False
            else:
                split = line.strip().split()
                if split[2] not in statusDict:
                    statusDict[split[2]] = {}
                else:
                    print (f"duplicate {split[2]}")
                if pop not in statusDict[split[2]]:
                    if "TYPED" in line:
                        statusDict[split[2]][pop] = 0
                    else:
                        statusDict[split[2]][pop] = 1
    return statusDict

def getPositionCol(headerLine):
    dictFields = {}
    split = headerLine.strip().split()

    interest =  ["ID", "Chrom", "Pos", "Allele1", "Allele2", "Allele2", "P", "SE", "OR", "U95", "L95", "N"]

    for i in range(0, len(split)):
        if split[i] in interest:
            print(f"Col {split[i]} -> index {i}")
            dictFields[split[i]] = i

    return dictFields

def prepareInputGWAMAOR(fileDict, freqDict, statusDict, folder, name, sex):
    fileToGWAMA = open(f"{folder}/input{name}.in", 'w')
    missingSNPs = open(f"{folder}/missing_snps.txt", 'w')

    for pop in fileDict:
        fileOut = open(f"{folder}/input{pop}.in", 'w')
        if sex:
            fileToGWAMA.write(f"{folder}/input{pop}.in\t{fileDict[pop]['sex']}\n")
        else:
            fileToGWAMA.write(f"{folder}/input{pop}.in\n")

        fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tOR\tOR_95L\tOR_95U\n")

        hybrid = fileDict[pop]["hybrid"]
        print(f"Open {hybrid} file")

        fileHybrid = open(hybrid)
        header = True

        for line in fileHybrid:
            if header:
                dictColHeader = getPositionCol(line)
                header = False
            else:
                split = line.strip().split()
                if split[dictColHeader["OR"]] != "NA":
                    ID = split[dictColHeader["ID"]]

                    # Skip missing SNPs while logging them
                    if ID not in statusDict or pop not in statusDict[ID]:
                        missingSNPs.write(f"{ID}\t{pop}\n")
                        continue

                    CHR = split[dictColHeader["Chrom"]]
                    POS = split[dictColHeader["Pos"]]
                    IMPUTED = statusDict[ID][pop]
                    N = split[dictColHeader["N"]]
                    EA = split[dictColHeader["Allele2"]]
                    NEA = split[dictColHeader["Allele1"]] if EA == split[dictColHeader["Allele2"]] else split[dictColHeader["Allele2"]]

                    if ID in freqDict and pop in freqDict[ID]:
                        if EA == freqDict[ID][pop]["ALT"]:
                            EAF = freqDict[ID][pop]["ALT_FREQ"]
                        else:
                            EAF = 1 - float(freqDict[ID][pop]["ALT_FREQ"])
                    else:
                        missingSNPs.write(f"{ID}\t{pop}\n")
                        continue

                    OR = split[dictColHeader["OR"]]
                    OR_95L = split[dictColHeader["L95"]]
                    OR_95U = split[dictColHeader["U95"]]

                    fileOut.write(f"{ID}\t{CHR}\t{POS}\t{IMPUTED}\t{N}\t{EA}\t{NEA}\t{EAF}\t{OR}\t{OR_95L}\t{OR_95U}\n")

        fileOut.close()

    missingSNPs.close()
    return f"{folder}/input{name}.in"


def prepareInputGWAMABeta(fileDict, freqDict, statusDict, folder, name, sex):
    fileToGWAMA = open(f"{folder}/input{name}.in", 'w')
    missingSNPs = open(f"{folder}/missing_snps.txt", 'w')

    for pop in fileDict:
        fileOut = open(f"{folder}/input{pop}.in", 'w')
        if sex:
            fileToGWAMA.write(f"{folder}/input{pop}.in\t{fileDict[pop]['sex']}\n")
        else:
            fileToGWAMA.write(f"{folder}/input{pop}.in\n")

        fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tBETA\tSE\n")

        hybrid = fileDict[pop]["hybrid"]
        print(f"Open {hybrid} file")

        fileHybrid = open(hybrid)
        header = True

        for line in fileHybrid:
            if header:
                dictColHeader = getPositionCol(line)
                header = False
            else:
                split = line.strip().split()
                if split[dictColHeader["OR"]] != "NA":
                    ID = split[dictColHeader["ID"]]

                    # Skip missing SNPs while logging them
                    if ID not in statusDict or pop not in statusDict[ID]:
                        missingSNPs.write(f"{ID}\t{pop}\n")
                        continue

                    CHR = split[dictColHeader["Chrom"]]
                    POS = split[dictColHeader["Pos"]]
                    IMPUTED = statusDict[ID][pop]
                    N = split[dictColHeader["N"]]
                    EA = split[dictColHeader["Allele2"]]
                    NEA = split[dictColHeader["Allele1"]] if EA == split[dictColHeader["Allele2"]] else split[dictColHeader["Allele2"]]

                    if ID in freqDict and pop in freqDict[ID]:
                        if EA == freqDict[ID][pop]["ALT"]:
                            EAF = freqDict[ID][pop]["ALT_FREQ"]
                        else:
                            EAF = 1 - float(freqDict[ID][pop]["ALT_FREQ"])
                    else:
                        missingSNPs.write(f"{ID}\t{pop}\n")
                        continue

                    BETA = np.log(float(split[dictColHeader["OR"]]))
                    SE = split[dictColHeader["SE"]]

                    fileOut.write(f"{ID}\t{CHR}\t{POS}\tIMPUTED\tN\tEA\tNEA\tEAF\tBETA\tSE\n")

        fileOut.close()

    missingSNPs.close()
    return f"{folder}/input{name}.in"


def runGWAMABeta(fileGWAMA, gwama, folder, name, randomEffectCorrection, genomicControl, sex):
    line = f"{gwama} -i {fileGWAMA} -o {folder}/{name} -qt"
    if randomEffectCorrection:
        line = f"{line} -r"
    if genomicControl:
        line = f"{line} -gc"
    if sex:
        line = f"{line} --sex"

    execute(line)

def runGWAMAOR(fileGWAMA, gwama, folder, name, randomEffectCorrection, genomicControl, sex):
    line = f"{gwama} -i {fileGWAMA} -o {folder}/{name}"
    if randomEffectCorrection:
        line = f"{line} -r"
    if genomicControl:
        line = f"{line} -gc"
    if sex:
        line = f"{line} --sex"


    execute(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GWAMA automatic')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-l', '--list', help='List of hybrid files. Format: <NAME> <HYBRID> <PGEN prefix file>',
                          required=True)
    required.add_argument('-f', '--folder', help='Output folder name', required=True)
    required.add_argument('-n', '--name', help='Name to output file', required=True)

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-o', '--odds', help='Use OR instead BETA (Default: False)', default=False, action="store_true")
    optional.add_argument('-F', '--firth', help='Signalize that we used firth on PLINK2', default=False,
                          action="store_true")
    optional.add_argument('-s', '--sex', help='Run gender-differentiated and gender-heterogeneity analysis (Default: False)', default=False,
                          action="store_true")
    optional.add_argument('-r', '--random', help='Random effect correction from GWAMA (Default: False)', default=False,
                          action="store_true")
    optional.add_argument('-g', '--genomic', help='Use genomic control for adjusting studies result files from GWAMA (Default: False)', default=False,
                          action="store_true")

    programs = parser.add_argument_group("Programs")
    programs.add_argument('-G', '--gwama', help='Path to gwama', required=True)
    programs.add_argument('-P', '--plink2', help='Path to PLINK2', required=True)
    args = parser.parse_args()

    execute(f"mkdir {args.folder}")

    fileDict = openListFile(args.list, args.sex)
    freqDict = calculateFrq(fileDict, args.plink2, args.folder)
    statusDict = getStatus(fileDict)
    if args.odds:
        fileGWAMA = prepareInputGWAMAOR(fileDict, freqDict, statusDict, args.folder, args.name, args.sex)
        runGWAMAOR(fileGWAMA, args.gwama, args.folder, args.name, args.random, args.genomic, args.sex)
    else:
        fileGWAMA = prepareInputGWAMABeta(fileDict, freqDict, statusDict, args.folder, args.name, args.sex)
        runGWAMABeta(fileGWAMA, args.gwama, args.folder, args.name, args.random, args.genomic, args.sex)

```
</details>




Particularly the script for running SAIGE meta-analysis has some changes when compared to the scripts needed to run ATT and TRACTOR meta-analyses. These specific changes account for the issue of missing rsIDs for some specific variants in the sumstats of SAIGE results. Given the recent changes of ID labeling of variants by the imputation server TopMed, a bunch of variants do not have an rsID yet, and still have the old identification label. So, when applying the formatting for GWAMA, they appear to be missing. 

In order to handle this issue, an aditional step is added inside the script: Meta-analysis_SAIGE.py. This adjustment, helps to identify the variants with the old identification label and saves them in the missing_snps.txt output file. 

Furthermore, it was necessary to inspect if some of those variants were statistically significant (because they are excluded from the analysis by GWAMA). Using the script: annotate_missing_variants_SAIGE.py, the p-value for each of those variants comming from the original sumstats is appended, so one can visually inspect if some significant variants were ineed excluded from the meta-analysis. After running the script, the variants are exported and sorted by p-value in a new text file (called annotated_missing_IDs.txt) like this:

    ID	Population	P
    chr6:2860729:C:T	P2	0.0005699132
    chr19:11654035:T:C	P2	0.0007023961
    chr19:11605131:T:C	P2	0.0007710347
    chr19:11747510:C:T	P2	0.001188544
    chr12:132435728:A:G	P2	0.00146866
    chr7:6619559:C:T	P1	0.001721198
    chr15:98118198:A:T	P2	0.001874586
    chr7:6616664:G:A	P1	0.002014589
    chr12:132459846:G:A	P2	0.002282474

In this way, is evident that none of the variants in question have a genome-wide significant p-value, so we can keep working with the output SAIGE GWAMA sumstats for plotting.
    
### GWAMA using ATT and TRACTOR sumstats
The workflow for running GWAMA with these results is the same as outlined with the SAIGE meta-analysis example. However, the step for checking the significance of variants with old identification IDs from TopMed can be avoided because within the ATT and TRACTOR workflow, this problem is not present inside the plink files generated to include local ancestry. 

For either ATT and TRACTOR meta-analysis you could use this script: Meta-analysis_admixkit.py. 
Just remember to append the adequate sample size for the UNO dataset (1419 for phase 1 and 4086 for phase 2). And, do not forget to update the names and paths to be included in the input list. Finally, use the desired flags when running the python script. 

    
<details>
    <summary>Meta-analysis_admixkit.py</summary>
            
```python=
import numpy as np
import argparse
import os

def openListFile(fileName, sex):
    dictFiles = {}

    listFile = open(fileName)

    for fileLine in listFile:
        split = fileLine.strip().split()
        if not sex:
            if len(split) == 3:
                dictFiles[split[0]] = {}
                dictFiles[split[0]]["hybrid"] = split[1]
                dictFiles[split[0]]["pgen"] = split[2]
            else:
                print(f"Ignoring line {fileLine} because it has only ({len(split)} columns)")
        else:
            if len(split) == 4:
                dictFiles[split[0]] = {}
                dictFiles[split[0]]["hybrid"] = split[1]
                dictFiles[split[0]]["pgen"] = split[2]
                dictFiles[split[0]]["sex"] = split[3]
            else:
                print(f"Ignoring line {fileLine} because it has only ({len(split)} columns)")

    return dictFiles

def execute(command):
    print(command)
    os.system(command)

def calculateFrq(fileDict, plink2, folder):
    freqDict = {}
    for pop in fileDict:
        pgenFile = fileDict[pop]["pgen"]

        command = f"{plink2} --pfile {pgenFile} --freq --out {folder}/{pop}_freq"
        execute(command)

        file = open(f"{folder}/{pop}_freq.afreq")
        for line in file:
            split = line.strip().split()
            if split[1] not in freqDict:
                freqDict[split[1]] = {}
            if pop not in freqDict[split[1]]:
                freqDict[split[1]][pop] = {}
                freqDict[split[1]][pop]["ALT"] = split[3]
                freqDict[split[1]][pop]["ALT_FREQ"] = split[4]
    return freqDict

def getStatus(fileDict):
    statusDict = {}
    for pop in fileDict:
        pgenFile = fileDict[pop]["pgen"]
        pvarFile = open(f"{pgenFile}.pvar")

        header = True
        for line in pvarFile:
            if header:
                if "#CHROM" in line:
                    header = False
            else:
                split = line.strip().split()
                if split[2] not in statusDict:
                    statusDict[split[2]] = {}
                if pop not in statusDict[split[2]]:
                    if "TYPED" in line:
                        statusDict[split[2]][pop] = 0
                    else:
                        statusDict[split[2]][pop] = 1
    return statusDict

def getPositionCol(headerLine):
    dictFields = {}
    split = headerLine.strip().split()

    interest =  ["ID", "Chrom", "Pos", "Allele1", "Allele2", "Allele2", "P", "SE", "OR", "U95", "L95", "N"]

    for i in range(0, len(split)):
        if split[i] in interest:
            print(f"Col {split[i]} -> index {i}")
            dictFields[split[i]] = i

    return dictFields

def prepareInputGWAMAOR(fileDict, freqDict, statusDict, folder, name, sex):
    fileToGWAMA = open(f"{folder}/input{name}.in", 'w')

    for pop in fileDict:
        fileOut = open(f"{folder}/input{pop}.in", 'w')
        if sex:
            fileToGWAMA.write(f"{folder}/input{pop}.in\t{fileDict[pop]['sex']}\n")
        else:
            fileToGWAMA.write(f"{folder}/input{pop}.in\n")
        #fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tBETA\tSE\n")
        fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tOR\tOR_95L\tOR_95U\n")
        #fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tOR\tOR_95L\tOR_95U\n")

        hybrid = fileDict[pop]["hybrid"]

        print(f"Open {hybrid} file")

        fileHybrid = open(hybrid)
        header = True

        for line in fileHybrid:
            if header:
                dictColHeader = getPositionCol(line)
                header = False
            else:
                split = line.strip().split()
                if split[dictColHeader["OR"]] != "NA":
                    ID = split[dictColHeader["ID"]]
                    CHR = split[dictColHeader["Chrom"]]
                    POS = split[dictColHeader["Pos"]]
                    IMPUTED = statusDict[ID][pop]
                    N = split[dictColHeader["N"]]
                    EA = split[dictColHeader["Allele2"]]
                    if EA == split[dictColHeader["Allele2"]]:
                        NEA = split[dictColHeader["Allele1"]]
                    else:
                        NEA = split[dictColHeader["Allele2"]]

                    if EA == freqDict[ID][pop]["ALT"]:
                        EAF = freqDict[ID][pop]["ALT_FREQ"]
                    else:
                        EAF = 1-float(freqDict[ID][pop]["ALT_FREQ"])

                    OR = split[dictColHeader["OR"]]
                    OR_95L = split[dictColHeader["L95"]]
                    OR_95U = split[dictColHeader["U95"]]

                    fileOut.write(f"{ID}\t{CHR}\t{POS}\t{IMPUTED}\t{N}\t{EA}\t{NEA}\t{EAF}\t{OR}\t{OR_95L}\t{OR_95U}\n")
        fileOut.close()
    return f"{folder}/input{name}.in"


def prepareInputGWAMABeta(fileDict, freqDict, statusDict, folder, name, sex):
    fileToGWAMA = open(f"{folder}/input{name}.in", 'w')

    for pop in fileDict:
        fileOut = open(f"{folder}/input{pop}.in", 'w')
        if sex:
            fileToGWAMA.write(f"{folder}/input{pop}.in\t{fileDict[pop]['sex']}\n")
        else:
            fileToGWAMA.write(f"{folder}/input{pop}.in\n")
        fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tBETA\tSE\n")
        #fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tOR\tOR_95L\tOR_95U\n")
        #fileOut.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tOR\tOR_95L\tOR_95U\n")

        hybrid = fileDict[pop]["hybrid"]

        print(f"Open {hybrid} file")

        fileHybrid = open(hybrid)
        header = True

        for line in fileHybrid:
            if header:
                dictColHeader = getPositionCol(line)
                header = False
            else:
                split = line.strip().split()
                if split[dictColHeader["OR"]] != "NA":
                    ID = split[dictColHeader["ID"]]
                    CHR = split[dictColHeader["Chrom"]]
                    POS = split[dictColHeader["Pos"]]
                    IMPUTED = statusDict[ID][pop]
                    N = split[dictColHeader["N"]]
                    EA = split[dictColHeader["Allele2"]]
                    if EA == split[dictColHeader["Allele2"]]:
                        NEA = split[dictColHeader["Allele1"]]
                    else:
                        NEA = split[dictColHeader["Allele2"]]

                    if EA == freqDict[ID][pop]["ALT"]:
                        EAF = freqDict[ID][pop]["ALT_FREQ"]
                    else:
                        EAF = 1-float(freqDict[ID][pop]["ALT_FREQ"])

                    BETA = np.log(float(split[dictColHeader["OR"]]))
                    SE = split[dictColHeader["SE"]]
                    fileOut.write(f"{ID}\t{CHR}\t{POS}\t{IMPUTED}\t{N}\t{EA}\t{NEA}\t{EAF}\t{BETA}\t{SE}\n")
        fileOut.close()
    return f"{folder}/input{name}.in"

def runGWAMABeta(fileGWAMA, gwama, folder, name, randomEffectCorrection, genomicControl, sex):
    line = f"{gwama} -i {fileGWAMA} -o {folder}/{name} -qt"
    if randomEffectCorrection:
        line = f"{line} -r"
    if genomicControl:
        line = f"{line} -gc"
    if sex:
        line = f"{line} --sex"

    execute(line)

def runGWAMAOR(fileGWAMA, gwama, folder, name, randomEffectCorrection, genomicControl, sex):
    line = f"{gwama} -i {fileGWAMA} -o {folder}/{name}"
    if randomEffectCorrection:
        line = f"{line} -r"
    if genomicControl:
        line = f"{line} -gc"
    if sex:
        line = f"{line} --sex"


    execute(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GWAMA automatic')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-l', '--list', help='List of hybrid files. Format: <NAME> <HYBRID> <PGEN prefix file>',
                          required=True)
    required.add_argument('-f', '--folder', help='Output folder name', required=True)
    required.add_argument('-n', '--name', help='Name to output file', required=True)

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-o', '--odds', help='Use OR instead BETA (Default: False)', default=False, action="store_true")
    optional.add_argument('-F', '--firth', help='Signalize that we used firth on PLINK2', default=False,
                          action="store_true")
    optional.add_argument('-s', '--sex', help='Run gender-differentiated and gender-heterogeneity analysis (Default: False)', default=False,
                          action="store_true")
    optional.add_argument('-r', '--random', help='Random effect correction from GWAMA (Default: False)', default=False,
                          action="store_true")
    optional.add_argument('-g', '--genomic', help='Use genomic control for adjusting studies result files from GWAMA (Default: False)', default=False,
                          action="store_true")

    programs = parser.add_argument_group("Programs")
    programs.add_argument('-G', '--gwama', help='Path to gwama', required=True)
    programs.add_argument('-P', '--plink2', help='Path to PLINK2', required=True)
    args = parser.parse_args()

    execute(f"mkdir {args.folder}")

    fileDict = openListFile(args.list, args.sex)
    freqDict = calculateFrq(fileDict, args.plink2, args.folder)
    statusDict = getStatus(fileDict)
    if args.odds:
        fileGWAMA = prepareInputGWAMAOR(fileDict, freqDict, statusDict, args.folder, args.name, args.sex)
        runGWAMAOR(fileGWAMA, args.gwama, args.folder, args.name, args.random, args.genomic, args.sex)
    else:
        fileGWAMA = prepareInputGWAMABeta(fileDict, freqDict, statusDict, args.folder, args.name, args.sex)
        runGWAMABeta(fileGWAMA, args.gwama, args.folder, args.name, args.random, args.genomic, args.sex)

```
</details>


## Meta-analysis with Meta-regression: SAS PD GWASs and LPD   

In order to assess the independency of LPD lead signal inside ITPKB gene, perform a random-effects meta-analysis coupled with MR-MEGA. Using all the availabe SAS PD GWASs and both LPD cohorts. 
    
To begin with, create the working directory to run the formal analysis

    mkdir /home/Meta_LPD_SAS

And create a sub-folder to parse the summary statistics for running GWAMA and MR-MEGA:

    mkdir reformat_sumstats

### Reformat summary statistics
    
Afterwards, inside the reformat_sumstats directory, download the summary statistics from the SAS PD GWASs by Kishore et al., 2025 and Andrews et al., 2024 from public repositories. And for LPD, use the reformated sumstats from SAIGE that were generated during the meta-analysis of LPD phase 1 and phase 2 (inputP1.in and inputP2.in)

All summary statistics are in hg38
    
    cp /home/LPD_Meta-analysis/GWAMA_SAIGE/inputP*.in ./

    
Use the scripts below to parse each summary statistics accordignly.

Please considert that Andrews et al., 2024 summary statistics lack a specific ID fo each variant. However, for the analysis a consistent marker ID is needed. Hence, a homogeneous ID for each marker across all the summary statistics needs to be constructed.
    
To parse Kishore et al., sumstats (append OR, 95%CI and change header names)

<details>
    <summary>Parse_Kishore_SAS_PD_GWAS.py</summary>
            
```python=
import pandas as pd
import numpy as np
import gzip

 
INPUT_FILE = "Kishore_SAS_PD_GWAS.txt"  # or .gz 
OUTPUT_FILE = "Kishore_SAS_PD_GWAS_modified.txt"
EXCLUDED_ROWS_FILE = "excluded_rows_kishore_gwas.txt" # must be empty, just to double check
EXCLUSION_SUMMARY_FILE = "exclusion_summary_kishore_gwas.txt"  # must be empty, just to double check


if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
else:
    df = pd.read_csv(INPUT_FILE, sep="\t")

# Add sample size column
df["N"] = 11170


# Compute new OR and 95% CI 
df["OR"] = np.exp(df["BETA"])
df["OR_95U"] = np.exp(df["BETA"] + 1.96 * df["SE"])
df["OR_95L"] = np.exp(df["BETA"] - 1.96 * df["SE"])

#  Rename columns 
df = df.rename(columns={
    "SNP": "MARKERNAME",
    "A1": "EA",
    "A2": "NEA"
})

#  Reorder columns 
final_cols = ["MARKERNAME", "CHR", "POS", "N", "EA", "NEA", "EAF", "OR", "OR_95L", "OR_95U"]
df = df[final_cols]

# Handle missing values 
is_missing = df.isnull()
excluded_rows = df[is_missing.any(axis=1)].copy()
included_rows = df.dropna()


with open(EXCLUDED_ROWS_FILE, "w") as excl_file:
    excl_file.write("ID\tMissing_Columns\n")
    for idx, row in excluded_rows.iterrows():
        missing_cols = is_missing.loc[idx]
        missing_fields = ",".join(missing_cols[missing_cols].index)
        excl_file.write(f"{row['ID']}\t{missing_fields}\n")


missing_summary = is_missing.sum()
missing_summary = missing_summary[missing_summary > 0]
with open(EXCLUSION_SUMMARY_FILE, "w") as summary_file:
    summary_file.write("Column\tMissing_Entries\n")
    for col, count in missing_summary.items():
        summary_file.write(f"{col}\t{count}\n")


included_rows.to_csv(OUTPUT_FILE, sep="\t", index=False)



```
</details>

To construct the new variant ID specifically for LPD and Kishore sumstats

<details>
    <summary>Create_newVariantID_LPD_Kishore.py</summary>
            
```python=

import pandas as pd
from pathlib import Path

# 1) List the files you want to process
files = [
    "inputP1.in",
    "inputP2.in",
    "Kishore_SAS_PD_GWAS_modified.txt"
]

for filepath in files:
    path = Path(filepath)

    # 2) Read in the file 
    df = pd.read_csv(path, sep="\t")

    # 3) Rebuild the MARKERNAME column
    #    This overwrites the existing column in place
    df["MARKERNAME"] = (
        "chr"
        + df["CHR"].astype(str)  
        + ":"
        + df["POS"].astype(str)  
        + ":"
        + df["NEA"]             
        + ":"
        + df["EA"]         
    )

    # 4) (Optional) Re-order so MARKERNAME is first
    cols = ["MARKERNAME"] + [c for c in df.columns if c != "MARKERNAME"]
    df = df[cols]

    # 5) Write out to a new file
    out = path.with_name(path.stem + "_new_snpID" + path.suffix)
    df.to_csv(out, sep="\t", index=False)
   

```
</details>

To parse Andrews et al., 2024 sumstats (append OR, 95%CI, create variant ID and change header names):
    
<details>
    <summary>Parse_Andrews_SAS_PD_GWAS.py</summary>
            
```python=
import pandas as pd
import numpy as np
import gzip


INPUT_FILE = "Andrews_SAS_PD_GWAS.txt"  # or .gz #the zipped file should not output any excluded rows because the one that has missing values is the plain .txt file 
OUTPUT_FILE = "Andrews_SAS_PD_GWAS_modified.txt"
EXCLUDED_ROWS_FILE = "excluded_rows_other_indian_gwas_with_P.txt" #must be empty, just to double check
EXCLUSION_SUMMARY_FILE = "exclusion_summary_other_indian_gwas_with_P.txt" #must be empty, jut to double check

if INPUT_FILE.endswith(".gz"):
    with gzip.open(INPUT_FILE, 'rt') as f:
        # any run of spaces or tabs splits into columns
        df = pd.read_csv(f, delim_whitespace=True)
else:
    df = pd.read_csv(INPUT_FILE, delim_whitespace=True)

#  Add sample size column 
df["N"] = 1878

#  Compute new OR and 95% CI 
df["OR"] = np.exp(df["BETA"])
df["OR_95U"] = np.exp(df["BETA"] + 1.96 * df["SE"])
df["OR_95L"] = np.exp(df["BETA"] - 1.96 * df["SE"])

#  Rename columns 
df = df.rename(columns={
    'BP': 'POS',
    'A1': 'EA',
    'A2': 'NEA',
    'MAF': 'EAF'
})

df['MARKERNAME'] = (
    'chr' 
    + df['CHR'].astype(str) 
    + ':' 
    + df['POS'].astype(str) 
    + ':' 
    + df['NEA'] 
    + ':' 
    + df['EA']
)

#  Reorder columns 
final_cols = ["MARKERNAME", "CHR", "POS", "P", "N", "EA", "NEA", "EAF", "OR", "OR_95L", "OR_95U"]
df = df[final_cols]

#  Handle missing values 
is_missing = df.isnull()
excluded_rows = df[is_missing.any(axis=1)].copy()
included_rows = df.dropna()


with open(EXCLUDED_ROWS_FILE, "w") as excl_file:
    excl_file.write("ID\tMissing_Columns\n")
    for idx, row in excluded_rows.iterrows():
        missing_cols = is_missing.loc[idx]
        missing_fields = ",".join(missing_cols[missing_cols].index)
        excl_file.write(f"{row['ID']}\t{missing_fields}\n")


missing_summary = is_missing.sum()
missing_summary = missing_summary[missing_summary > 0]
with open(EXCLUSION_SUMMARY_FILE, "w") as summary_file:
    summary_file.write("Column\tMissing_Entries\n")
    for col, count in missing_summary.items():
        summary_file.write(f"{col}\t{count}\n")


included_rows.to_csv(OUTPUT_FILE, sep="\t", index=False)

```
</details>

### Run GWAMA with LPDs and SAS PD GWASs

Create the input file for GWAMA with the path to the summary statistics 
    
    > gwama.in
    
    /home/Meta_LPDs_SAS/reformat_sumstats_inputP1.in
    /home/Meta_LPDs_SAS/reformat_sumstats_inputP2.in
    /home/Meta_LPDs_SAS/Kishore_SAS_PD_GWAS_modified.txt
    /home/Meta_LPDs_SAS/Andrews_SAS_PD_GWAS_modified.txt
    
And run a random-effects meta-analysis with double genomic control (we decided this based on significant inflation when genomic control was not applied)
    
    /home/programs/GWAMA -i gwama.in -o Meta-random_LPDs_SASs -r -gc -gco

### Run MR-MEGA with LPDs and SAS PD GWASs

Create the input file for MR-MEGA with the path to the summary statistics 
    
    > MR-MEGA.in
    
    /home/Meta_LPDs_SAS/reformat_sumstats_inputP1.in_new_snpID.txt
    /home/Meta_LPDs_SAS/reformat_sumstats_inputP2.in_new_snpID.txt
    /home/Meta_LPDs_SAS/Kishore_SAS_PD_GWAS_modified_new_snpID.txt
    /home/Meta_LPDs_SAS/Andrews_SAS_PD_GWAS_modified.txt

Then run MR-MEGA with double genomic control as well and specifying the names of the chromosome and position columns (MR-MEGA assumes different strings here than GWAMA)
    
    /home/programs/MR-MEGA -i MR-MEGA_input.in --pc 1 -o MR-MEGA_LPDs_SASs --name_chr CHR --name_pos POS --debug --gc --gco
    
## Fine Mapping 

Leveraging the bayes factor that MR-MEGA outputs, we performed bayesian fine mapping building the 95% and 99% credible sets for the genomic risk region of ITPKB.

First, create the working directory
    
    mkdir /home/Fine_Mapping

Then after performing MR-MEGA and uploading it to FUMA, download the results from FUMA portal

To perform bayesian fine mapping through the construction of 95% and 99% credible sets based on the posterior probability of each snp being the putative causal variant. 
    
    python fine_map.py --fuma-dir ./ --mrmega-file /home/Meta_LPD_SAS/run_MR-MEGA_LPDs_SAS/MR-MEGA_LPDs_SASs --out-dir ./fine_mapping_results/

<details>
    <summary>fine_map.py</summary>
            
```python=

import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2


def load_fuma_and_summary(fuma_dir: str, mrmega_file: str):
    """
    Load FUMA risk-locus definitions and MR-MEGA summary statistics.
    """
    loci_path = os.path.join(fuma_dir, "GenomicRiskLoci.txt")
    if not os.path.exists(loci_path):
        raise FileNotFoundError(f"Cannot find FUMA loci file at {loci_path}")
    loci = pd.read_csv(loci_path, delim_whitespace=True)

    if not os.path.exists(mrmega_file):
        raise FileNotFoundError(f"Cannot find MR-MEGA file at {mrmega_file}")
    sumstat = pd.read_csv(
        mrmega_file,
        delim_whitespace=True,
        compression='infer',
        engine='c'
    )
    return loci, sumstat


def finemap_locus(loci_df, sumstat_df):
    """
    Compute 95% and 99% credible sets for each locus based on Bayes Factors.
    Returns report, cs95_df, cs99_df.
    """
    # Convert lnBF to BF
    sumstat_df['BF'] = np.exp(sumstat_df['lnBF'])
    report_rows = []
    cs95_list = []
    cs99_list = []
    grouped = {chrom: df for chrom, df in sumstat_df.groupby('Chrom')}

    for idx, row in loci_df.iterrows():
        chrom = row['chr']
        start = row['start']
        end = row['end']
        df_chr = grouped.get(chrom)
        if df_chr is None:
            continue
        window = df_chr[(df_chr['Pos'] >= start) & (df_chr['Pos'] <= end)].copy()
        window = window.sort_values('BF', ascending=False).reset_index(drop=True)
        total_bf = window['BF'].sum()

        # 95% credible set
        cum = 0.0
        sel95 = []
        for i, bf in enumerate(window['BF']):
            cum += bf
            sel95.append(i)
            if cum / total_bf >= 0.95:
                break
        cs95 = window.loc[sel95].copy()
        cs95['PP'] = cs95['BF'] / total_bf
        cs95['Locus'] = idx + 1
        cs95_list.append(cs95)

        # 99% credible set
        cum = 0.0
        sel99 = []
        for i, bf in enumerate(window['BF']):
            cum += bf
            sel99.append(i)
            if cum / total_bf >= 0.99:
                break
        cs99 = window.loc[sel99].copy()
        cs99['PP'] = cs99['BF'] / total_bf
        cs99['Locus'] = idx + 1
        cs99_list.append(cs99)

        # Collect SNP lists per locus
        snps95 = window.loc[sel95, 'MarkerName'] if 'MarkerName' in window.columns else window.loc[sel95, 'rs_number']
        snps99 = window.loc[sel99, 'MarkerName'] if 'MarkerName' in window.columns else window.loc[sel99, 'rs_number']

        report_rows.append({
            'Locus': idx + 1,
            'NumSNPs': window.shape[0],
            'NumSNPs_in_95%': len(sel95),
            'NumSNPs_in_99%': len(sel99),
            'SNPs_in_95%': ';'.join(snps95.astype(str)),
            'SNPs_in_99%': ';'.join(snps99.astype(str))
        })

    report = pd.DataFrame(report_rows)
    cs95_df = pd.concat(cs95_list, ignore_index=True) if cs95_list else pd.DataFrame()
    cs99_df = pd.concat(cs99_list, ignore_index=True) if cs99_list else pd.DataFrame()
    return report, cs95_df, cs99_df


def annotate_with_fuma(cs_df, fuma_dir: str):
    """
    Annotate a credible-set DataFrame with FUMA snp annotations.
    """
    snp_file = os.path.join(fuma_dir, 'snps.txt')
    if not os.path.exists(snp_file):
        raise FileNotFoundError(f"Cannot find FUMA SNP file at {snp_file}")
    fuma = pd.read_csv(snp_file, delim_whitespace=True)
    fuma = fuma.rename(columns={'chr': 'CHR', 'pos': 'BP'})
    annotated = cs_df.merge(
        fuma,
        left_on=['Chrom', 'Pos'],
        right_on=['CHR', 'BP'],
        how='left'
    )
    return annotated


def main():
    parser = argparse.ArgumentParser(description="Fine-map MR-MEGA results using FUMA loci.")
    parser.add_argument('--fuma-dir', required=True, help='Path to FUMA output directory')
    parser.add_argument('--mrmega-file', required=True, help='Path to MR-MEGA summary stats (.tsv.gz)')
    parser.add_argument('--out-dir', default='output', help='Directory for fine-mapping outputs')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    loci, sumstat = load_fuma_and_summary(args.fuma_dir, args.mrmega_file)
    report, cs95, cs99 = finemap_locus(loci, sumstat)

    # Save reports
    report_file = os.path.join(args.out_dir, 'finemap_report.tsv')
    cs95_file   = os.path.join(args.out_dir, 'finemap_cs_95.tsv')
    cs99_file   = os.path.join(args.out_dir, 'finemap_cs_99.tsv')
    report.to_csv(report_file, sep='\t', index=False)
    cs95.to_csv(cs95_file, sep='\t', index=False)
    cs99.to_csv(cs99_file, sep='\t', index=False)

    # Annotate credible sets
    annotated95 = annotate_with_fuma(cs95, args.fuma_dir)
    annotated99 = annotate_with_fuma(cs99, args.fuma_dir)
    annotated95_file = os.path.join(args.out_dir, 'finemap_cs_95_annot.tsv')
    annotated99_file = os.path.join(args.out_dir, 'finemap_cs_99_annot.tsv')
    annotated95.to_csv(annotated95_file, sep='\t', index=False)
    annotated99.to_csv(annotated99_file, sep='\t', index=False)

    # Identify singleton loci and their rsIDs
    single95_loci = report[report['NumSNPs_in_95%'] == 1]['Locus']
    single95_snps = cs95[cs95['Locus'].isin(single95_loci)][['Locus', 'rs_number']]
    single95_snps = single95_snps.rename(columns={'rs_number': 'rsID'})
    single95_file = os.path.join(args.out_dir, 'singleton_loci_95.tsv')
    single95_snps.to_csv(single95_file, sep='\t', index=False)

    single99_loci = report[report['NumSNPs_in_99%'] == 1]['Locus']
    single99_snps = cs99[cs99['Locus'].isin(single99_loci)][['Locus', 'rs_number']]
    single99_snps = single99_snps.rename(columns={'rs_number': 'rsID'})
    single99_file = os.path.join(args.out_dir, 'singleton_loci_99.tsv')
    single99_snps.to_csv(single99_file, sep='\t', index=False)

if __name__ == "__main__":
    main()

```
</details>
    

The meta-analysis, MR-MEGA and fine mapping was performed by Felipe Duarte Zambrano (felipe.dz.101@gmail.com).

# Folder structure
For visual clarity, we just listed files that are mentioned in any command line. A lot of steps produces multiples intermediary files that was ommited.
```bash
home/
├── input/
|    ├── idat/
|    ├── calling/
|    |    ├── recluster_09092022.egt 
|    |    ├── sampleSheet.csv
|    |    ├── NeuroBooster_20042459_A4.bpm 
|    |    └── exclude.txt
|    ├── Demographic/
|    |    ├── Country.txt
|    |    ├── Covar.txt
|    |    ├── CovarToSAIGE.txt
|    |    └── CovarToAK.txt
|    ├── Imputed/
|    |    ├── chr<chrom>.dose.vcf.gz.csi
|    |    └── chr<chrom>.dose.vcf.gz
|    ├── Phased/
|    |    ├── chr<chrom>.reheaderphased.vcf.gz.csi
|    |    └── chr<chrom>.reheaderphased.vcf.gz 
|    └── LA/
|         ├── configBest.yaml
|         └── G_chr<chrom>.b38.gmap
├── output/
|    ├── genotyping/
|    |    ├── NeuroBooster_20042459_A4.ped/map
|    |    ├── <Sentrix_Barcode>_<Sentrix_Position>.ped
|    |    ├── NB.bed/bim/fam
|    |    └── NB_NoLift.bed/bim/fam
|    ├── QCed/
|    |    ├── main.py
|    |    ├── LARGE_NAToRA_Removal.txt 
|    |    ├── Test.vcf.gz
|    |    ├── FixRef.vcf.gz 
|    |    ├── LARGE_chr<chrom>_ToTOPMED.vcf.gz
|    |    └── LARGE_HWECase.pvar/pgen/psam
|    ├── Admixture/
|    |    ├── prepareAdmixture.py 
|    |    ├── LD_Bfile1And2.bed/bim/fam/pop
|    |    ├── LD_Bfile1And2.P/Q
|    |    └── Common_Bfile2.bed/bim/fam
|    ├── PCA/
|    |    ├── covarProjectedPCA.py
|    |    ├── selectModel.R
|    |    ├── Projected/
|    |    |   ├── outlier.py
|    |    |   ├── toRemovePC1-10.txt
|    |    |   └── AutosomalPCA_TargetPC.proj.eigenvec
|    |    ├── NonProjected/
|    |    |   └── AutosomalPCA_PCs.eigenvec
|    |    └── NonProjectedNonOutiler/
|    |        └── AutosomalPCA_PCs.eigenvec
|    ├── LA/
|    |    ├── Correspondence.txt
|    |    ├── Typed/
|    |    |   └── CHR<chrom>_OnlyTyped.vcf.gz
|    |    ├── Panels/
|    |    |   └── VCF
|    |    |        ├── Ref
|    |    |        |   └── REFERENCE_chr<chrom>.vcf.gz
|    |    |        └── Target
|    |    |            └── TARGET_chr<chrom>.vcf.gz
|    |    └── Inference/
|    |        └── LPD_chrom<chrom>/
|    |            └── query_results.msp
|    └── Association/
|         ├── SAIGE/
|         |   ├── GRM_relatednessCutoff_0.088_500000_randomMarkersUsed.sparseGRM.mtx
|         |   ├── GRM_relatednessCutoff_0.088_500000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
|         |   ├── Model_LPD_Phase2_PC1-10.rda
|         |   ├── Model_LPD_Phase2_PC1-10.varianceRatio.txt
|         |   ├── PLINK
|         |   |    └── LPD_Phase2.bed/bim.fam
|         |   └── Output
|         |        └──Result_LPD_Phase2_Chrom<chrom>
|         |        └──Result_LPD_Phase2_Chrom<chrom>.index
|         ├── AK/
|         |   ├── removeToAK.txt
|         |   ├── PGEN
|         |   |   ├── MSP_<chrom>.msp
|         |   |   └── LPD_<chrom>.pvar/.pgen/.psam/.lanc
|         |   └──  Results
|         |       ├── LPD_ATT_<chrom>.ATT.assoc
|         |       ├── LPD_TRACTOR_<chrom>.TRACTOR.assoc
|         |       └──  LPD_ADM_<chrom>.ADM.assoc
|         └── AM/
|             └── Phase2/
|                 └── AM_Results
|                     └── AM_Phase2_CHROM<chrom>_<anc>.csv
|
|
|
|
|── LPD_meta-analysis/
|    ├── SAIGE_P1_ALL_1-10_N.txt
|    ├── SAIGE_P2_ALL_1-10_N.txt
|    ├── ATT_P2_ALL_1-10_N.txt 
|    ├── ATT_P2_ALL_1-10_N.txt
|    ├── TRACTOR_NAT_P2_ALL_1-10_N.txt
|    ├── TRACTOR_NAT_P2_ALL_1-10_N.txt
|    ├── TRACTOR_AFR_P2_ALL_1-10_N.txt
|    ├── TRACTOR_AFR_P2_ALL_1-10_N.txt
|    ├── TRACTOR_EUR_P2_ALL_1-10_N.txt
|    ├── TRACTOR_EUR_P2_ALL_1-10_N.txt
|    ├── GWAMA_SAIGE
|    |    └── GWAMA_SAIGE_results
|    |    └── inputP1.in
|    |    └── inputP2.in
|    ├── GWAMA_ATT
|    |    └── GWAMA_ATT_results
|    |    └── inputP1.in
|    |    └── inputP2.in
|    ├── GWAMA_TRACTOR_NAT
|    |    └── GWAMA_TRACTOR_NAT_results
|    |    └── inputP1.in
|    |    └── inputP2.in
|    ├── GWAMA_TRACTOR_AFR
|    |    └── GWAMA_TRACTOR_AFR_results
|    |    └── inputP1.in
|    |    └── inputP2.in
|    ├── GWAMA_TRACTOR_EUR
|    |    └── GWAMA_TRACTOR_EUR_results
|    |    └── inputP1.in
|    |    └── inputP2.in
|
|
|── Meta_LPDs_SAS/
|    ├── run_GWAMA_LPDs_SAS
|    ├── run_MR-MEGA_LPDs_SAS
|    ├── reformat_sumstats
|    |    └── inputP1.in
|    |    └── inputP2.in
|    |    └── inputP1.in_new_snpID.txt
|    |    └── inputP2.in_new_snpID.txt
|    |    └── Kishore_SAS_PD_GWAS_modified.txt
|    |    └── Kishore_SAS_PD_GWAS_modified_new_snpID.txt
|    |    └── Andrews_SAS_PD_GWAS.txt
|    |    └── Andrews_SAS_PD_GWAS_modified.txt
|
|
|── Fine_Mapping/
|    ├── FUMA_jobID.zip
|    ├── fine_mapping_results
|    |      └── finemap_report.tsv
|    
|
├── references/
|    ├── AllRef.txt
|    ├── Parentals
|    |    ├── AFR.txt
|    |    ├── AMR.txt
|    |    ├── EAS.txt
|    |    ├── EUR.txt
|    |    └── SAS.txt
|    ├── 1KGP
|    |    ├── VCF
|    |    |    └── /home/references/1KGP/PGEN/1KGP30x_chrom<chrom>.vcf.gz
|    |    └── PGEN
|    |        ├── OneThousand_chrom<chrom>.pvar/pgen/psam
|    |        ├── ListOne.txt
|    |        └── OneThousand_All.pvar/pgen/psam
|    └── hg38
|         └── Homo_sapiens_assembly38.fasta
└── programs/
     ├── NAToRA_Public.py
     ├── plink
     ├── plink2
     ├── king
     ├── admixture
     ├── gcta
     ├── bcftools
     ├── GWAMA
     ├── MR-MEGA
     └── gnomix
          └── gnomix.py
```

# Software
| Software           | Version | Resource URL                                                                                                                                 | RRID | Notes                                                                                                      |
| ------------------ | ------- | -------------------------------------------------------------------------------------------------------------------------------------------- | ---- | ---------------------------------------------------------------------------------------------------------- |
| Python             | \> 3.8  | [http://www.python.org/](http://www.python.org/)                                                                                             |  RRID:SCR_008394    | pandas, numpy, matplotlib, statsmodel, networkX, argparse. Used for general analysis and data manipulation |
| R                  | \> 4    | [http://www.r-project.org/](http://www.r-project.org/)                                                                                       |   RRID:SCR_001905 | MASS, tidyverse, dplyr, tidyr, ggplot, GENESIS, data.table. Used for statistical analysis                  |
| IAAP               |         | [https://emea.support.illumina.com/downloads/iaap-genotyping-cli.html](https://emea.support.illumina.com/downloads/iaap-genotyping-cli.html) |      | Used for genotyping calling                                                                                |
| BCFTOOLS           |    RRID:SCR_005227     | [https://samtools.github.io/bcftools/howtos/install.html](https://samtools.github.io/bcftools/howtos/install.html)                           |      | Used for VCF manipulation                                                                                  |
| PLINK              | 1.9     | [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)                                                                   |   RRID:SCR_001757   | Used for genetic analyses                                                                                  |
| PLINK2             | 2       | [https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/)                                                               |      | Used for genetic analyses                                                                                  |
| GCTA               |         | [https://yanglab.westlake.edu.cn/software/gcta/#Overview](https://yanglab.westlake.edu.cn/software/gcta/#Overview)                           |      | Used for PCA inference                                                                                     |
| ADMIXTURE          |     | [https://dalexander.github.io/admixture/](https://dalexander.github.io/admixture/)                                                           |    RRID:SCR_001263  | Used for Individual ancestry inference                                                                     |
| GWAMA          | 2.2.2    | [https://genomics.ut.ee/en/tools](https://genomics.ut.ee/en/toolshttps://dalexander.github.io/admixture/)                                                           |    RRID:SCR_006624 | Used for Meta-analysis                                                                     |
| Gnomix             |         | [https://github.com/AI-sandbox/gnomix](https://github.com/AI-sandbox/gnomix)                                                                 |      | Used for Local Ancestry inference                                                                          |
| Admix-kit          |         | [https://github.com/KangchengHou/admix-kit](https://github.com/KangchengHou/admix-kit)                                                       |      | Used for TRACTOR/ATT analysis                                                                              |
| NAToRA             |         | [https://github.com/ldgh/NAToRA_Public](https://github.com/ldgh/NAToRA_Public)                                                               |      | Used for genetic relationship removal                                                                      |
| Mata Lab GWASQC    |         | [https://github.com/MataLabCCF/GWASQC](https://github.com/MataLabCCF/GWASQC)                                                                 |      | Used for genetic qc                                                                                        |
| Mata Lab PCA       |         | [https://github.com/MataLabCCF/ProjectedPCAAndModelSelection](https://github.com/MataLabCCF/ProjectedPCAAndModelSelection)                   |      | Automated pipeline to perform PCA                                                                          |
| Mata Lab Admixture |         | [https://github.com/MataLabCCF/AdmixtureWithReference](https://github.com/MataLabCCF/AdmixtureWithReference)                                 |      | Automated pipeline to perform supervised ADMIXTURE                                                         |
| Mata Lab AM        |         | [https://github.com/MataLabCCF/GenesisAM](https://github.com/MataLabCCF/GenesisAM)                                                           |      | Automated pipeline to perform Admixture Mapping                                                            |


## Acknowledgements
This work is supported by NIH Grant R01 1R01NS112499-01A1, MJFF Grant ID: 18298, ASAP-GP2, Parkinson’s Foundation, and  Veterans Affair Puget Sound Healthcare System (5I01ABX005978-2).

