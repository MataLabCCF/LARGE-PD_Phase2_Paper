# LARGE-PD_Phase2_Paper
This repository contains a detailed description of all analyses performed as part of the LARGE-PD Phase 2 manuscript. The Latin American Research Consortium on the Genetics of Parkinson’s Disease (LARGE-PD) is an international collaborative effort aimed at understanding the genetic architecture of Parkinson's disease in Latin American populations, with a particular focus on admixed and underrepresented groups.

The analyses documented here include:

🧬 Genotyping Quality Control (QC)
🌍 Global and Local Ancestry Inference
📊 Genome-Wide Association Studies (GWAS)
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
