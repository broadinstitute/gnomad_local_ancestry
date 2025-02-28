# gnomAD Local Ancestry Inference (LAI) Pipeline

This repository provides a streamlined pipeline for performing **Local Ancestry Inference (LAI)** using gnomAD samples. Due to the large scale of gnomAD, we implemented this pipeline using the **Hail Batch Python module**. The pipeline leverages **Eagle for phasing**, **RFMix v2 for ancestry painting**, **Tractor for extracting ancestry-specific allele frequencies**, and **generate\_output\_vcf.py** for generating a joint VCF with ancestry-specific allele count (AC), allele number (AN), and allele frequency (AF) annotations. Below, we outline the steps to run this pipeline on your dataset using Hail/Python scripts.

---

## **Getting Started**

To run the LAI pipeline, set up your environment and install the required dependencies. Follow the step-by-step instructions below.

---

## **1. Pipeline Overview**

This pipeline processes genomic data by:

- **Phasing** haplotypes using **Eagle**
- **Inferring local ancestry** using **RFMix v2**
- **Extracting ancestry-specific allele frequencies** from phased and painted data using **Tractor**
- **Generating a joint VCF with ancestry-specific calls** using **generate\_output\_vcf.py**

---

## **2. Installation & Setup**

### **Step 1: Clone the Repository**

```bash
git clone https://github.com/broadinstitute/gnomad_local_ancestry.git
cd gnomad_local_ancestry/batch
```

### **Step 2: Install Dependencies**

The pipeline requires **Python 3**, **Hail**, and several additional tools. Install the necessary dependencies using:

```bash
pip install hail
pip install numpy pandas
```

Make sure you have Eagle and RFMix v2 installed. You can find installation instructions and a toy dataset in the **[Tractor Tutorial](https://github.com/Atkinson-Lab/Tractor-tutorial/tree/main)**.

---

## **3. Running the LAI Pipeline**

### **Step 1: Phasing**

To phase your genotype data using **Eagle**, run:

```bash
eagle --vcf input_data.vcf.gz --out phased_data.vcf.gz
```

Refer to the **[Tractor wiki](https://github.com/Atkinson-Lab/Tractor-tutorial/tree/main)** for a detailed guide on phasing.

### **Step 2: Local Ancestry Inference using RFMix**

After phasing, run **RFMix v2** to infer local ancestry:

```bash
rfmix \
  -f phased_data.vcf.gz \
  -r reference_panel.vcf.gz \
  -m samplemap.txt \
  -g geneticmap.txt \
  -o painted_lai \
  --chromosome=22
```

See the **[Tractor wiki](https://github.com/Atkinson-Lab/Tractor-tutorial/tree/main)** for additional instructions on reference panels and sample maps.

### **Step 3: Extracting Ancestry-Specific Allele Frequencies**

Once local ancestry inference is complete, extract **ancestry-specific allele frequencies** using **Tractor**:

```bash
python3 extract_tracts.py \
  --vcf phased_data.vcf.gz \
  --msp painted_lai.msp.tsv \
  --num-ancs 2
```

### **Step 4: Generating a Joint VCF with Ancestry-Specific Annotations**

The `generate_lai_vcf` function calls **generate\_output\_vcf.py**, a standalone Python/Hail script. **Batch processing cannot be used** for this function, but `generate_lai_vcf` remains useful for extracting LAI-informed VCF outputs. This script outputs an annotated VCF containing ancestry-specific allele frequency data:

```bash
python3 generate_output_vcf.py \
  --msp-file painted_lai.msp.tsv \
  --tractor-output tractor_output_path \
  --output-path output_lai \
  --is-zipped \
  --mt-path-for-adj pipeline_input.mt \
  --add-gnomad-af
```

---

## **4. Additional Resources**

For detailed explanations of phasing, local ancestry painting, and extracting tracts, refer to:

- [Tractor Tutorial](https://github.com/Atkinson-Lab/Tractor-tutorial/tree/main)

For our Hail Batch Python pipeline, refer to:

- [gnomAD Local Ancestry Pipeline](https://github.com/broadinstitute/gnomad_local_ancestry/blob/main/batch/lai_batch_pipeline.py)

---

## **5. Citation**

If you use this pipeline in your research, please cite[:](https://www.biorxiv.org/content/10.1101/2024.10.30.620961v1)

> [Kore, P., Wilson, M. et al., Improved Allele Frequencies in gnomAD through Local Ancestry Inference](https://www.biorxiv.org/content/10.1101/2024.10.30.620961v1).

Please direct questions to [pragati.kore@bcm.edu](mailto\:pragati.kore@bcm.edu) or [mwilson@broadinstitute.org](mailto\:mwilson@broadinstitute.org).
