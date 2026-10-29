#!/usr/bin/env python3
import pandas as pd
from cyvcf2 import VCF
from io import StringIO
import subprocess
from tempfile import NamedTemporaryFile
import os
import logging
import argparse
import shutil


logging.basicConfig(
    format='[%(levelname)s] %(message)s',
    level=logging.INFO
)


TEMP_FILE_PREFIX = "altaf_tmp_"

# Example mappings for inheritance interpretation
INHERITANCE_DICT_TRIO = {
    "1-0-0": "De novo",
    "1-1-0": "Maternal",
    "1-0-1": "Paternal",
    "1-1-1": "Shared"
}

INHERITANCE_DICT_DUO_MOTHER = {
    "1-1": "Shared with mother",
    "1-0": "Unique to index"
}

INHERITANCE_DICT_DUO_FATHER = {
    "1-1": "Shared with father",
    "1-0": "Unique to index"
}


def delete_vcfs(vcf_path):
    """Delete VCF and its tabix index if they exist."""
    for ext in ["", ".tbi"]:
        path = vcf_path + ext
        if os.path.exists(path):
            os.remove(path)
            logging.info(f"Deleted temporary file: {path}")


def create_vcf_tbi(vcf_file):
    """Create tabix index (.tbi) for VCF file if not present."""
    if not os.path.isfile(vcf_file + ".tbi"):
        try:
            cmd = f"tabix {vcf_file}"
            subprocess.check_output(cmd, shell=True)
            return True
        except subprocess.CalledProcessError:
            logging.warning("tabix failed — ensure VCF is bgzipped and formatted correctly.")
            return False
    return True


def collect_snv_inheritance(vcf_dict):
    """Collect SNV inheritance information depending on available family members."""
    cmd = None
    domain_setting = None

    if vcf_dict.get("index") and vcf_dict.get("mother") and vcf_dict.get("father"):
        inheritance_dict = INHERITANCE_DICT_TRIO
        cmd = f"bcftools isec -n +1 {vcf_dict['index']} {vcf_dict['mother']} {vcf_dict['father']}"
        domain_setting = "trio"

    elif vcf_dict.get("index") and vcf_dict.get("mother"):
        inheritance_dict = INHERITANCE_DICT_DUO_MOTHER
        cmd = f"bcftools isec -n +1 {vcf_dict['index']} {vcf_dict['mother']}"
        domain_setting = "duo_mat"

    elif vcf_dict.get("index") and vcf_dict.get("father"):
        inheritance_dict = INHERITANCE_DICT_DUO_FATHER
        cmd = f"bcftools isec -n +1 {vcf_dict['index']} {vcf_dict['father']}"
        domain_setting = "duo_pat"

    if not cmd:
        raise ValueError("Invalid VCF dictionary: missing index/mother/father combination.")

    logging.info(f"Running bcftools isec for {domain_setting} mode...")
    isec_result = StringIO(subprocess.check_output(cmd, shell=True).decode('utf-8'))
    df_isec = pd.read_csv(isec_result, sep="\t", header=None, dtype=str,
                          names=["chr", "pos", "ref", "alt", "snv_occurence"])
    df_isec["pos"] = df_isec["pos"].astype(int)
    df_isec["snv_occurence"] = df_isec["snv_occurence"].replace(inheritance_dict, regex=True)
    
    return df_isec, domain_setting


def read_vcf_file(vcf_file):
    """Read a VCF file using cyvcf2 and extract variant-level information."""
    if vcf_file is None or not os.path.isfile(vcf_file):
        raise FileNotFoundError("No VCF file found or path invalid.")

    df_vcf_variants = pd.DataFrame(columns=["chr", "start", "end", "altAF", "quality", "reads"])
    li_chr, li_start, li_end, li_qual, li_af, li_reads = [], [], [], [], [], []

    try:
        for variant in VCF(vcf_file):
            if not variant.ALT == []:
                ad = variant.format('AD')
                dp = variant.format('DP')
                
                if (ad is None) or (dp is None):
                    altaf = 0
                else:
                    altaf = ad[0,1] / dp[0,0]

                li_chr.append(variant.CHROM)
                li_start.append(variant.POS)
                li_end.append(variant.POS + len(variant.REF))
                li_qual.append(variant.QUAL)
                li_af.append(altaf)
                li_reads.append(dp[0,0] if dp is not None else 0)
        
        df_vcf_variants["chr"] = li_chr
        df_vcf_variants["start"] = li_start
        df_vcf_variants["end"] = li_end
        df_vcf_variants["quality"] = li_qual
        df_vcf_variants["altAF"] = li_af
        df_vcf_variants["reads"] = li_reads
        df_vcf_variants["start"] = df_vcf_variants["start"].astype(int)

        return df_vcf_variants

    except Exception as e:
        delete_vcfs(vcf_file)
        raise RuntimeError(f"VCF parsing failed: {repr(e)}") from e


def detect_roh(vcf_file):
    """Detect ROH regions using bcftools roh."""
    cmd = f"bcftools roh --AF-dflt 0.4 -I {vcf_file} | awk '$1==\"RG\"{{print $0}}'"
    try:
        bcf_roh_results = StringIO(subprocess.check_output(cmd, shell=True).decode('utf-8'))
        df_roh_rg = pd.read_csv(bcf_roh_results, sep="\t",
                                names=["RG", "sample", "chr", "start", "end", "length", "number_of_markers", "quality"])

        if df_roh_rg.empty:
            cmd = f"bcftools roh --AF-dflt 0.4 -G 30 -I {vcf_file} | awk '$1==\"RG\"{{print $0}}'"
            bcf_roh_results = StringIO(subprocess.check_output(cmd, shell=True).decode('utf-8'))
            df_roh_rg = pd.read_csv(bcf_roh_results, sep="\t",
                                    names=["RG", "sample", "chr", "start", "end", "length", "number_of_markers", "quality"])
        
        df_roh_rg["chr"] = df_roh_rg["chr"].astype(str).str.replace("chr", "")
        return df_roh_rg

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"bcftools roh failed: {repr(e)}") from e


def save_temporary_file(vcf_file_in):
    """Save uploaded file temporarily and create tabix index."""
    with NamedTemporaryFile("wb", prefix=TEMP_FILE_PREFIX, suffix=".vcf.gz", delete=False) as vcf_file:
        vcf_file.write(vcf_file_in.getvalue())

    if not create_vcf_tbi(vcf_file.name):
        delete_vcfs(vcf_file.name)
        raise RuntimeError("tabix indexing failed.")

    return vcf_file.name

def main():
    parser = argparse.ArgumentParser(description="CLI version of ALTAF VCF processing tool.")
    parser.add_argument("--vcf", required=True, help="Path to input VCF file")
    parser.add_argument("--mode", choices=["roh", "inheritance", "read"], required=True, help="Analysis mode")
    parser.add_argument("--mother", help="Path to mother's VCF (optional)")
    parser.add_argument("--father", help="Path to father's VCF (optional)")
    parser.add_argument("--output", default="output.tsv", help="Output file name")

    args = parser.parse_args()

    try:
        if args.mode == "roh":
            df_roh = detect_roh(args.vcf)
            df_roh.to_csv(args.output, sep="\t", index=False)
            logging.info(f"ROH analysis complete → {args.output}")

        elif args.mode == "inheritance":
            vcf_dict = {"index": args.vcf, "mother": args.mother, "father": args.father}
            df_inherit, mode = collect_snv_inheritance(vcf_dict)
            output_file = f"inheritance_{mode}.tsv"
            df_inherit.to_csv(output_file, sep="\t", index=False)
            logging.info(f"Inheritance analysis ({mode}) complete → {output_file}")

        elif args.mode == "read":
            df_vcf = read_vcf_file(args.vcf)
            df_vcf.to_csv(args.output, sep="\t", index=False)
            logging.info(f"VCF parsed successfully → {args.output}")

    except Exception as e:
        logging.error(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
