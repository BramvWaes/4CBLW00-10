import os
import json
from rdkit import Chem
from rdkit.Chem import Fragments

# 1. Read the compound_map.json file into a dictionary
with open("compound_map.json", "r", encoding="utf-8") as f:
    # Remove the first line if it is a comment (// ...)
    lines = f.readlines()
    if lines[0].strip().startswith("//"):
        lines = lines[1:]
    compound_map = json.loads("".join(lines))

functional_group_smarts = {
    # "Alkane":           ["[CX4]", "[CX4;H0,H1,H2]"],
    #"Alkene": "[$([CX2]=[CX2])]",
    "Alkene": "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
    # "Alkyne":             "[$([CX2]#C)]",
    "Arene": "[c]",
    "Ketone": "[#6][CX3](=O)[#6]",
    "Ester": "[#6][CX3](=O)[OX2H0][#6]",
    "Amide": "[NX3][CX3](=O)[#6]",
    "Carboxylic acid": "[CX3](=O)[OX2H1]",
    # "Alcohol":            "[CHX4][OX2H]",
    # "Amine":     "[NX3;H2,H1;!$(NC=O)]",
    # "Nitrile":           "[NX1]#[CX2]",
    # "Alkyl halide":      "[CX4][F,Cl,Br,I]",
    # "Acyl halide":      "[CX3](=O)[F,Cl,Br,I]",
    # "Ether":            "[OD2]([#6])[#6]",
    # "Nitro":             "[$([NX3](==O)==O),$([NX3+](==O)[O-])][!#8]",
    # "Methyl":            "[CH3X4]",
}

# Compile those SMARTS into RDKit Mol objects:
compiled = {}
for name, smarts in functional_group_smarts.items():
    # Nitro has two SMARTS, so handle list vs. single string:
    if isinstance(smarts, (list, tuple)):
        compiled[name] = [Chem.MolFromSmarts(s) for s in smarts]
    else:
        compiled[name] = [Chem.MolFromSmarts(smarts)]

# Read all file names in data/IR and extract ids
ir_folder = "data/IR"
ids = {}
titles = {}
fnames = {}
if os.path.isdir(ir_folder):
    for fname in os.listdir(ir_folder):
        if fname.endswith(".jdx") and "_" in fname:
            id_part = fname.split("_", 1)[0]
            ids[id_part] = fname
            fnames[id_part] = fname
            # Read the first line to get the title
            with open(
                os.path.join(ir_folder, fname), encoding="ascii", errors="ignore"
            ) as f:
                first_line = f.readline().strip()
                if first_line.startswith("##TITLE="):
                    title = first_line[len("##TITLE=") :]
                    titles[id_part] = title

# Point to your .sdf file
supplier = Chem.SDMolSupplier("data/nist_mol3D.sdf", removeHs=False)

# supplier is a lazy sequence of Mol or None
mols = [m for m in supplier if m is not None]

results = []

no_ir = 0
missing_cas = []
missing_ir = []

for mol in mols:
    if mol.HasProp("WEBBOOK.ID") == 0:
        casNumber = mol.GetProp("CAS.NUMBER")
        if casNumber in compound_map:
            webbook_id = compound_map[casNumber]
            mol.SetProp("WEBBOOK.ID", webbook_id)
        else:
            missing_cas.append(casNumber)
            continue

    webbook_id = mol.GetProp("WEBBOOK.ID")
    if webbook_id not in ids:
        no_ir += 1
        missing_ir.append(webbook_id)
        continue

    # 3) Extract the groups
    counts = {
        # ##########
        # # Oxygens
        # ##########
        # "Al_OH": Fragments.fr_Al_OH(mol),
        # "Ar_OH": Fragments.fr_Ar_OH(mol),
        # "methoxy": Fragments.fr_methoxy(mol),
        # "oxime": Fragments.fr_oxime(mol),
        # "ester": Fragments.fr_ester(mol),
        # "Al_COO": Fragments.fr_Al_COO(mol),
        # "Ar_COO": Fragments.fr_Ar_COO(mol),
        # "COO": Fragments.fr_COO(mol),
        # "COO2": Fragments.fr_COO2(mol),
        # "ketone": Fragments.fr_ketone(mol),
        # "ether": Fragments.fr_ether(mol),
        # "phenol": Fragments.fr_phenol(mol),
        # "aldehyde": Fragments.fr_aldehyde(mol),
        # ############
        # # Nitrogens
        # ############
        # "quatN": Fragments.fr_quatN(mol),
        # "NH2": Fragments.fr_NH2(mol),
        # "NH1": Fragments.fr_NH1(mol),
        # "NH0": Fragments.fr_NH0(mol),
        # "Ar_N": Fragments.fr_Ar_N(mol),
        # "Ar_NH": Fragments.fr_Ar_NH(mol),
        # "aniline": Fragments.fr_aniline(mol),
        # "Imine": Fragments.fr_Imine(mol),
        # "nitrile": Fragments.fr_nitrile(mol),
        # "hdrzine": Fragments.fr_hdrzine(mol),
        # "hdrzone": Fragments.fr_hdrzone(mol),
        # "nitroso": Fragments.fr_nitroso(mol),
        # "nitro": Fragments.fr_nitro(mol),
        # "azo": Fragments.fr_azo(mol),
        # "diazo": Fragments.fr_diazo(mol),
        # "azide": Fragments.fr_azide(mol),
        # "amide": Fragments.fr_amide(mol),
        # "priamide": Fragments.fr_priamide(mol),
        # "amidine": Fragments.fr_amidine(mol),
        # "guanido": Fragments.fr_guanido(mol),
        # "Nhpyrrole": Fragments.fr_Nhpyrrole(mol),
        # "imide": Fragments.fr_imide(mol),
        # "isocyan": Fragments.fr_isocyan(mol),
        # "isothiocyan": Fragments.fr_isothiocyan(mol),
        # "thiocyan": Fragments.fr_thiocyan(mol),
        # ###########
        # # Halogens
        # ###########
        # "halogen": Fragments.fr_halogen(mol),
        # "alkyl_halide": Fragments.fr_alkyl_halide(mol),
        # ##########
        # # Sulfurs
        # ##########
        # "sulfide": Fragments.fr_sulfide(mol),
        # "SH": Fragments.fr_SH(mol),
        # "sulfone": Fragments.fr_sulfone(mol),
        # "sulfonamd": Fragments.fr_sulfonamd(mol),
        # "prisulfonamd": Fragments.fr_prisulfonamd(mol),
        # ##################################
        # # Miscellaneous Functional Groups
        # ##################################
        # "barbitur": Fragments.fr_barbitur(mol),
        # "urea": Fragments.fr_urea(mol),
        # "term_acetylene": Fragments.fr_term_acetylene(mol),
        # "imidazole": Fragments.fr_imidazole(mol),
        # "furan": Fragments.fr_furan(mol),
        # "thiophene": Fragments.fr_thiophene(mol),
        # "thiazole": Fragments.fr_thiazole(mol),
        # "oxazole": Fragments.fr_oxazole(mol),
        # "pyridine": Fragments.fr_pyridine(mol),
        # "piperdine": Fragments.fr_piperdine(mol),
        # "piperzine": Fragments.fr_piperzine(mol),
        # "morpholine": Fragments.fr_morpholine(mol),
        # "lactam": Fragments.fr_lactam(mol),
        # "lactone": Fragments.fr_lactone(mol),
        # "tetrazole": Fragments.fr_tetrazole(mol),
        # "epoxide": Fragments.fr_epoxide(mol),
        # "unbrch_alkane": Fragments.fr_unbrch_alkane(mol),
        # "bicyclic": Fragments.fr_bicyclic(mol),
        # "benzene": Fragments.fr_benzene(mol),
        # #############
        # # Phosphates
        # #############
        # "phos_acid": Fragments.fr_phos_acid(mol),
        # "phos_ester": Fragments.fr_phos_ester(mol),
        # #####################
        # # Topliss Metabolism
        # #####################
        # "nitro_arom": Fragments.fr_nitro_arom(mol),
        # "nitro_arom_nonortho": Fragments.fr_nitro_arom_nonortho(mol),
        # "dihydropyridine": Fragments.fr_dihydropyridine(mol),
        # "phenol_noOrthoHbond": Fragments.fr_phenol_noOrthoHbond(mol),
        # "Al_OH_noTert": Fragments.fr_Al_OH_noTert(mol),
        # "benzodiazepine": Fragments.fr_benzodiazepine(mol),
        # "para_hydroxylation": Fragments.fr_para_hydroxylation(mol),
        # "allylic_oxid": Fragments.fr_allylic_oxid(mol),
        # "aryl_methyl": Fragments.fr_aryl_methyl(mol),
        # "Ndealkylation1": Fragments.fr_Ndealkylation1(mol),
        # "Ndealkylation2": Fragments.fr_Ndealkylation2(mol),
        # "alkyl_carbamate": Fragments.fr_alkyl_carbamate(mol),
        # "ketone_Topliss": Fragments.fr_ketone_Topliss(mol),
        # "ArN": Fragments.fr_ArN(mol),
        # "HOCCN": Fragments.fr_HOCCN(mol),
        # "Alkene": 1 if mol.HasSubstructMatch(alkene_pat) else 0,
    }
    for name, patterns in compiled.items():
        counts[name] = 0
        for pattern in patterns:
            counts[name] += mol.HasSubstructMatch(pattern)

    #if any(count > 0 for count in counts.values()):
    result = {
        "webbook_id": webbook_id,
        "title": titles.get(webbook_id, ""),
        "filename": fnames.get(webbook_id, ""),
        "groups": counts,
    }
    results.append(result)

with open("groups.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2)

with open("missing_cas.json", "w", encoding="utf-8") as f:
    json.dump(missing_cas, f, indent=2)

with open("missing_ir.json", "w", encoding="utf-8") as f:
    json.dump(missing_ir, f, indent=2)


print(f"Read {len(mols)} valid molecules out of {len(supplier)} entries")
print(f"Found {no_ir} molecules without matching IR file")
print(f"Found {len(missing_cas)} without matched CAS number")
print(f"Found {len(results)} molecules with functional groups")
