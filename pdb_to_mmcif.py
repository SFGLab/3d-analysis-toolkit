import time
def convert(fileName):
    header = """data_3dnome
# 
_entry.id 3dnome
# 
_audit_conform.dict_name       mmcif_pdbx.dic 
_audit_conform.dict_version    5.296 
_audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic 
#
loop_
_atom_site.group_PDB 
_atom_site.id 
_atom_site.type_symbol 
_atom_site.label_atom_id 
_atom_site.label_alt_id 
_atom_site.label_comp_id 
_atom_site.label_asym_id 
_atom_site.label_entity_id 
_atom_site.label_seq_id 
_atom_site.pdbx_PDB_ins_code 
_atom_site.Cartn_x 
_atom_site.Cartn_y 
_atom_site.Cartn_z 
_atom_site.occupancy 
_atom_site.B_iso_or_equiv 
_atom_site.auth_asym_id
"""
    atoms = list()
    atoms.append(header)
    with open(fileName, 'r') as f: #open the file
        lines = f.readlines()
        for line in lines:
            values = line.split()
            i = 1
            if(values[4] != "A"):
                i = 0
                num = values[4].split("A")[1]
            else:
                num = values[5]
            
            new_line = "ATOM " + values[1] + " C " + values[2] + " . " + values[3] + " A  1 " + num + " ? " + values[5+i] + " " + values[6+i] + " " \
            + values[7+i] + " " + values[8+i] + " " + values[9+i] + " A\n"
            atoms.append(new_line)

    with open("conv.mmcif", 'w') as f:
        f.writelines(atoms)

def main():
    start_time_total = time.time()
    convert("/mnt/raid/ctcf_prediction_anal/converter_hcm_pdb/hcm_files/top100k_chr8.pdb")
    print("--- Executed in %s seconds ---" % (time.time() - start_time_total))

if __name__ == "__main__":
    main()
