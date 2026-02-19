import os
import re
import csv
import tempfile
from concurrent.futures import ProcessPoolExecutor

# --- CONFIGURATION ---
BASE_DIR = "folder/"                                            # Root directory for fpocket output folders. A.K.A where is your "protein_out/ directory?"                     
ORIG_DIR = "folder/"                                            # your original pdb file directory, A.K.A where is your residue/pocket corrected, APoc compatible pdb files?                        
OUT_SUFFIX = "_out"                                             # if you changed the output suffix during Fpocket run, please indicate here.    
DRUG_THRESHOLD = 0.5                                            # druggability score cutoff (0 - 1.0). Pockets below this value will be ignored, regardless of VOLUME.          
VOL_LIMITS = (400, 2500)                                        # Pocket volume acceptence interval. (int- int). Pockets outside of the range will be ignored regardless of DRUG THRESHOLD.
STATS_FILE = "pocket_screening_master_data.csv"                 # Name of the statistics results CSV.
READY_DIR = os.path.join(BASE_DIR, "APoc_Ready_Pockets")        # Filtered and accepted pockets (pdb files) will be placed here under BASE_DIR
CORES = 4                                                       # Number of cores (or instances if $nproc < CORES)  to use during run. Default is 4.   

def process_single_protein(folder_name):
    """Simplified worker: No shared counter to avoid AttributeError."""
    prot_name = folder_name.replace(OUT_SUFFIX, "")
    info_path = os.path.join(BASE_DIR, folder_name, f"{prot_name}_info.txt")
    orig_pdb_path = os.path.join(ORIG_DIR, f"{prot_name}.pdb")
    
    if not os.path.exists(info_path) or not os.path.exists(orig_pdb_path):
        return None

    results_for_protein = []
    try:
        with open(orig_pdb_path, 'r') as f:
            orig_lines = [line for line in f if line.startswith("ATOM")]

        with open(info_path, 'r') as f:
            content = f.read()
            pockets = content.split("Pocket ")[1:]
            
            total_pockets = len(pockets)
            # Efficiently count high druggable pockets
            high_druggable_count = 0
            for p in pockets:
                d_match = re.search(r"Druggability Score\s*:\s*([\d.]+)", p)
                if d_match and float(d_match.group(1)) >= DRUG_THRESHOLD:
                    high_druggable_count += 1
            
            for p_data in pockets:
                try:
                    d_match = re.search(r"Druggability Score\s*:\s*([\d.]+)", p_data)
                    g_match = re.search(r"Score\s*:\s*([\d.]+)", p_data)
                    v_match = re.search(r"Volume\s*:\s*([\d.]+)", p_data)
                    
                    if not (d_match and g_match and v_match): continue
                    
                    d_score = float(d_match.group(1))
                    g_score = float(g_match.group(1))
                    vol = float(v_match.group(1))
                    
                    if d_score >= DRUG_THRESHOLD and VOL_LIMITS[0] <= vol <= VOL_LIMITS[1]:
                        p_num = re.match(r"(\d+)", p_data.strip()).group(1)
                        pocket_pdb = os.path.join(BASE_DIR, folder_name, "pockets", f"pocket{p_num}_atm.pdb")
                        
                        if not os.path.exists(pocket_pdb): continue
                        
                        res_to_keep = set()
                        with open(pocket_pdb, 'r') as pf:
                            for line in pf:
                                if line.startswith("ATOM"):
                                    res_to_keep.add((line[21], line[22:26].strip()))
                        
                        out_name = f"{prot_name}_p{p_num}_fixed.pdb"
                        final_path = os.path.join(READY_DIR, out_name)
                        
                        # Atomic write to prevent corruption
                        with tempfile.NamedTemporaryFile('w', delete=False, dir=READY_DIR) as tf:
                            for line in orig_lines:
                                if (line[21], line[22:26].strip()) in res_to_keep:
                                    tf.write(line)
                            tf.write("TER\nEND\n")
                            temp_name = tf.name
                        os.replace(temp_name, final_path)
                        
                        results_for_protein.append([
                            prot_name, d_score, g_score, vol, p_num, high_druggable_count, total_pockets
                        ])
                except: continue
            return results_for_protein
    except:
        return None

if __name__ == "__main__":
    os.makedirs(READY_DIR, exist_ok=True)
    all_folders = [f for f in os.listdir(BASE_DIR) if f.endswith(OUT_SUFFIX)]
    total_len = len(all_folders)
    
    print(f"--- HPC NODE STARTING ---")
    print(f"Proteins to process: {total_len}")
    print(f"Cores allocated: {CORES}")
    print(f"Output Directory: {READY_DIR}")

    final_csv_rows = []
    
    # Process everything. 
    # Python 3.9's ProcessPoolExecutor works best when kept simple.
    with ProcessPoolExecutor(max_workers=CORES) as executor:
        # We use a list to force the execution and allow us to track chunks
        chunk_size = 1000
        for i in range(0, total_len, chunk_size):
            chunk = all_folders[i:i + chunk_size]
            results = list(executor.map(process_single_protein, chunk))
            
            for res in results:
                if res: final_csv_rows.extend(res)
            
            print(f" > Processed {min(i + chunk_size, total_len)}/{total_len} proteins...")

    # Write the Master CSV
    with open(STATS_FILE, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Protein_ID", "Druggability_Score", "General_Score", "Volume", "Pocket_No", "High_Druggable_Count", "Total_Pockets"])
        writer.writerows(final_csv_rows)
        
    print(f"\n--- SUCCESS ---")
    print(f"Total pockets extracted: {len(final_csv_rows)}")
    print(f"Data saved to: {STATS_FILE}")
