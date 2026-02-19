import os
import re
import csv
import tempfile
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager

# --- CONFIGURATION ---
BASE_DIR = "folder/"              
ORIG_DIR = "folder/"              
OUT_SUFFIX = "_out"               
DRUG_THRESHOLD = 0.5              
VOL_LIMITS = (400, 2000)          
STATS_FILE = "pocket_screening_master_data.csv"
READY_DIR = os.path.join(BASE_DIR, "APoc_Ready_Pockets")
CORES = 100  
PROGRESS_INTERVAL = 1000 # Print update every 1000 proteins

def process_single_protein(folder_name, counter, total):
    """Function executed by each core."""
    prot_name = folder_name.replace(OUT_SUFFIX, "")
    info_path = os.path.join(BASE_DIR, folder_name, f"{prot_name}_info.txt")
    orig_pdb_path = os.path.join(ORIG_DIR, f"{prot_name}.pdb")
    
    # Update progress safely across cores
    with counter.get_lock():
        counter.value += 1
        if counter.value % PROGRESS_INTERVAL == 0:
            print(f" [Progress] {counter.value}/{total} proteins processed ({(counter.value/total)*100:.1f}%)")

    if not os.path.exists(info_path) or not os.path.exists(orig_pdb_path):
        return None

    pocket_data_list = []
    try:
        with open(orig_pdb_path, 'r') as f:
            orig_lines = [line for line in f if line.startswith("ATOM")]

        with open(info_path, 'r') as f:
            content = f.read()
            pockets = content.split("Pocket ")[1:]
            
            for p_data in pockets:
                try:
                    d_score = float(re.search(r"Druggability Score\s*:\s*([\d.]+)", p_data).group(1))
                    vol = float(re.search(r"Volume\s*:\s*([\d.]+)", p_data).group(1))
                    
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
                        
                        with tempfile.NamedTemporaryFile('w', delete=False, dir=READY_DIR) as tf:
                            for line in orig_lines:
                                if (line[21], line[22:26].strip()) in res_to_keep:
                                    tf.write(line)
                            tf.write("TER\nEND\n")
                            temp_name = tf.name
                        
                        os.replace(temp_name, final_path)
                        pocket_data_list.append([prot_name, p_num, d_score, vol])
                except: continue
            
            return pocket_data_list
    except:
        return None

if __name__ == "__main__":
    os.makedirs(READY_DIR, exist_ok=True)
    all_folders = [f for f in os.listdir(BASE_DIR) if f.endswith(OUT_SUFFIX)]
    total_folders = len(all_folders)
    
    # Manager for cross-process communication
    manager = Manager()
    counter = manager.Value('i', 0)
    
    print(f"--- HPC NODE INITIALIZED ---")
    print(f"Target count: {total_folders} proteins")
    print(f"Cores: {CORES} | Filter: Drug Score > {DRUG_THRESHOLD}, Volume {VOL_LIMITS}")
    print(f"-----------------------------")

    final_csv_rows = []
    with ProcessPoolExecutor(max_workers=CORES) as executor:
        # Wrap the function to include the counter
        results = list(executor.map(process_single_protein, all_folders, 
                                    [counter]*total_folders, 
                                    [total_folders]*total_folders))
        
        for res in results:
            if res: final_csv_rows.extend(res)

    # Save Stats
    with open(STATS_FILE, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Protein_ID", "Pocket_Number", "Druggability_Score", "Volume_A3"])
        writer.writerows(final_csv_rows)

    print(f"\n--- TASK COMPLETE ---")
    print(f"Total pockets extracted: {len(final_csv_rows)}")
    print(f"Master data saved to: {STATS_FILE}")
