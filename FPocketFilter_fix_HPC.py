import os
import re
import csv
import tempfile
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager

# --- CONFIGURATION ---
BASE_DIR = "mybigassstorage/fpocketscreen/"                                   # Root directory for fpocket output folders. A.K.A where is your "protein_out/ directory?"       
ORIG_DIR = "mybigassstorage/fpocketscreen/proteindataset/"                    # your original pdb file directory, A.K.A where is your residue/pocket corrected, APoc compatible pdb files?              
OUT_SUFFIX = "_out"                                                           # if you changed the output suffix during Fpocket run, please indicate here.
DRUG_THRESHOLD = 0.5                                                          # druggability score cutoff (0 - 1.0). Pockets below this value will be ignored, regardless of VOLUME.          
VOL_LIMITS = (400, 2000)                                                      # Pocket volume acceptence interval. (int- int). Pockets outside of the range will be ignored regardless of DRUG THRESHOLD.
STATS_FILE = "pocket_screening_master_data.csv"                               # Name of the statistics results CSV.
READY_DIR = os.path.join(BASE_DIR, "APoc_Ready_Pockets")                      # Filtered and accepted pockets (pdb files) will be placed here under BASE_DIR
CORES = 4                                                                     # Number of cores (or instances if $nproc < CORES)  to use during run. Default is 4.   
PROGRESS_INTERVAL = 1000                                                      # # of files processed until next progress output

def process_single_protein(folder_name, counter, total):
    prot_name = folder_name.replace(OUT_SUFFIX, "")
    info_path = os.path.join(BASE_DIR, folder_name, f"{prot_name}_info.txt")
    orig_pdb_path = os.path.join(ORIG_DIR, f"{prot_name}.pdb")
    
    with counter.get_lock():
        counter.value += 1
        if counter.value % PROGRESS_INTERVAL == 0:
            print(f" [Progress] {counter.value}/{total} proteins ({(counter.value/total)*100:.1f}%)")

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
            high_druggable_count = 0
            
            # Temporary storage to count high druggable pockets first
            temp_pockets = []
            for p_data in pockets:
                try:
                    d_score = float(re.search(r"Druggability Score\s*:\s*([\d.]+)", p_data).group(1))
                    g_score = float(re.search(r"Score\s*:\s*([\d.]+)", p_data).group(1))
                    vol = float(re.search(r"Volume\s*:\s*([\d.]+)", p_data).group(1))
                    p_num = re.match(r"(\d+)", p_data.strip()).group(1)
                    
                    is_high_druggable = d_score >= DRUG_THRESHOLD
                    if is_high_druggable:
                        high_druggable_count += 1
                    
                    # Only keep if passes screening filters
                    if is_high_druggable and VOL_LIMITS[0] <= vol <= VOL_LIMITS[1]:
                        temp_pockets.append((p_num, d_score, g_score, vol))
                except: continue

            # Second pass: write files and finalize rows with total counts
            for p_num, d_score, g_score, vol in temp_pockets:
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
                
                results_for_protein.append([
                    prot_name, d_score, g_score, vol, p_num, high_druggable_count, total_pockets
                ])
            
            return results_for_protein
    except:
        return None

if __name__ == "__main__":
    os.makedirs(READY_DIR, exist_ok=True)
    all_folders = [f for f in os.listdir(BASE_DIR) if f.endswith(OUT_SUFFIX)]
    manager = Manager()
    counter = manager.Value('i', 0)
    
    final_csv_rows = []
    with ProcessPoolExecutor(max_workers=CORES) as executor:
        results = list(executor.map(process_single_protein, all_folders, [counter]*len(all_folders), [len(all_folders)]*len(all_folders)))
        for res in results:
            if res: final_csv_rows.extend(res)

    with open(STATS_FILE, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Protein_ID", "Druggability_Score", "General_Score", "Volume", "Pocket_No", "High_Druggable_Count", "Total_Pockets"])
        writer.writerows(final_csv_rows)
    print(f"DONE. Output in {STATS_FILE}")
