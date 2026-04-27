import streamlit as st
import csv
import itertools
from Bio import Restriction
from Bio.Seq import Seq

# ==========================================
# 1. PAGE CONFIG & STYLING
# ==========================================
st.set_page_config(page_title="Plasmid Shredder", page_icon="🧬", layout="wide")

st.title("🧬 Plasmid Shredder & Cocktail Optimizer")
st.markdown("""
Upload your sequence, set your constraints, and find the most compatible NEB enzyme cocktails 
to get the fragment sizes you need.
""")

# ==========================================
# 2. HELPER FUNCTIONS (Your Logic)
# ==========================================
@st.cache_data # This prevents re-loading the CSV every time you click a button
def load_neb_data(filepath):
    db = {}
    try:
        with open(filepath, 'r', encoding='cp1252', errors='replace') as f:
            reader = csv.reader(f)
            next(reader); next(reader)
            for row in reader:
                if not row or not row[0].strip(): continue
                clean_name = row[0].strip().replace('®', '').replace('™', '').strip()
                buffers = {'r1.1': row[4], 'r2.1': row[5], 'r3.1': row[6], 'rCutSmart': row[7]}
                parsed_buffers = {k: {'activity': int(v.replace('*','').replace('<','')) if v.strip().replace('*','').replace('<','').isdigit() else 0, 
                                     'star': '*' in v} for k, v in buffers.items()}
                db[clean_name] = {'buffers': parsed_buffers, 'temp': row[9].strip()}
        return db
    except FileNotFoundError:
        return None

def calculate_fragments(all_cuts, seq_len):
    if not all_cuts: return []
    if len(all_cuts) == 1: return [seq_len]
    fragments = [all_cuts[i+1] - all_cuts[i] for i in range(len(all_cuts) - 1)]
    fragments.append((seq_len - all_cuts[-1]) + all_cuts[0])
    return sorted(fragments)

def score_cocktail(combo, db):
    missing = [enz for enz in combo if enz not in db]
    if missing: return -999, "Missing Data"
    
    temps = set(db[enz]['temp'] for enz in combo)
    best_score = -100
    best_buf = ""
    
    for b in ['r1.1', 'r2.1', 'r3.1', 'rCutSmart']:
        activity = min(db[enz]['buffers'][b]['activity'] for enz in combo)
        star = any(db[enz]['buffers'][b]['star'] for enz in combo)
        current_score = activity - (10 if star else 0)
        if current_score > best_score:
            best_score, best_buf, star_risk = current_score, b, star

    score = 100 - (50 if len(temps) > 1 else 0) - (100 - min(db[enz]['buffers'][best_buf]['activity'] for enz in combo)) - (10 if star_risk else 0)
    info = f"{best_buf} ({min(db[enz]['buffers'][best_buf]['activity'] for enz in combo)}% Act), {', '.join(temps)}"
    return score, info

# ==========================================
# 3. SIDEBAR INPUTS
# ==========================================
with st.sidebar:
    st.header("Parameters")
    min_size = st.number_input("Min Fragment Size (bp)", value=30)
    max_size = st.number_input("Max Fragment Size (bp)", value=500)
    max_enzymes = st.slider("Max Enzymes in Cocktail", 1, 3, 2)
    
    st.divider()
    st.info("Ensure 'neb_buffers.csv' is in the app directory.")

# ==========================================
# 4. MAIN INTERFACE
# ==========================================
raw_seq = st.text_area("Paste Plasmid Sequence", height=200, placeholder="ATGC...")

if st.button("Run Shredder Analysis"):
    if not raw_seq:
        st.error("Please provide a sequence!")
    else:
        neb_db = load_neb_data('neb_buffers.csv')
        if not neb_db:
            st.error("CSV file not found!")
        else:
            clean_seq = "".join(raw_seq.split()).upper()
            seq_obj = Seq(clean_seq)
            seq_len = len(clean_seq)
            
            # Pre-compute cuts
            neb_batch = Restriction.RestrictionBatch(first=[], suppliers=['N'])
            enzyme_cuts = {str(e): e.search(seq_obj, linear=False) for e in neb_batch if e.search(seq_obj, linear=False)}
            active_enzymes = list(enzyme_cuts.keys())
            
            results = []
            with st.spinner("Crunching combinations..."):
                for r in range(1, max_enzymes + 1):
                    for combo in itertools.combinations(active_enzymes, r):
                        all_cuts = sorted(list(set().union(*(enzyme_cuts[e] for e in combo))))
                        frags = calculate_fragments(all_cuts, seq_len)
                        
                        if frags and all(min_size <= f <= max_size for f in frags):
                            score, info = score_cocktail(combo, neb_db)
                            results.append({
                                "Score": score,
                                "Cocktail": ", ".join(combo),
                                "Fragments (bp)": str(frags),
                                "Compatibility Info": info
                            })
            
            if results:
                st.success(f"Found {len(results)} valid combinations!")
                # Display results in a nice table
                import pandas as pd
                df = pd.DataFrame(results).sort_values(by="Score", ascending=False)
                st.dataframe(df, use_container_width=True, hide_index=True)
            else:
                st.warning("No combinations found matching those size criteria.")
