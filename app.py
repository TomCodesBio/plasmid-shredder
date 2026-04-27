import streamlit as st
import csv
import itertools
from Bio import Restriction
from Bio.Seq import Seq
import pandas as pd

# ==========================================
# 1. PAGE CONFIG & STYLING
# ==========================================
st.set_page_config(page_title="Plasmid Shredder", page_icon="🧬", layout="wide")

st.title("🧬 Plasmid Shredder & Cocktail Optimiser")
st.markdown("""
Upload your sequence, set your constraints, and find the most compatible **NEB enzyme cocktails** to get the fragment sizes you need. You can now also pre-linearize your plasmid using an enzyme from **any supplier**.
""")

# ==========================================
# 2. HELPER FUNCTIONS
# ==========================================
@st.cache_data
def load_neb_data(filepath):
    db = {}
    try:
        with open(filepath, 'r', encoding='cp1252', errors='replace') as f:
            reader = csv.reader(f)
            next(reader); next(reader)
            for row in reader:
                if not row or not row[0].strip(): continue
                clean_name = row[0].strip().replace('®', '').replace('™', '').strip()
                try:
                    buffers = {'r1.1': row[4], 'r2.1': row[5], 'r3.1': row[6], 'rCutSmart': row[7]}
                    parsed_buffers = {}
                    for b, val in buffers.items():
                        v_str = val.strip()
                        num_val = v_str.replace('*','').replace('<','')
                        parsed_buffers[b] = {
                            'activity': int(num_val) if num_val.isdigit() else 0,
                            'star': '*' in v_str
                        }
                    db[clean_name] = {'buffers': parsed_buffers, 'temp': row[9].strip()}
                except IndexError:
                    pass
        return db
    except FileNotFoundError:
        return None

def calculate_fragments(all_cuts, seq_len, is_circular):
    if not all_cuts: 
        return [seq_len] if not is_circular else []
        
    fragments = []
    if is_circular:
        if len(all_cuts) == 1: return [seq_len]
        for i in range(len(all_cuts) - 1):
            fragments.append(all_cuts[i+1] - all_cuts[i])
        fragments.append((seq_len - all_cuts[-1]) + all_cuts[0])
    else: 
        prev_cut = 0 
        for cut in all_cuts:
            fragments.append(cut - prev_cut)
            prev_cut = cut
        fragments.append(seq_len - prev_cut)
        
    return sorted(fragments)

def score_cocktail(combo, db):
    missing = [enz for enz in combo if enz not in db]
    if missing: return -999, "Missing Data"
    
    temps = set(db[enz]['temp'] for enz in combo)
    best_score, best_buf, star_risk = -100, "", False
    
    for b in ['r1.1', 'r2.1', 'r3.1', 'rCutSmart']:
        activity = min(db[enz]['buffers'][b]['activity'] for enz in combo)
        star = any(db[enz]['buffers'][b]['star'] for enz in combo)
        current_score = activity - (10 if star else 0)
        if current_score > best_score:
            best_score, best_buf, star_risk = current_score, b, star

    min_best_act = min(db[enz]['buffers'][best_buf]['activity'] for enz in combo)
    score = 100 - (50 if len(temps) > 1 else 0) - (100 - min_best_act) - (10 if star_risk else 0)
    info = f"{best_buf} ({min_best_act}% Act), Temp: {', '.join(temps)}"
    if star_risk: info += " ⚠️ Star Risk"
    return score, info

# ==========================================
# 3. SIDEBAR INPUTS
# ==========================================
with st.sidebar:
    st.header("1. Fragment Constraints")
    min_size = st.number_input("Min Fragment Size (bp)", value=30)
    max_size = st.number_input("Max Fragment Size (bp)", value=500)
    max_enzymes = st.slider("Max Enzymes in Cocktail", 1, 3, 2)
    
    st.divider()
    
    st.header("2. Topology / Linearization")
    is_linearized = st.checkbox("Plasmid will be pre-linearized")
    lin_enzyme_name = None
    if is_linearized:
        lin_enzyme_name = st.text_input("Linearizing Enzyme Name", value="XhoI", help="Must be an exact commercial name (e.g., EcoRI). Case-sensitive.")
    
    st.divider()
    st.info("Ensure 'neb_buffers.csv' is in the app directory.")

# ==========================================
# 4. MAIN INTERFACE
# ==========================================
raw_seq = st.text_area("Paste Plasmid Sequence", height=200, placeholder="ATGC...")

if st.button("Run Shredder Analysis", type="primary"):
    if not raw_seq:
        st.error("Please provide a sequence!")
        st.stop()
        
    neb_db = load_neb_data('neb_buffers.csv')
    if not neb_db:
        st.error("CSV file not found! Make sure 'neb_buffers.csv' is uploaded.")
        st.stop()

    clean_seq = "".join(raw_seq.split()).upper()
    seq_obj = Seq(clean_seq)
    seq_len = len(clean_seq)
    is_circular = not is_linearized

    # --- TOPOLOGY / LINEARIZATION LOGIC ---
    if is_linearized:
        if not lin_enzyme_name:
            st.error("Please enter a linearizing enzyme.")
            st.stop()
            
        try:
            # Check entire commercial database
            lin_enz = getattr(Restriction, lin_enzyme_name)
        except AttributeError:
            st.error(f"❌ '{lin_enzyme_name}' not found in the Biopython database. Check spelling/capitalization.")
            st.stop()

        cuts = lin_enz.search(seq_obj, linear=False)
        if len(cuts) == 0:
            st.error(f"❌ {lin_enzyme_name} does not cut this plasmid.")
            st.stop()
        elif len(cuts) > 1:
            st.error(f"❌ {lin_enzyme_name} cuts {len(cuts)} times! A linearizing enzyme must cut exactly once.")
            st.stop()
        
        # Shift sequence so cut site is at index 0
        cut_pos = cuts[0]
        clean_seq = clean_seq[cut_pos-1:] + clean_seq[:cut_pos-1]
        seq_obj = Seq(clean_seq)
        st.success(f"✅ Plasmid successfully linearized at position {cut_pos} using {lin_enzyme_name}.")

    # --- SHREDDER ANALYSIS ---
    with st.spinner("Crunching NEB combinations..."):
        # Load ONLY NEB enzymes for the shredding cocktails
        neb_batch = Restriction.RestrictionBatch(first=[], suppliers=['N'])
        
        # Find cuts for all NEB enzymes on the (possibly shifted) sequence
        enzyme_cuts = {}
        for enz in neb_batch:
            enz_str = str(enz)
            if is_linearized and enz_str == lin_enzyme_name:
                continue # Don't use the linearizing enzyme in the shredder cocktail
                
            cuts = enz.search(seq_obj, linear=not is_circular)
            if cuts:
                enzyme_cuts[enz_str] = cuts
                
        active_enzymes = list(enzyme_cuts.keys())
        results = []
        
        for r in range(1, max_enzymes + 1):
            for combo in itertools.combinations(active_enzymes, r):
                all_cuts = sorted(list(set().union(*(enzyme_cuts[e] for e in combo))))
                frags = calculate_fragments(all_cuts, seq_len, is_circular)
                
                if frags and all(min_size <= f <= max_size for f in frags):
                    score, info = score_cocktail(combo, neb_db)
                    results.append({
                        "Score": score,
                        "Cocktail": ", ".join(combo),
                        "Fragments (bp)": str(frags),
                        "Compatibility Info": info
                    })
    
    # --- DISPLAY RESULTS ---
    if results:
        st.success(f"Found {len(results)} valid NEB combinations!")
        df = pd.DataFrame(results).sort_values(by="Score", ascending=False)
        st.dataframe(df, use_container_width=True, hide_index=True)
    else:
        st.warning("No NEB combinations found matching those size criteria.")
