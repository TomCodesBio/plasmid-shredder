import streamlit as st
import csv
import itertools
from Bio import Restriction
from Bio.Seq import Seq
import pandas as pd

# ==========================================
# 1. PAGE CONFIG & STYLING
# ==========================================
st.set_page_config(page_title="Plasmid Shredder Pro", page_icon="🧬", layout="wide")

st.title("🧬 Plasmid Shredder: Mass Analytics Edition")
st.markdown("""
Find the best NEB enzyme cocktails to fragment your plasmid, and calculate high-precision 
monoisotopic and average masses for both strands of the resulting fragments.
""")

# ==========================================
# 2. ANALYTICAL FUNCTIONS (High Precision)
# ==========================================
def calculate_dna_masses(sequence_obj):
    """
    Calculates precise Average and Monoisotopic masses.
    Uses standard DNA residue masses.
    A restriction fragment has a 5' Phosphate and 3' OH. 
    Because standard residue masses already contain 1 phosphate per base,
    summing them and adding H2O perfectly models this structure!
    """
    seq_str = str(sequence_obj).upper()
    
    # 1. EXACT Monoisotopic masses of DNA residues (already includes 1 phosphate each)
    mono_residue = {
        'A': 313.0576,
        'C': 289.0464,
        'G': 329.0525,
        'T': 304.0460
    }
    
    # 2. AVERAGE masses of DNA residues
    avg_residue = {
        'A': 313.21,
        'C': 289.18,
        'G': 329.21,
        'T': 304.20
    }
    
    total_mono = sum(mono_residue.get(base, 0) for base in seq_str)
    total_avg = sum(avg_residue.get(base, 0) for base in seq_str)
    
    if len(seq_str) > 0:
        # Add H2O to cap the 5' and 3' ends
        # Since the residues already provide exactly N phosphates, this yields a 5'-P, 3'-OH fragment.
        total_mono += 18.0106
        total_avg += 18.0153
        
    return total_avg, total_mono

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
                    parsed_buffers = {
                        b: {
                            'activity': int(v.replace('*','').replace('<','')) if v.strip().replace('*','').replace('<','').isdigit() else 0, 
                            'star': '*' in v
                        } for b, v in buffers.items()
                    }
                    db[clean_name] = {'buffers': parsed_buffers, 'temp': row[9].strip()}
                except Exception:
                    pass
        return db
    except FileNotFoundError: 
        return None

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
            
    min_act = min(db[enz]['buffers'][best_buf]['activity'] for enz in combo)
    score = 100 - (50 if len(temps) > 1 else 0) - (100 - min_act) - (10 if star_risk else 0)
    info = f"{best_buf} ({min_act}% Act), Temp: {', '.join(temps)}"
    if star_risk: info += " ⚠️ Star Risk"
    return score, info

# ==========================================
# 3. SIDEBAR & INPUT
# ==========================================
with st.sidebar:
    st.header("1. Fragment Constraints")
    min_size = st.number_input("Min Size (bp)", value=30)
    max_size = st.number_input("Max Size (bp)", value=500)
    max_enz = st.slider("Max Enzymes in Cocktail", 1, 3, 2)
    
    st.divider()
    
    st.header("2. Topology / Linearization")
    is_linearized = st.checkbox("Plasmid will be pre-linearized")
    lin_enz_name = st.text_input("Linearizing Enzyme", "XhoI") if is_linearized else None

raw_seq = st.text_area("Paste Plasmid Sequence", height=150, placeholder="ATGC...")

# ==========================================
# 4. EXECUTION LOGIC (SEARCH)
# ==========================================
if st.button("Step 1: Analyze & Shred", type="primary"):
    if not raw_seq:
        st.error("Please paste a sequence first.")
        st.stop()

    neb_db = load_neb_data('neb_buffers.csv')
    if not neb_db:
        st.error("CSV file not found! Make sure 'neb_buffers.csv' is in the directory.")
        st.stop()

    clean_seq = "".join(raw_seq.split()).upper()
    seq_obj = Seq(clean_seq)
    seq_len = len(clean_seq)

    # --- Linearization logic ---
    if is_linearized:
        try:
            lin_enz = getattr(Restriction, lin_enz_name)
            cuts = lin_enz.search(seq_obj, linear=False)
            if len(cuts) != 1: 
                st.error(f"Linearizing enzyme must cut exactly once. Found {len(cuts)} sites.")
                st.stop()
            cut_pos = cuts[0]
            # Shift sequence to start at the cut site (Biopython sites are 1-indexed)
            clean_seq = clean_seq[cut_pos-1:] + clean_seq[:cut_pos-1]
            seq_obj = Seq(clean_seq)
            st.success(f"✅ Linearized at site {cut_pos} using {lin_enz_name}")
        except AttributeError:
            st.error(f"❌ '{lin_enz_name}' not found in the database. Check spelling.")
            st.stop()

    # --- Pre-search NEB enzymes ---
    with st.spinner("Finding combinations..."):
        neb_batch = Restriction.RestrictionBatch(first=[], suppliers=['N'])
        enzyme_cuts = {}
        for enz in neb_batch:
            enz_str = str(enz)
            if is_linearized and enz_str == lin_enz_name:
                continue # Don't shred with the linearizing enzyme
                
            cuts = enz.search(seq_obj, linear=not (not is_linearized))
            if cuts:
                enzyme_cuts[enz_str] = cuts
                
        active_enzymes = list(enzyme_cuts.keys())
        results = []
        
        # --- Find Combinations ---
        for r in range(1, max_enz + 1):
            for combo in itertools.combinations(active_enzymes, r):
                all_cuts = sorted(list(set().union(*(enzyme_cuts[e] for e in combo))))
                
                # Calculate fragment lengths
                if not is_linearized:
                    frags = [all_cuts[i+1]-all_cuts[i] for i in range(len(all_cuts)-1)] + [(seq_len-all_cuts[-1])+all_cuts[0]] if all_cuts else []
                else:
                    prev = 0
                    frags = []
                    for c in all_cuts:
                        frags.append(c - prev)
                        prev = c
                    frags.append(seq_len - prev)
                
                if frags and all(min_size <= f <= max_size for f in frags):
                    score, info = score_cocktail(combo, neb_db)
                    results.append({
                        "Score": score, 
                        "Cocktail": ", ".join(combo), 
                        "Info": info, 
                        "Fragments (bp)": str(sorted(frags)),
                        "Cuts": all_cuts
                    })

    if results:
        df_res = pd.DataFrame(results).sort_values("Score", ascending=False)
        st.success(f"Found {len(results)} valid cocktails!")
        
        # Save to session state so it doesn't disappear when interacting with the dropdown below
        st.session_state['results_df'] = df_res
        st.session_state['final_seq'] = seq_obj
        st.session_state['is_circ'] = not is_linearized
    else:
        st.warning("No matches found matching your size criteria.")

# ==========================================
# 5. DETAILED SELECTION & MASS ANALYSIS
# ==========================================
if 'results_df' in st.session_state and not st.session_state['results_df'].empty:
    st.divider()
    st.subheader("Step 2: Fragment Mass Analysis")
    st.markdown("Select a cocktail from your results below to view precise strand masses.")
    
    df_res = st.session_state['results_df']
    st.dataframe(df_res.drop(columns=["Cuts"]), use_container_width=True, hide_index=True)
    
    options = df_res['Cocktail'].tolist()
    selected_name = st.selectbox("Select a cocktail to analyze:", options)
    
    if selected_name:
        try:
            # Retrieve the selected cocktail's data
            row = df_res[df_res['Cocktail'] == selected_name].iloc[0]
            cuts = row['Cuts']
            
            # 1. BULLETPROOF TYPE HANDLING
            if isinstance(cuts, str):
                import ast
                cuts = ast.literal_eval(cuts)
                
            # Force the sequence back to a plain string
            full_seq_str = str(st.session_state['final_seq']).upper()
            
            # Slice the actual physical sequence fragments
            frag_seqs = []
            if not cuts:
                frag_seqs.append(full_seq_str)
            elif st.session_state['is_circ']:
                for i in range(len(cuts)-1): 
                    frag_seqs.append(full_seq_str[cuts[i]-1 : cuts[i+1]-1])
                # Handle the circular wrap-around junction
                frag_seqs.append(full_seq_str[cuts[-1]-1:] + full_seq_str[:cuts[0]-1])
            else:
                prev = 0
                for c in cuts:
                    frag_seqs.append(full_seq_str[prev : c-1])
                    prev = c - 1
                frag_seqs.append(full_seq_str[prev:])

            # Process masses for both strands
            mass_data = []
            for i, f_seq_str in enumerate(frag_seqs, 1):
                if len(f_seq_str) == 0: continue
                
                # 2. SAFE BIO-CONVERSION
                sense_seq = Seq(f_seq_str)
                anti_seq = sense_seq.reverse_complement()
                
                sense_avg, sense_mono = calculate_dna_masses(sense_seq)
                anti_avg, anti_mono = calculate_dna_masses(anti_seq)
                
                # Split out Sense and Antisense into full untruncated sequences
                mass_data.append({
                    "Fragment": i,
                    "Length (bp)": len(f_seq_str),
                    "Sense Sequence (5' -> 3')": str(sense_seq),
                    "Sense Mono (Da)": f"{sense_mono:.4f}",
                    "Sense Avg (Da)": f"{sense_avg:.2f}",
                    "Antisense Sequence (5' -> 3')": str(anti_seq),
                    "Antisense Mono (Da)": f"{anti_mono:.4f}",
                    "Antisense Avg (Da)": f"{anti_avg:.2f}"
                })
            
            st.write(f"**Mass Data for: {selected_name}**")
            
            # Use st.dataframe for automatic scrolling of long sequences
            st.dataframe(pd.DataFrame(mass_data), use_container_width=True, hide_index=True)
            
        except Exception as e:
            st.error(f"An error occurred while calculating masses: {e}")
