# gene_melody_web.py
import streamlit as st
import pandas as pd
from collections import Counter
import mido
from mido import MidiFile, MidiTrack, Message
import base64
import io

# -----------------------
# Constants
# -----------------------
BASE_NOTE_MAP = {'A': 60, 'T': 62, 'C': 64, 'G': 67}
INSTRUMENTS = {
    0: "Acoustic Grand Piano", 40: "Violin", 41: "Viola", 56: "Trumpet",
    60: "French Horn", 73: "Flute", 81: "Lead Synth", 89: "New Age Pad"
}
DNA_MEANINGS = {

    # TATA box variants
    "TATAAA": "TATA box (variant: TATAAA) — promoter element near TSS.",
    "TATATA": "TATA box (variant: TATATA) — promoter element near TSS.",
    "TATATT": "TATA box (variant: TATATT) — promoter element near TSS.",
    "TATAAT": "TATA box (variant: TATAAT) — promoter element near TSS.",
    "TATGAA": "TATA-like box (variant: TATGAA) — promoter element near TSS.",

    # CAAT box variants (often ~50–200 bp upstream of TSS)
    "CAAT":   "CAAT box — TF-binding site (often ~50–200 bp upstream of TSS).",
    "CCAAT":  "CAAT box (CCAAT) — TF-binding site (~50–200 bp upstream of TSS).",
    "CAATT":  "CAAT box (CAATT) — TF-binding site (~50–200 bp upstream of TSS).",
    "CCATT":  "CAAT box (CCATT) — TF-binding site (~50–200 bp upstream of TSS).",

    # GC box variants (Sp1-like binding, can be multiple copies)
    "GGGCGG": "GC box (GGGCGG) — GC-rich TF-binding site.",
    "GGGAGG": "GC box (GGGAGG) — GC-rich TF-binding site.",
    "GGGCCG": "GC box (GGGCCG) — GC-rich TF-binding site.",
    "GCGGGG": "GC box (GCGGGG) — GC-rich TF-binding site.",

    # Palindromic restriction sites
    "GAATTC": "EcoRI restriction site (GAATTC).",
    "GGATCC": "BamHI restriction site (GGATCC).",
    "AAGCTT": "HindIII restriction site (AAGCTT).",
    "CTGCAG": "PstI restriction site (CTGCAG).",
    "CCCGGG": "SmaI restriction site (CCCGGG).",
    "GCGGCCGC": "NotI restriction site (GCGGCCGC).",
    "TCTAGA": "XbaI restriction site (TCTAGA).",
    "GCTAGC": "NheI restriction site (GCTAGC).",
    "GGTACC": "KpnI restriction site (GGTACC).",
    "GAGCTC": "SacI restriction site (GAGCTC).",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}
AVG_NT_DA = 330.0

# -----------------------
# Helper functions
# -----------------------
def clean_dna(raw: str) -> str:
    s = raw.upper()
    return ''.join([c for c in s if c in BASE_NOTE_MAP])

def base_counts_and_stats(dna: str):
    counts = Counter(dna)
    total = len(dna) if dna else 1
    gc = counts.get('G', 0) + counts.get('C', 0)
    at = counts.get('A', 0) + counts.get('T', 0)
    gc_percent = (gc / total) * 100
    at_percent = (at / total) * 100
    at_gc_ratio = (at / gc) if gc > 0 else float('inf')
    mw_da = len(dna) * AVG_NT_DA
    mw_kda = mw_da / 1000.0
    return {
        "counts": counts,
        "total": len(dna),
        "gc_percent": gc_percent,
        "at_percent": at_percent,
        "at_gc_ratio": at_gc_ratio,
        "mw_da": mw_da,
        "mw_kda": mw_kda
    }

def find_motifs(dna: str, motifs: dict):
    found = []
    for motif, meaning in motifs.items():
        start = 0
        while True:
            idx = dna.find(motif, start)
            if idx == -1:
                break
            found.append({"motif": motif, "start": idx + 1, "meaning": meaning})
            start = idx + 1
    return sorted(found, key=lambda x: x["start"])

def find_orfs_all_frames(dna: str):
    orfs = []
    n = len(dna)
    for frame in range(3):
        i = frame
        while i + 3 <= n:
            codon = dna[i:i+3]
            if codon == "ATG":
                j = i
                while j + 3 <= n:
                    cod = dna[j:j+3]
                    if cod in STOP_CODONS:
                        length_nt = (j + 3) - i
                        length_aa = length_nt // 3
                        orfs.append({
                            "frame": frame,
                            "start": i + 1,
                            "end": j + 3,
                            "length_nt": length_nt,
                            "length_aa": length_aa
                        })
                        break
                    j += 3
                i += 3
            else:
                i += 3
    return orfs

def summarize_orfs(orfs):
    if not orfs:
        return "No ORFs found."
    longest = max(orfs, key=lambda o: o["length_nt"])
    summary = [f"Total ORFs: {len(orfs)}",
               f"Longest ORF: frame {longest['frame']} start {longest['start']} end {longest['end']} "
               f"({longest['length_nt']} nt / {longest['length_aa']} aa)"]
    return "\n".join(summary)

def dna_to_midi_file(dna_sequence, bpm=120, program=81):
    midi = MidiFile()
    track = MidiTrack()
    midi.tracks.append(track)
    track.append(Message('program_change', program=program, time=0))
    track.append(mido.MetaMessage('set_tempo', tempo=mido.bpm2tempo(bpm)))
    for base in dna_sequence:
        note = BASE_NOTE_MAP.get(base)
        if note:
            track.append(Message('note_on', note=note, velocity=100, time=0))
            track.append(Message('note_off', note=note, velocity=100, time=200))
    buf = io.BytesIO()
    midi.save(file=buf)
    return buf.getvalue()

# -----------------------
# Streamlit UI
# -----------------------
st.set_page_config(page_title="GeneMelody Web", page_icon="🎵", layout="centered")
st.title("🎵 GeneMelody - DNA to Music")
st.write("Convert DNA sequences into musical compositions with biological annotations.")

# DNA input
dna_input = st.text_area("Enter DNA Sequence (A, T, C, G):", height=150)
uploaded_file = st.file_uploader("Or upload a DNA/FASTA file", type=["txt", "fasta"])
if uploaded_file:
    dna_input = ''.join(line.strip() for line in uploaded_file if not line.startswith(">"))
    dna_input = clean_dna(dna_input)

# Tempo and instrument
bpm = st.number_input("Tempo (BPM):", min_value=1, value=120)
instrument_list = [f"{v} ({k})" for k, v in INSTRUMENTS.items()]
instrument_choice = st.selectbox("Instrument:", instrument_list)
program = int(instrument_choice.split("(")[-1].strip(")"))

# Process
if st.button("Generate Music & Analysis"):
    if not dna_input:
        st.error("Please enter or upload a DNA sequence.")
    else:
        dna = clean_dna(dna_input)
        stats = base_counts_and_stats(dna)
        motifs_found = find_motifs(dna, DNA_MEANINGS)
        orfs = find_orfs_all_frames(dna)

        # Display Stats
        st.subheader("Base Composition & Physical Properties")
        st.write({
            "Length (bp)": stats["total"],
            "GC%": f"{stats['gc_percent']:.2f}%",
            "AT%": f"{stats['at_percent']:.2f}%",
            "AT/GC ratio": "∞" if stats['at_gc_ratio'] == float('inf') else f"{stats['at_gc_ratio']:.3f}",
            "Molecular weight": f"{stats['mw_da']:.0f} Da ({stats['mw_kda']:.2f} kDa)"
        })
        st.bar_chart(pd.DataFrame.from_dict(stats["counts"], orient='index', columns=['Count']))

        # Motifs
        st.subheader("Motifs Found")
        if motifs_found:
            st.dataframe(pd.DataFrame(motifs_found))
        else:
            st.write("No motifs found.")

        # ORF Summary
        st.subheader("ORF Analysis")
        st.text(summarize_orfs(orfs))

        # Generate MIDI file
        midi_data = dna_to_midi_file(dna, bpm, program)
        st.download_button("Download MIDI File", midi_data, file_name="dna_music.mid", mime="audio/midi")

        # Generate analysis text
        analysis_text = f"""
DNA Sequence Analysis
Length: {stats['total']} bp
GC%: {stats['gc_percent']:.2f}%
AT%: {stats['at_percent']:.2f}%
AT/GC ratio: {"∞" if stats['at_gc_ratio'] == float('inf') else f"{stats['at_gc_ratio']:.3f}"}
Molecular weight: {stats['mw_da']:.0f} Da ({stats['mw_kda']:.2f} kDa)

Motifs Found:
{pd.DataFrame(motifs_found).to_string(index=False) if motifs_found else "None"}

ORF Summary:
{summarize_orfs(orfs)}
"""
        st.download_button("Download Analysis TXT", analysis_text, file_name="analysis.txt")

