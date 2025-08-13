import streamlit as st
from mido import MidiFile, MidiTrack, Message, MetaMessage
import os
import tempfile

# Constants
BASE_NOTE_MAP = {'A': 60, 'T': 62, 'C': 64, 'G': 67}
INSTRUMENTS = {
    0: "Acoustic Grand Piano", 40: "Violin", 41: "Viola",
    56: "Trumpet", 60: "French Horn", 73: "Flute",
    81: "Lead Synth", 89: "New Age Pad"
}

# UI setup with favicon
st.set_page_config(
    page_title="GeneMelody",
    page_icon="favicon.ico",  # Local icon file
    layout="centered"
)


st.title("🎼 GeneMelody – DNA Music Generator")

# File uploader or manual input
uploaded_file = st.file_uploader("Upload DNA/FASTA file (.txt, .fasta)", type=["txt", "fasta"])
if uploaded_file:
    file_contents = uploaded_file.read().decode("utf-8")
    dna_sequence = ''.join([line.strip() for line in file_contents.splitlines() if not line.startswith(">")])
else:
    dna_sequence = ""

text_input = st.text_area("Or enter DNA Sequence (A, T, C, G only):", value=dna_sequence, height=150)
dna_sequence = text_input.upper().strip()

# Validation
valid_sequence = ''.join([c for c in dna_sequence if c in BASE_NOTE_MAP])
invalid_chars = set(dna_sequence) - set(BASE_NOTE_MAP)
if invalid_chars:
    st.warning(f"Ignored invalid characters: {', '.join(invalid_chars)}")

# Controls
tempo = st.number_input("Tempo (BPM):", min_value=40, max_value=300, value=120)
instrument_label = st.selectbox("Choose Instrument:", [f"{name} ({num})" for num, name in INSTRUMENTS.items()])
instrument_program = int(instrument_label.split('(')[-1].strip(')'))

# Generate button
generate_btn = st.button("🎼 Generate MIDI")

# MIDI generation
def save_midi(sequence, bpm, program, filename):
    midi = MidiFile()
    track = MidiTrack()
    midi.tracks.append(track)
    track.append(Message('program_change', program=program, time=0))
    track.append(MetaMessage('set_tempo', tempo=int(60000000 / bpm)))
    for base in sequence:
        note = BASE_NOTE_MAP.get(base)
        if note:
            track.append(Message('note_on', note=note, velocity=100, time=0))
            track.append(Message('note_off', note=note, velocity=100, time=200))
    midi.save(filename)

# Main action
if generate_btn and valid_sequence:
    with tempfile.TemporaryDirectory() as tmpdir:
        midi_path = os.path.join(tmpdir, "gene.mid")
        save_midi(valid_sequence, tempo, instrument_program, midi_path)

        # MIDI download only
        with open(midi_path, "rb") as f:
            st.success("MIDI file generated!")
            st.download_button("⬇️ Download MIDI", f, file_name="gene_music.mid", mime="audio/midi")
