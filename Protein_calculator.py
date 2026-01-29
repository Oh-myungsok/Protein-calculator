from collections import Counter
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter


# 아미노산별 분자량 (Da)
aa_weights = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10,
    'C': 121.16, 'E': 147.13, 'Q': 146.15, 'G': 75.07,
    'H': 155.16, 'I': 131.18, 'L': 131.18, 'K': 146.19,
    'M': 149.21, 'F': 165.19, 'P': 115.13, 'S': 105.09,
    'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}

# pKa 값 (단순화된 버전)
pKa = {
    'Cterm': 3.1, 'Nterm': 8.0,
    'D': 3.9, 'E': 4.1, 'C': 8.3, 'Y': 10.1,
    'H': 6.0, 'K': 10.5, 'R': 12.5
}

# 280nm 흡광계수
aa_extinction = {'W': 5500, 'Y': 1490, 'C': 125}

def calc_mw(seq):
    return sum(aa_weights[aa] for aa in seq)

def calc_extinction(seq):
    counts = Counter(seq)
    return sum(aa_extinction.get(aa, 0) * counts[aa] for aa in counts)

def net_charge(seq, pH):
    counts = Counter(seq)
    charge = 0
    # N-term
    charge += 1 / (1 + 10**(pH - pKa['Nterm']))
    # C-term
    charge -= 1 / (1 + 10**(pKa['Cterm'] - pH))
    # Side chains
    for aa in ['D','E','C','Y']:
        charge -= counts[aa] / (1 + 10**(pKa[aa] - pH))
    for aa in ['H','K','R']:
        charge += counts[aa] / (1 + 10**(pH - pKa[aa]))
    return charge

def calc_pI(seq):
    pHs = np.linspace(0,14,1401)
    charges = [net_charge(seq,pH) for pH in pHs]
    idx = np.argmin(np.abs(charges))
    return pHs[idx]

# Streamlit UI
st.title("Protein Calculator")

seq = st.text_input("Enter protein sequence (single-letter code):", "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE")

if seq:
    mw = calc_mw(seq)
    ext = calc_extinction(seq)
    pI = calc_pI(seq)

    st.write(f"**Molecular Weight:** {mw:.2f} Da")
    st.write(f"**Extinction Coefficient (280nm):** {ext} M^-1 cm^-1")
    st.write(f"**Isoelectric Point (pI):** {pI:.2f}")

    # Charge vs pH plot
    pHs = np.linspace(0,14,200)
    charges = [net_charge(seq,pH) for pH in pHs]
    fig, ax = plt.subplots()
    ax.plot(pHs, charges)
    ax.axhline(0, color='gray', linestyle='--')
    ax.set_xlabel("pH")
    ax.set_ylabel("Net Charge")
    ax.set_title("Charge vs pH")
    st.pyplot(fig)
