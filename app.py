from collections import Counter
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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
    n = len(seq)
    total = sum(aa_weights.get(aa, 0) for aa in seq)
    if n > 1:
        total -= (n - 1) * 18.015
    return total

def calc_extinction(seq):
    counts = Counter(seq)
    return sum(aa_extinction.get(aa, 0) * counts[aa] for aa in counts)

def net_charge(seq, pH):
    counts = Counter(seq)
    charge = 0
    charge += 1 / (1 + 10**(pH - pKa['Nterm']))
    charge -= 1 / (1 + 10**(pKa['Cterm'] - pH))
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

# ---- Streamlit UI ----
st.title("Protein Calculator")

# 세션 상태 초기화 (처음 실행 시 기본 서열 설정)
if "sequence" not in st.session_state:
    st.session_state.sequence = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"

# 입력창 (세션 상태 값 사용, key 지정)
seq_input = st.text_area(
    "Enter protein sequence (single-letter code):",
    st.session_state.get("sequence", ""),   # clear 후에는 ""로 표시됨
    height=300,
    key="seq_input"
)

# 버튼 두 개
col1, col2 = st.columns(2)
with col1:
    submit = st.button("Submit")
with col2:
    reset = st.button("Reset")

# Reset 버튼 → 세션 상태 전체 초기화 + 입력창 비우기
if reset:
    st.session_state.clear()   # 모든 세션 상태 초기화
    st.session_state.sequence = ""   # sequence도 빈 문자열로 설정
    st.session_state.seq_input = ""  # 입력창도 빈 문자열로 설정
    st.rerun()

# Submit 버튼 → 계산 실행
if submit and seq_input:
    seq = seq_input.strip().upper().replace(" ", "").replace("\n", "")
    st.session_state.sequence = seq  # 입력값을 세션 상태에 저장

    mw = calc_mw(seq)
    ext = calc_extinction(seq)
    pI = calc_pI(seq)

    st.write(f"**Molecular Weight:** {mw:.2f} Da")
    st.write(f"**Extinction Coefficient (280nm):** {ext} M^-1 cm^-1")
    st.write(f"**Isoelectric Point (pI):** {pI:.2f}")

    calibrated_con = ext / mw if mw > 0 else 0
    st.write(f"**Calibrated Con (A280nm):** {calibrated_con:.6f}")

    st.caption("Calibrated Con 값은 실제 A280nm 흡광도 측정 결과에 대해 "
               "과대 또는 과소 추정될 수 있으므로, 반드시 실측 결과 값을 "
               "Calibrated Con으로 나누어 해석해야 합니다.")
    st.caption("The Calibrated Con value may be over-estimated or under-estimated "
               "relative to actual A280nm absorbance measurements, "
               "and should always be interpreted by dividing the measured result "
               "by the Calibrated Con.")

    ph_values = np.arange(4.0, 10.5, 0.5)
    results = []
    for ph in ph_values:
        charge = net_charge(seq, ph)
        results.append({"pH": ph, "Net Charge": charge})
    df = pd.DataFrame(results)

    st.subheader("Net Charge Table (pH 4.0 ~ 10.0, step 0.5)")
    st.table(df)

    fig, ax = plt.subplots()
    ax.plot(df["pH"], df["Net Charge"], marker="o", linestyle="-", color="blue")
    ax.axhline(0, color='gray', linestyle='--')
    ax.set_xlabel("pH")
    ax.set_ylabel("Net Charge")
    ax.set_title("Net Charge vs pH (4.0 ~ 10.0)")
    st.pyplot(fig)

    pHs = np.linspace(0,14,200)
    charges = [net_charge(seq,pH) for pH in pHs]
    fig2, ax2 = plt.subplots()
    ax2.plot(pHs, charges)
    ax2.axhline(0, color='gray', linestyle='--')
    ax2.set_xlabel("pH")
    ax2.set_ylabel("Net Charge")
    ax2.set_title("Charge vs pH (0 ~ 14)")
    st.pyplot(fig2)

