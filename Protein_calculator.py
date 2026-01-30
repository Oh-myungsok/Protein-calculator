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
    """분자량 계산: 아미노산 개별 분자량 합 - (펩타이드 결합 수 × 물 분자 질량)"""
    n = len(seq)
    total = sum(aa_weights.get(aa, 0) for aa in seq)
    if n > 1:
        total -= (n - 1) * 18.015  # 펩타이드 결합 형성 시 빠지는 물(H2O)
    return total

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

seq = st.text_area("Enter protein sequence (single-letter code):", 
                   "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE", height=300)

if seq:
    mw = calc_mw(seq)
    ext = calc_extinction(seq)
    pI = calc_pI(seq)

    st.write(f"**Molecular Weight:** {mw:.2f} Da")
    st.write(f"**Extinction Coefficient (280nm):** {ext} M^-1 cm^-1")
    st.write(f"**Isoelectric Point (pI):** {pI:.2f}")

    # 추가: Calibrated Con (A280nm) = Extinction Coefficient / Molecular Weight 
    calibrated_con = ext / mw if mw > 0 else 0
    st.write(f"**Calibrated Con (A280nm):** {calibrated_con:.6f}")

    # 설명 문장 caption으로 추가 (한국어 + 영어) 
    st.caption("Calibrated Con 값은 실제 A280nm 흡광도 측정 결과에 대해 "
               "과대 또는 과소 추정될 수 있으므로, 반드시 실측 결과 값을 "
               "Calibrated Con으로 나누어 해석해야 합니다.")
    st.caption("The Calibrated Con value may be over-estimated or under-estimated "
               "relative to actual A280nm absorbance measurements, "
               "and should always be interpreted by dividing the measured result "
               "by the Calibrated Con.")

    # pH 4.0 ~ 10.0, 0.5 간격으로 Net Charge 테이블 출력
    ph_values = np.arange(4.0, 10.5, 0.5)
    results = []
    for ph in ph_values:
        charge = net_charge(seq, ph)
        results.append({"pH": ph, "Net Charge": charge})
    df = pd.DataFrame(results)

    st.subheader("Net Charge Table (pH 4.0 ~ 10.0, step 0.5)")
    st.table(df)

    # 테이블과 같은 범위의 그래프 추가
    fig, ax = plt.subplots()
    ax.plot(df["pH"], df["Net Charge"], marker="o", linestyle="-", color="blue")
    ax.axhline(0, color='gray', linestyle='--')
    ax.set_xlabel("pH")
    ax.set_ylabel("Net Charge")
    ax.set_title("Net Charge vs pH (4.0 ~ 10.0)")
    st.pyplot(fig)

    # 기존 Charge vs pH 전체 범위 그래프도 유지
    pHs = np.linspace(0,14,200)
    charges = [net_charge(seq,pH) for pH in pHs]
    fig2, ax2 = plt.subplots()
    ax2.plot(pHs, charges)
    ax2.axhline(0, color='gray', linestyle='--')
    ax2.set_xlabel("pH")
    ax2.set_ylabel("Net Charge")
    ax2.set_title("Charge vs pH (0 ~ 14)")
    st.pyplot(fig2)
