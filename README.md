# Protein Calculator 🧬

단백질 서열(1-letter amino acid code)을 입력하면 분자량(Molecular Weight), 
흡광계수(Extinction Coefficient), 등전점(pI), Net Charge 등을 계산해주는 
Streamlit 기반 웹 애플리케이션입니다.

## Features
- 단백질 서열 입력 및 분석
- 분자량(Molecular Weight) 계산
- 280nm 흡광계수(Extinction Coefficient) 계산
- 등전점(Isoelectric Point, pI) 추정
- pH에 따른 Net Charge 테이블 및 그래프 시각화
- Calibrated Con 값 계산 및 안내

## Demo
앱은 Streamlit Cloud에 배포되어 있으며, 아래 링크에서 바로 실행할 수 있습니다:

👉 [Protein Calculator App](https://protein-calculator.streamlit.app)

또는 GitHub Pages 버전(간단 계산기, JS 변환 버전)은 여기서 확인할 수 있습니다:

👉 [Protein Calculator (GitHub Pages)](https://oh-myungsok.github.io/Protein-calculator/)

## Installation (Local 실행)
로컬 환경에서 실행하려면 다음 단계를 따르세요:

```bash
# 저장소 클론
git clone https://github.com/Oh-myungsok/Protein-calculator.git
cd Protein-calculator

# 필요한 라이브러리 설치
pip install -r requirements.txt

# 앱 실행
streamlit run app.py



_________________________________________________________________________________________________________________

Protein Calculator 🧬 (English Version)
A Streamlit-based web application that calculates protein properties such as
Molecular Weight, Extinction Coefficient, Isoelectric Point (pI), and Net Charge
from a given amino acid sequence (single-letter code).

Features
Input and analyze protein sequences

Calculate Molecular Weight

Calculate Extinction Coefficient at 280nm

Estimate Isoelectric Point (pI)

Display Net Charge table and graphs across pH values

Provide Calibrated Con value with interpretation notes

Demo
The app is deployed on Streamlit Cloud and can be accessed here:

👉 Protein Calculator (Streamlit Cloud)

A simplified GitHub Pages version (JavaScript-based calculator) is also available:

👉 Protein Calculator (GitHub Pages)

Installation (Local Run)
To run locally:

bash
# Clone repository
git clone https://github.com/Oh-myungsok/Protein-calculator.git
cd Protein-calculator

# Install dependencies
pip install -r requirements.txt

# Run app
streamlit run app.py
Requirements
streamlit

numpy

matplotlib

pandas

(optional) biopython

Repository Structure
코드
Protein-calculator/
 ├── app.py              # Main Streamlit app code
 ├── requirements.txt    # Dependencies
 ├── README.md           # Project description
 └── (other files)
Author
Created by Oh Myungsok

GitHub: Oh-myungsok (github.com in Bing)
