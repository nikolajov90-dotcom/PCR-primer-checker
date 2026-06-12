import streamlit as st
import pandas as pd
import base64

from core import (
    analyze_primers,
    get_amplicon_sequence
)


# -------------------------------
# Background image
# -------------------------------

def set_bg(img_file):

    with open(img_file, "rb") as f:
        data = f.read()

    encoded = base64.b64encode(
        data
    ).decode()

    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image:
            url("data:image/jpg;base64,{encoded}");
            background-size: cover;
            background-attachment: fixed;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )


set_bg("background.jpg")


# -------------------------------
# Styling
# -------------------------------

st.markdown("""
<style>

[data-testid="stMainBlockContainer"] {
    background-color:
    rgba(255,255,255,0.85);
    padding: 2rem;
    border-radius: 20px;
    box-shadow:
    0px 0px 30px rgba(0,0,0,0.3);
}

[data-testid="stMainBlockContainer"] * {
    color: black !important;
}

@media (prefers-color-scheme: dark) {

  [data-testid="stMainBlockContainer"] {
    background-color:
    rgba(30,30,30,0.85);
  }

  [data-testid="stMainBlockContainer"] * {
    color: white !important;
  }

}

</style>
""", unsafe_allow_html=True)


# -------------------------------
# Header
# -------------------------------

st.markdown("""
# PCR Primer Specificity Checker

### Laboratorija za Molekularnu biologiju

Departman za biologiju i ekologiju

Prirodno-matematički fakultet

Univerzitet u Nišu
""")


# -------------------------------
# User input
# -------------------------------

forward = st.text_input(
    "Forward primer"
)

reverse = st.text_input(
    "Reverse primer"
)

uploaded_file = st.file_uploader(
    "Upload referentnog genoma (FASTA)",
    type=[
        "fasta",
        "fa",
        "fna"
    ]
)

col1, col2 = st.columns(2)

with col1:

    max_mismatches = st.slider(
        "Maksimum \"mismatch\" pozicija",
        min_value=0,
        max_value=5,
        value=0
    )

with col2:

    max_amplicon_length = st.number_input(
        "Maksimalna dužina amplikona (bp)",
        min_value=50,
        value=1000,
        step=50
    )


# -------------------------------
# Analysis
# -------------------------------

if st.button("Analyze"):

    if not forward:

        st.error(
            "Molim unesite \"forward\" prajmer."
        )

        st.stop()

    if not reverse:

        st.error(
            "Molim unesite \"reverse\" prajmer."
        )

        st.stop()

    if uploaded_file is None:

        st.error(
            "Molim upload FASTA fajl-a."
        )

        st.stop()

    with st.spinner(
            "Pretraga mesta vezivanja \"prajmera\"..."):

        try:

            genomes, amplicons = analyze_primers(
                uploaded_file,
                forward,
                reverse,
                max_mismatches=max_mismatches,
                max_length=max_amplicon_length
            )       

        except Exception as e:

            st.error(
                f"Error: {e}"
            )

            st.stop()


    # ---------------------------
    # Summary
    # ---------------------------

    total_amplicons = sum(
        len(products)
        for products
        in amplicons.values()
    )

    st.success(
        f"Analiza gotova. "
        f"Pronađeno {total_amplicons} "
        f"potencijalnih amplikona."
    )

    # ---------------------------
    # Results table
    # ---------------------------

    rows = []

    for chromosome, products in amplicons.items():

        for (
            start,
            end,
            size,
            strand,
            f_mm,
            r_mm,
            total_mm
        ) in products:

            strand_label = (
                "Plus (+)"
                if strand == "+"
                else "Minus (-)"
            )

            rows.append({

                "Hromozom":
                    chromosome,

                "Forward pozicija":
                    start,

                "Reverse pozicija":
                    end,

                "Matrični lanac":
                    strand_label,

                "Duzina amplikona (bp)":
                    size,

                "Broj forward \"mismatch\" pozicija":
                    f_mm,

                "Broj reverse \"mismatch\" pozicija":
                    r_mm,

                "Ukupan broj \"mismatch\" pozicija":
                    total_mm
            })

    if rows:

        df = pd.DataFrame(rows)

        st.subheader(
            "Amplikoni"
        )

        st.dataframe(
            df,
            use_container_width=True
        )

    else:

        st.warning(
            "Bez pronađenih amplikona."
        )

    # ---------------------------
    # Detailed view
    # ---------------------------

    st.subheader(
        "Detalji produkata"
    )

    for chromosome, products in amplicons.items():

        with st.expander(
                f"{chromosome} "
                f"({len(products)} amplicons)"):

            for (
                start,
                end,
                size,
                strand,
                f_mm,
                r_mm,
                total_mm
            ) in products:

                strand_label = (
                "Forward → Reverse"
                if strand == "+"
                else "Reverse → Forward"
            )


                st.markdown(
                    f"""
                    **Forward:** {start}

                    **Reverse:** {end}

                    **Veličina:** {size} bp

                    **Matrični lanac:** {strand_label}

                    **Broj Forward mismatch pozicija:** {f_mm}

                    **Broj Reverse mismatch pozicija:** {r_mm}

                    **Ukupno mismatch pozicija:** {total_mm}
                    """
                )

                seq = (
                    get_amplicon_sequence(
                        genomes,
                        chromosome,
                        start,
                        end,
                        len(reverse)
                    )
                )

                st.code(
                    seq[:500]
                )

                st.divider()
