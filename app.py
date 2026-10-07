import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re
import os
from PIL import Image

# Sayfa Ayarları
st.set_page_config(page_title="SNP-Spesifik CAPS Tasarlayıcı", layout="wide", initial_sidebar_state="collapsed")

# Premium Dark Mode ve Siberpunk CSS (Hatalar Giderildi)
st.markdown("""
    <style>
    /* Arka planı logoya uygun koyu lacivert/siyah yap */
    .stApp {background-color: #030712;}
    
    /* Temel metinler beyaz (div ve span çıkarıldı, çakışma önlendi) */
    h1, h2, h3, p, label {color: #ffffff !important;}
    
    /* 🧬 Dizi Giriş Kutusu (Text Area) Özel Tasarımı - Biyoinformatik Terminal Havası */
    .stTextArea textarea {
        background-color: #020617 !important;
        color: #00ffcc !important; /* Dizi harfleri neon camgöbeği parlayacak */
        font-family: 'Courier New', Consolas, monospace !important; /* Kod fontu */
        font-size: 1.05rem !important;
        border: 2px solid #1e293b !important;
        border-radius: 8px;
    }
    .stTextArea textarea:focus {
        border: 2px solid #10b981 !important;
        box-shadow: 0 0 10px rgba(16, 185, 129, 0.5);
    }

    /* Havalı Buton Tasarımı */
    .stButton>button {
        background-color: transparent;
        color: #10b981 !important;
        font-weight: 800;
        font-size: 1.1rem;
        border-radius: 8px;
        border: 2px solid #10b981;
        box-shadow: 0 0 10px rgba(16, 185, 129, 0.3);
        transition: all 0.3s ease-in-out;
    }
    .stButton>button:hover {
        background-color: #10b981;
        color: #000000 !important;
        box-shadow: 0 0 20px rgba(16, 185, 129, 0.8);
        border: 2px solid #10b981;
    }
    
    /* Geliştirici Künyesi */
    .credit-text {
        font-size: 1.2rem;
        color: #00e5ff !important;
        font-weight: 800;
        border-left: 4px solid #00e5ff;
        padding-left: 15px;
        margin-top: 5px;
        letter-spacing: 0.5px;
    }
    </style>
""", unsafe_allow_html=True)

# Üst Başlık ve Logo Alanı
col_logo, col_title = st.columns([1, 4])

with col_logo:
    # Sunucuya geçtiğinde buraya direkt link verebilirsin. Örn: st.image("https://...", use_container_width=True)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    img_path = os.path.join(current_dir, "genom.jpg")
    try:
        img = Image.open(img_path)
        st.image(img, use_container_width=True)
    except FileNotFoundError:
        st.warning("[Logo Sunucudan Çekilecek]")

with col_title:
    st.markdown("<h1>🧬 SNP-Spesifik CAPS Markör Tasarlayıcı v1.0</h1>", unsafe_allow_html=True)
    st.markdown("### **Genom Okulu** | Biyoinformatik ve Moleküler Islah Platformu")
    st.markdown("<div class='credit-text'>Sistem Mimarı ve Baş Geliştirici: Doç. Dr. İbrahim Çelik</div>", unsafe_allow_html=True)

st.divider()

# Mini Restriksiyon Enzim Kütüphanesi
ENZYMES = {
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", 
    "TaqI": "TCGA", "Tsp509I": "AATT", "AluI": "AGCT", 
    "HaeIII": "GGCC", "RsaI": "GTAC", "MseI": "TTAA", 
    "SspI": "AATATT", "DraI": "TTTAAA", "EcoRV": "GATATC",
    "ApoI": "[AG]AATT[CT]", "HinfI": "GA[ATGC]TC"
}

def get_fragments(seq, enz_motif):
    cuts = [0]
    for match in re.finditer(enz_motif, seq):
        cut_pos = match.start() + len(match.group()) // 2
        cuts.append(cut_pos)
    cuts.append(len(seq))
    cuts.sort()
    frags = [cuts[i] - cuts[i-1] for i in range(1, len(cuts))]
    return sorted(frags, reverse=True)

default_seq = "ATCGTCAAGCTAGCCTAGCATCCGATTACGCATGCTTACGGATCGACTACTACGATCGAATCGAGCCTAGCATCCGATTACGCATGCTTACGGATCGACTACTACGATCGAATCGAGCCTAGCATCCGATTACGCATGCTTACGGATCGACTACTACGATCGAATCGAGCCTAGCATCCGATTACGCATGCTCGATCGTG[A/G]ATTCGCTAGCTATTACGGATCGACTACTACGATCGAATCGAGCCTAGCATCCGATTACGCATGCTTACGGATCGACTACTACGATCGAATCGAGCCTAG"

st.markdown("### 🧬 1. Hedef Genomik Diziyi Girin")
seq_input = st.text_area("Hedef Dizi (Parantezli SNP Formatı: Örn. [A/G])", default_seq, height=120)

if st.button("🚀 Algoritmayı Çalıştır ve Enzimleri Tara", use_container_width=True):
    clean_seq = seq_input.replace("\n", "").replace(" ", "").upper()
    match = re.search(r'(.*)\[(.)/(.)\](.*)', clean_seq)
    
    if not match:
        st.error("Sözdizimi Hatası: Lütfen [A/G] formatında bir SNP noktası içeren geçerli bir dizi girin.")
    else:
        prefix, a1, a2, suffix = match.groups()
        alel1_seq = prefix + a1 + suffix
        alel2_seq = prefix + a2 + suffix
        total_len = len(alel1_seq)
        
        caps_enzymes = {}
        
        for e_name, e_motif in ENZYMES.items():
            frags1 = get_fragments(alel1_seq, e_motif)
            frags2 = get_fragments(alel2_seq, e_motif)
            
            if frags1 != frags2:
                caps_enzymes[e_name] = {"Alel1": frags1, "Alel2": frags2}
        
        if not caps_enzymes:
            st.warning("Bu SNP noktası için kütüphanede ayırt edici bir polimorfik enzim profili bulunamadı.")
            st.session_state['caps_results'] = None
        else:
            st.session_state['caps_results'] = caps_enzymes
            st.session_state['total_len'] = total_len
            st.session_state['alleles'] = (a1, a2)
            st.success(f"Sistem Taraması Tamamlandı: {len(caps_enzymes)} adet uygun CAPS enzimi izole edildi.")

st.divider()

if 'caps_results' in st.session_state and st.session_state['caps_results']:
    caps_results = st.session_state['caps_results']
    total_len = st.session_state['total_len']
    a1, a2 = st.session_state['alleles']
    
    st.markdown("### 🖥️ 2. Dijital Jel Elektroforez Simülasyonu")
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("Aşağıdaki listeden analiz etmek istediğiniz enzimi seçin:")
        selected_enzyme = st.selectbox("Aday Enzimler:", list(caps_results.keys()))
        
        frags1 = caps_results[selected_enzyme]["Alel1"]
        frags2 = caps_results[selected_enzyme]["Alel2"]
        
        st.info(f"**{selected_enzyme} Restriksiyon Profili:**\n\n"
                f"🧫 **Alel 1 ({a1}):** {frags1} bp\n\n"
                f"🧫 **Alel 2 ({a2}):** {frags2} bp")
    
    with col2:
        bg_color = '#020617' 
        fig, ax = plt.subplots(figsize=(6, 8), facecolor=bg_color)
        ax.set_facecolor(bg_color)
        
        max_bp = ((total_len // 100) + 1) * 100 if total_len % 100 != 0 else total_len
        marker_bands = sorted(list(range(100, max_bp + 100, 100)), reverse=True)
        
        lanes = {
            "DNA\nMarker": marker_bands,
            f"Alel 1\n({a1}{a1})": frags1,
            f"Alel 2\n({a2}{a2})": frags2,
            f"Melez\n({a1}{a2})": sorted(list(set(frags1 + frags2)), reverse=True)
        }
        
        well_y = max_bp + 50
        for i in range(1, len(lanes) + 1):
            well = patches.Rectangle((i - 0.25, well_y), 0.5, max_bp * 0.05, fill=False, edgecolor='#334155', lw=2)
            ax.add_patch(well)
            ax.text(i, well_y + (max_bp * 0.06), "KUYU", color='#94a3b8', ha='center', fontsize=9, fontweight='bold')
            
        for i, (lane_name, bands) in enumerate(lanes.items()):
            x = i + 1
            for band in bands:
                alpha_core = 0.9 if lane_name != "DNA\nMarker" else 0.5
                neon_color = '#10b981' 
                ax.hlines(band, x - 0.25, x + 0.25, colors=neon_color, linewidth=18, alpha=alpha_core * 0.25)
                ax.hlines(band, x - 0.2, x + 0.2, colors='#ffffff', linewidth=4, alpha=alpha_core) 
                if lane_name == "DNA\nMarker":
                    ax.text(x - 0.35, band, f"{band}", color='#cbd5e1', va='center', ha='right', fontsize=10)

        ax.set_xticks(list(range(1, len(lanes) + 1)))
        ax.set_xticklabels(lanes.keys(), color='#ffffff', fontsize=11, fontweight='bold')
        
        ax.set_ylim(0, max_bp + (max_bp * 0.15))
        
        ax.tick_params(colors='#ffffff', axis='y', labelsize=10)
        for spine in ax.spines.values():
            spine.set_color('#1e293b')
        ax.set_ylabel("DNA Büyüklüğü (bp)", color='#ffffff', fontsize=12, fontweight='bold')
        
        st.pyplot(fig)