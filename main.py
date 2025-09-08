import streamlit as st
import numpy as np
import plotly.graph_objects as go
from PIL import Image
from catiliver_beam import CantileverBeamApp
from simply_supported import SimplySupportedBeamApp
from overhang_beam import OverhangBeamApp

# Set page config
st.set_page_config(page_title="Beam Analysis Tool", page_icon="🔧")

# Main function to handle page navigation
def main():
    # st.set_page_config(page_title="Beam Analyzer", layout="centered")
    st.title("Statically Determinate Beam Analyzer")
    
    st.markdown("""
    This interactive tool allows you to generate **Shear Force Diagrams**, **Bending Moment Diagrams**, 
    and **Deflection Plots** for:
    -  Cantilever Beams  
    -  Simply Supported Beams  
    -  Overhang Beams  
    """)

    # Initialize session state
    if "page" not in st.session_state:
        st.session_state["page"] = "home"

    if st.session_state["page"] == "home":
        display_home()
    else:
        display_analysis()


def display_home():
    st.markdown("### 🛠️ Choose Beam Type")
    st.markdown("Click on any beam type below to begin your analysis:")

    col1, col2, col3 = st.columns(3)

    # Load images
    cantilever_img = Image.open("images/cantiliver.png")
    simply_supported_img = Image.open("images/simply_supported.png")
    overhang_img = Image.open("images/overhang.png")

    # Cantilever Beam
    with col1:
        st.image(cantilever_img, use_column_width=True)
        st.markdown("<div style='text-align:center;'> <b>Cantilever Beam</b></div>", unsafe_allow_html=True)
        st.button("Analyze Cantilever", use_container_width=True, on_click=set_page, args=("cantilever",))

    # Simply Supported Beam
    with col2:
        st.image(simply_supported_img, use_column_width=True)
        st.markdown("<div style='text-align:center;'> <b>Simply Supported</b></div>", unsafe_allow_html=True)
        st.button("Analyze Simply Supported", use_container_width=True, on_click=set_page, args=("simply_supported",))

    # Overhang Beam
    with col3:
        st.image(overhang_img, use_column_width=True)
        st.markdown("<div style='text-align:center;'> <b>Overhang Beam</b></div>", unsafe_allow_html=True)
        st.button("Analyze Overhang", use_container_width=True, on_click=set_page, args=("overhang",))

        
    
def set_page(page):
    st.session_state["page"] = page

def display_analysis():
    selected_page = st.session_state["page"]
    if selected_page == "cantilever":
        # cantilever_beam()
        CantileverBeamApp().render()
    elif selected_page == "simply_supported":
        # simply_supported_beam()
        SimplySupportedBeamApp().render()
    elif selected_page == "overhang":
        # overhang_beam()
        OverhangBeamApp().render()

    if st.button("Return to Beam Selection", on_click=set_page, args=("home",)):
        st.session_state["page"] = "home"


if __name__ == '__main__':
    main()
