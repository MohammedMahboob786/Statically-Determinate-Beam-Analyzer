# Statically Determinate Beam Analyzer

This Streamlit application provides an interactive tool for analyzing different types of statically determinate beams. Users can select from three types of beams - Cantilever Beam, Simply Supported Beam, and Overhang Beam - to view Shear Force Diagrams, Bending Moment Diagrams, and Deflection plots.

## Features

- **Cantilever Beam Analysis**
- **Simply Supported Beam Analysis**
- **Overhang Beam Analysis**
- Interactive plots for Shear Force, Bending Moment, and Deflection
- Easy navigation between different beam types

## File Structure
1. app.py: Main application script containing the Streamlit app logic.
2. cantiliver.png, simply_supported.png, overhang.png: Images used in the application for beam selection.
3. requirements.txt: List of required Python packages.

## Usage
1. Select a Beam Type: On the home page, click on the image or button corresponding to the beam type you want to analyze.
2. View Analysis: The application will display the Shear Force Diagram, Bending Moment Diagram, and Deflection plot for the selected beam type.
3. Navigate: Use the "Return to Beam Selection" button to go back to the home page and choose a different beam type.

## Installation

To run this application locally, follow these steps:

1. **Clone the repository:**
   ```sh
   git clone https://github.com/MohammedMahboob786/Statically-Determinate-Beam-Analyzer.git
   cd Statically-Determinate-Beam-Analyzer

2. **Create a virtual environment and activate it:**
    python -m venv venv
    ### On Windows
    venv\Scripts\activate
    ### On macOS/Linux
    source venv/bin/activate

3. **Install the required packages:**
    pip install -r requirements.txt

4. **Run the Streamlit app:**
    python -m streamlit run app.py





