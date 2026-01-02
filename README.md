# Aerospace Structure Instabilities Analysis Tool

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://YOUR-APP-URL.streamlit.app)

A Python-based verification tool for aerospace structural engineering, specifically focused on instability cases found in *StructAeroInstabilities.pdf*. This application uses **Streamlit** to provide an interactive UI with detailed step-by-step logging, LaTeX formula rendering, and stability visualization plots.

## 🚀 Features

The tool covers four distinct analytical cases from the course notes:

### 1. Flexural-Torsional Buckling (Pages 38-51)
* **Analysis:** Open thin-walled columns (C-channels).
* **Calculations:** Centroid, Shear Center, Warping Constant ($C_\Gamma$), Torsion Constant ($J$), and coupled critical loads.
* **Visuals:** Stability curve ($P_{cr}$ vs Length) comparing coupled vs. uncoupled modes.

### 2. Buckling of Thin Plates (Pages 60-74)
* **Analysis:** Critical stress and buckling coefficients ($k$) for plates under compression.
* **Calculations:** Supports SSSS (Simply Supported) and SSSF (One edge free) boundary conditions.
* **Visuals:** 3D Buckling Mode Shapes ($m, n$) and Buckling Coefficient charts ($k$ vs $a/b$).

### 3. Transversely Loaded Columns (Annex 1, Pages 76-82)
* **Analysis:** Beam-columns with axial compression and transverse loads.
* **Calculations:** Exact solutions for differential equations, integration constants, and magnification factors.
* **Visuals:** Deflection vs. Axial Load plot demonstrating the singularity at Euler buckling load.

### 4. Stiffener/Web Constructions (Annex 2, Pages 89-93)
* **Analysis:** Wagner beam theory for shear webs reinforced with stiffeners.
* **Calculations:** Diagonal tension angle ($\alpha$), flange stresses (direct + bending), and stiffener buckling checks.
* **Visuals:** Stiffener Stability Map showing safe vs. failure regions based on spacing.

## 🛠️ Installation & Usage

To run this application locally on your machine:

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/YOUR_USERNAME/structure-instability-analysis.git](https://github.com/YOUR_USERNAME/structure-instability-analysis.git)
    cd structure-instability-analysis
    ```

2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Run the App:**
    ```bash
    streamlit run app.py
    ```

## 📂 Project Structure

* `app.py`: The main application script containing all calculation logic and UI definitions.
* `requirements.txt`: List of Python dependencies (`streamlit`, `numpy`, `matplotlib`).
* `README.md`: Documentation.

## 📦 Dependencies

* [Streamlit](https://streamlit.io/) - For the web interface.
* [NumPy](https://numpy.org/) - For numerical calculations.
* [Matplotlib](https://matplotlib.org/) - For generating 2D and 3D plots.

## 📝 Reference

This tool is based on the *Aircraft Structures - Instabilities* lecture notes (University of Liège), specifically:
* **Flexural-Torsional example:** Page 38
* **Plate Buckling charts:** Page 61
* **Annex 1 (Beam-Column):** Page 76
* **Annex 2 (Wagner Web):** Page 90