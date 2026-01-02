import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# ==========================================
# UI HELPER FUNCTIONS
# ==========================================
def log_step(step_name, formula_latex, result_value, unit, page_ref, explanation=""):
    """
    Renders a consistent log entry for a calculation step with LaTeX support.
    """
    st.markdown(f"#### {step_name}")
    if explanation:
        st.write(explanation)
    
    # Display formula
    if formula_latex:
        st.latex(formula_latex)
    
    # Display result
    if result_value is not None:
        if isinstance(result_value, (list, tuple, np.ndarray)):
             if len(result_value) > 4:
                res_str = f"Array of shape {np.shape(result_value)}"
             else:
                res_str = ", ".join([f"{v:.4g}" for v in result_value])
             st.markdown(f"**Result:** `[{res_str}] {unit}`")
        else:
            st.markdown(f"**Result:** `{result_value:.4g} {unit}`")
            
    st.caption(f"📍 Reference: {page_ref}")
    st.divider()

def section_header(title, page_range):
    st.markdown(f"## {title}")
    st.markdown(f"*Covering content from Pages {page_range}*")
    st.divider()

# ==========================================
# CASE 1: FLEXURAL-TORSIONAL BUCKLING
# ==========================================
def case_1_flexural_torsional():
    section_header("1. Flexural-Torsional Buckling of Open Thin-Walled Columns", "38 - 51")

    # --- Inputs (Page 38) ---
    st.sidebar.header("Inputs (Ref: Page 38)")
    L_input = st.sidebar.number_input("Column Length (L) [m]", 1.0)
    E = st.sidebar.number_input("Young's Modulus (E) [Pa]", 70e9, format="%.2e")
    G = st.sidebar.number_input("Shear Modulus (mu) [Pa]", 30e9, format="%.2e")
    h = st.sidebar.number_input("Web Height (h) [m]", 0.1)
    b = st.sidebar.number_input("Flange Width (b) [m]", 0.1)
    t = st.sidebar.number_input("Thickness (t) [m]", 0.002, format="%.4f")

    # --- Calculations ---
    # Page 39: Centroid
    y_c_prime = (2 * b * t * (b / 2)) / (2 * b * t + h * t)
    
    # Page 39: Inertia
    I_yy = (t * h**3 / 12) * (1 + (6 * b) / h)
    term1 = 2 * (t * b**3) / 12
    term2 = 2 * t * b * ((b/2) - y_c_prime)**2
    term3 = h * t * y_c_prime**2
    I_zz = term1 + term2 + term3
    
    # Page 40: Shear Center
    # Analytical result for C-section relative to web
    y_s_prime_web = -(3*b**2)/(h + 6*b)
    y_S = y_s_prime_web - y_c_prime
    
    # Page 42: Torsion & Warping
    C_torsion = (G * t**3 / 3) * (2 * b + h)
    A_sect = 2 * b * t + h * t
    I_p_S = I_yy + I_zz + A_sect * (y_S**2)
    
    # Page 48: Warping Constant (Approx for C-section)
    C_gamma_geo = (h**2 * b**3 * t / 12) * ((3 * b + 2 * h) / (6 * b + h))
    C_gamma = C_gamma_geo * E

    # Solver for P_cr
    def solve_P_cr(length):
        P_y = (np.pi**2 * E * I_yy) / (length**2)
        P_z = (np.pi**2 * E * I_zz) / (length**2)
        P_theta = (A_sect / I_p_S) * ((C_gamma * np.pi**2 / length**2) + C_torsion)
        
        # Quadratic for Coupled Mode (Page 51)
        # Roots of: (P - Pz) * [ P^2((A/Ip)*ys^2 - 1) + P(P_theta + Py) - P_theta*Py ] = 0
        coeff_a = (A_sect / I_p_S) * y_S**2 - 1
        coeff_b = P_theta + P_y
        coeff_c = -P_theta * P_y
        
        discriminant = coeff_b**2 - 4 * coeff_a * coeff_c
        roots = [P_z]
        if discriminant >= 0:
            r1 = (-coeff_b + np.sqrt(discriminant)) / (2 * coeff_a)
            r2 = (-coeff_b - np.sqrt(discriminant)) / (2 * coeff_a)
            roots.extend([r1, r2])
        
        pos_roots = [r for r in roots if r > 0]
        return min(pos_roots) if pos_roots else 0.0

    P_critical_val = solve_P_cr(L_input)

    # --- UI Tabs ---
    tab1, tab2 = st.tabs(["Calculation Log", "Stability Plots"])

    with tab1:
        st.subheader("Step-by-Step Calculation Log")
        log_step("1. Centroid Position", r"y_C' = \frac{2bt(b/2)}{2bt + ht}", y_c_prime, "m", "Page 39")
        log_step("2. Second Moments of Area", r"I_{yy}, I_{zz}", [I_yy, I_zz], "m^4", "Page 39")
        log_step("3. Shear Center", r"y_S = y_S' - y_C'", y_S, "m", "Page 40")
        log_step("4. Torsion Properties", r"C, I_p^S", [C_torsion, I_p_S], "units", "Page 42")
        log_step("5. Warping Constant", r"C^{\Gamma}", C_gamma, "N*m^4", "Page 48")
        log_step("6. Critical Load (L={})".format(L_input), r"P_{CR}", P_critical_val, "N", "Page 51")

    with tab2:
        st.subheader("Buckling Load vs Column Length (Pages 49-51)")
        L_range = np.linspace(0.5, 5.0, 50)
        P_curve = [solve_P_cr(l) for l in L_range]
        P_z_curve = [(np.pi**2 * E * I_zz) / (l**2) for l in L_range]

        fig, ax = plt.subplots()
        ax.plot(L_range, np.array(P_curve)/1000, label="Coupled Mode (P_cr)", linewidth=2, color='blue')
        ax.plot(L_range, np.array(P_z_curve)/1000, '--', label="Pure Flexural (z) Mode", color='orange')
        
        ax.set_xlabel("Length L [m]")
        ax.set_ylabel("Critical Load P [kN]")
        ax.set_title("Stability Curve")
        ax.grid(True)
        ax.legend()
        ax.axvline(L_input, color='red', linestyle=':', label='Current L')
        st.pyplot(fig)

# ==========================================
# CASE 2: PLATE BUCKLING
# ==========================================
def case_2_plate_buckling():
    section_header("2. Buckling of Thin Plates & Secondary Buckling", "60 - 74")
    
    st.sidebar.header("Inputs")
    # Defaulting to values similar to the L-section flange example on Page 71 or general plate
    a = st.sidebar.number_input("Plate Length (a) [m]", 0.3)
    b = st.sidebar.number_input("Plate Width (b) [m]", 0.1)
    t = st.sidebar.number_input("Thickness (t) [m]", 0.002)
    E = st.sidebar.number_input("Young's Modulus (E) [Pa]", 70e9, format="%.2e")
    nu = st.sidebar.number_input("Poisson Ratio (nu)", 0.3)
    
    bc_type = st.sidebar.selectbox("Boundary Conditions", 
                                   ["Simply Supported (SSSS)", "One Edge Free (SSSF - L-Section Example)"])

    st.subheader("Analysis & Visualization")

    tab1, tab2, tab3 = st.tabs(["Calculation Log", "Buckling Coefficients (Page 61)", "3D Mode Shapes (Page 60)"])

    # --- Calculations ---
    # Critical Stress Base: D terms
    # Sigma_cr = k * (pi^2 E / 12(1-nu^2)) * (t/b)^2
    factor_D = (np.pi**2 * E) / (12 * (1 - nu**2)) * (t / b)**2
    aspect_ratio = a / b

    if bc_type == "Simply Supported (SSSS)":
        # k = (m/alpha + alpha/m)^2
        # We check m=1 to 5 to find min
        m_vals = range(1, 6)
        k_vals = [(m/aspect_ratio + aspect_ratio/m)**2 for m in m_vals]
        k_min = min(k_vals)
        m_crit = m_vals[k_vals.index(k_min)]
        ref_page = "Page 61"
    else:
        # L-section example (Page 71)
        # k approx 0.43 for long plates (a/b=3 in PDF example)
        k_min = 0.43
        m_crit = 1
        ref_page = "Page 71 (L-section example)"

    sigma_cr = k_min * factor_D
    P_cr_plate = sigma_cr * (b * t)

    with tab1:
        log_step("1. Aspect Ratio (a/b)", r"a/b", aspect_ratio, "-", "Page 60")
        log_step("2. Buckling Coefficient (k)", r"k_{min}", k_min, "-", ref_page, f"Critical m={m_crit}")
        log_step("3. Critical Stress", r"\sigma_{cr} = k \frac{\pi^2 E}{12(1-\nu^2)} (\frac{t}{b})^2", 
                 sigma_cr / 1e6, "MPa", "Page 61")
        log_step("4. Critical Load (Total)", r"P_{cr} = \sigma_{cr} bt", P_cr_plate, "N", "-")

    with tab2:
        st.subheader("Buckling Coefficient Chart (Ref: Page 61)")
        # Plot k vs a/b for m=1..4 (SSSS)
        ratios = np.linspace(0.1, 4.0, 200)
        fig, ax = plt.subplots()
        
        colors = ['blue', 'green', 'red', 'purple']
        for m in range(1, 5):
            k_curve = (m / ratios + ratios / m)**2
            ax.plot(ratios, k_curve, label=f"m={m}", color=colors[m-1])
        
        ax.set_ylim(0, 10)
        ax.set_xlim(0, 4)
        ax.set_xlabel("Aspect Ratio (a/b)")
        ax.set_ylabel("Buckling Coefficient k")
        ax.set_title("Buckling Coefficients (SSSS Plate)")
        ax.axvline(aspect_ratio, color='black', linestyle=':', label='Current a/b')
        if bc_type == "Simply Supported (SSSS)":
            ax.scatter([aspect_ratio], [k_min], color='black', zorder=5)
        
        ax.legend()
        ax.grid(True)
        st.pyplot(fig)

    with tab3:
        st.subheader(f"Buckling Mode Shape (m={m_crit}, n=1) (Ref: Page 60)")
        x = np.linspace(0, a, 50)
        y = np.linspace(0, b, 50)
        X, Y = np.meshgrid(x, y)
        
        # SSSS Shape: sin(m pi x / a) sin(n pi y / b)
        # For SSSF, the shape is more complex, but we approximate visual for "One Edge Free" as similar but shifting
        W = np.sin(m_crit * np.pi * X / a) * np.sin(1 * np.pi * Y / b)
        
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(X, Y, W, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.set_title(f"Mode Shape m={m_crit}")
        st.pyplot(fig)

# ==========================================
# CASE 3: TRANSVERSE COLUMNS
# ==========================================
def case_3_transverse_columns():
    section_header("3. Transversely Loaded Columns (Annex 1)", "76 - 82")
    
    st.sidebar.header("Inputs (Ref: Page 76)")
    L = st.sidebar.number_input("Length (L) [m]", 1.0)
    E = st.sidebar.number_input("Young's Modulus (E) [Pa]", 70e9, format="%.2e")
    I = st.sidebar.number_input("Inertia (I) [m^4]", 1.17e-6, format="%.2e")
    f_z = st.sidebar.number_input("Transverse Load (f_z) [N/m]", 1000.0)
    P_input = st.sidebar.number_input("Axial Load (P) [N]", 100000.0)

    # --- Calcs ---
    P_cr = (np.pi**2 * E * I) / (L**2)
    omega = np.sqrt(P_input / (E * I))
    
    # Page 80: Deflection
    sec_term = 1 / np.cos(omega * L / 2)
    u_max = (f_z / (omega**2 * P_input)) * (sec_term - 1) - (f_z * L**2) / (8 * P_input)
    
    # Page 82: Moment
    M_max = (f_z / omega**2) * (sec_term - 1)

    tab1, tab2 = st.tabs(["Calculation Log", "Beam-Column Plots"])
    
    with tab1:
        log_step("1. Euler Critical Load", r"P_{CR}", P_cr, "N", "Page 80")
        if P_input < P_cr:
            log_step("2. Max Deflection (x=L/2)", r"u_{max}", u_max, "m", "Page 80")
            log_step("3. Max Moment", r"M_{max}", M_max, "N*m", "Page 82")
        else:
            st.error("Axial Load Exceeds Critical Load! Structure Unstable.")

    with tab2:
        st.subheader("Magnification Effect (Ref: Page 81)")
        P_vals = np.linspace(100, 0.95 * P_cr, 50)
        u_vals = []
        for p in P_vals:
            w = np.sqrt(p / (E * I))
            val = (f_z / (w**2 * p)) * ((1 / np.cos(w * L / 2)) - 1) - (f_z * L**2) / (8 * p)
            u_vals.append(val)
            
        fig, ax = plt.subplots()
        ax.plot(P_vals/1000, np.array(u_vals)*1000, label="Max Deflection", color='blue')
        ax.axvline(P_cr/1000, color='red', linestyle='--', label="P_cr (Euler)")
        ax.axvline(P_input/1000, color='green', linestyle=':', label="Current Load")
        
        ax.set_xlabel("Axial Load P [kN]")
        ax.set_ylabel("Deflection [mm]")
        ax.set_title("Beam-Column Response")
        ax.legend()
        ax.grid(True)
        st.pyplot(fig)

# ==========================================
# CASE 4: STIFFENER/WEB
# ==========================================
def case_4_stiffener_web():
    section_header("4. Buckling of Stiffener/Web Constructions (Annex 2)", "89 - 93")
    
    st.sidebar.header("Inputs (Ref: Page 90)")
    d = st.sidebar.number_input("Web Depth (d) [m]", 0.4)
    b_stiff = st.sidebar.number_input("Stiffener Spacing (b) [m]", 0.3)
    t = st.sidebar.number_input("Web Thickness (t) [m]", 0.002)
    T_force = st.sidebar.number_input("Shear Force (T) [N]", 5000.0)
    E = st.sidebar.number_input("Young's Modulus (E) [Pa]", 70e9, format="%.2e")
    
    A_F = st.sidebar.number_input("Flange Area (A_F) [m^2]", 350e-6, format="%.2e")
    A_S = st.sidebar.number_input("Stiffener Area (A_S) [m^2]", 300e-6, format="%.2e")
    I_stiff = st.sidebar.number_input("Stiffener Inertia (I_xx) [m^4]", 2000e-12, format="%.2e")
    
    tab1, tab2 = st.tabs(["Calculation Log", "Stability Map"])
    
    # --- Calcs ---
    # Page 90: Diagonal Tension Angle
    num = 1 + (t * d) / (2 * A_F)
    den = 1 + (t * b_stiff) / A_S
    alpha_rad = np.arctan((num / den)**0.25)
    
    # Page 92: Stiffener Load
    P_stiff = -T_force * (b_stiff / d) * np.tan(alpha_rad)
    
    # Page 92/93: Effective Length & Buckling
    ratio = b_stiff / d
    if ratio < 1.5:
        le = d / np.sqrt(4 - 2 * ratio)
    else:
        le = d
    P_cr_stiff = (np.pi**2 * E * I_stiff) / (le**2)

    with tab1:
        log_step("1. Diagonal Tension Angle", r"\alpha", np.degrees(alpha_rad), "deg", "Page 90")
        log_step("2. Load on Stiffener", r"P_{stiff} = -T \frac{b}{d} \tan \alpha", abs(P_stiff), "N", "Page 92")
        log_step("3. Critical Buckling Load", r"P_{CR}", P_cr_stiff, "N", "Page 93")
        
        if abs(P_stiff) < P_cr_stiff:
            st.success("Stiffener is SAFE against buckling.")
        else:
            st.error("Stiffener FAILURE predicted.")

    with tab2:
        st.subheader("Stiffener Stability Map (Ref: Page 92)")
        b_range = np.linspace(0.1, 1.0, 50)
        P_cr_curve = []
        P_app_curve = []
        
        for b_val in b_range:
            # Re-calc local le
            if b_val < 1.5 * d:
                le_local = d / np.sqrt(4 - 2 * (b_val / d))
            else:
                le_local = d
            P_cr_curve.append((np.pi**2 * E * I_stiff) / (le_local**2))
            
            # Re-calc local alpha
            den_local = 1 + (t * b_val) / A_S
            alpha_local = np.arctan((num / den_local)**0.25)
            P_app_curve.append(abs(-T_force * (b_val / d) * np.tan(alpha_local)))

        fig, ax = plt.subplots()
        ax.plot(b_range, P_cr_curve, label="P_critical (Buckling)", color='green')
        ax.plot(b_range, P_app_curve, label="P_applied (Load)", color='red')
        
        # Shade failure region
        ax.fill_between(b_range, P_cr_curve, P_app_curve, where=(np.array(P_app_curve) > np.array(P_cr_curve)), 
                        color='red', alpha=0.1, label='Failure Zone')
        
        ax.set_xlabel("Stiffener Spacing b [m]")
        ax.set_ylabel("Load [N]")
        ax.set_title("Stiffener Stability Map")
        ax.axvline(b_stiff, color='blue', linestyle=':', label="Current b")
        ax.legend()
        ax.grid(True)
        st.pyplot(fig)

# ==========================================
# MAIN APP
# ==========================================
def main():
    st.set_page_config(page_title="Structure Instabilities Analysis", layout="wide")
    st.title("Aerospace Structure Instabilities Tool")
    st.markdown("""
    Interactive verification tool for **StructAeroInstabilities.pdf**.
    Select a case from the sidebar to begin analysis.
    """)
    st.divider()

    menu = [
        "1. Flexural-Torsional Buckling (Pg 38)", 
        "2. Plate Buckling (Pg 60)", 
        "3. Transversely Loaded Cols (Pg 76)", 
        "4. Stiffener/Web Buckling (Pg 90)"
    ]
    
    choice = st.sidebar.radio("Select Analysis Case:", menu)

    if "1." in choice:
        case_1_flexural_torsional()
    elif "2." in choice:
        case_2_plate_buckling()
    elif "3." in choice:
        case_3_transverse_columns()
    elif "4." in choice:
        case_4_stiffener_web()

if __name__ == "__main__":
    main()