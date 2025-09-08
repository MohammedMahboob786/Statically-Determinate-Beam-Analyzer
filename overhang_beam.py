# overhang.py
import ast
import numpy as np
import plotly.graph_objects as go
import streamlit as st

import numpy as np

def find_contraflexure_points(
    X, M,
    zero_tol: float = 1e-6,
    min_gap: float = 0.01,
    prefer: str = "bracket",   # "bracket" | "slope"
    min_amp: float = 0.0,      # ignore crossings whose bracket magnitude < min_amp
    return_details: bool = False
):
    """
    Returns de-duplicated x-positions where bending moment crosses 0, preferring
    the crossing with higher magnitude within each proximity cluster.

    Args:
        X, M: arrays
        zero_tol: treat |M| < zero_tol as zero (round-off guard)
        min_gap: merge roots closer than this distance (meters)
        prefer:  "bracket" -> keep root with largest max(|y1|,|y2|)
                 "slope"   -> keep root with largest |(y2 - y1)/(x2 - x1)|
        min_amp: ignore crossings whose bracket magnitude < min_amp
        return_details: if True, also returns a dict with metrics used

    Returns:
        roots (list[float]) or (roots, details) if return_details=True
    """
    X = np.asarray(X, dtype=np.float64)
    M = np.asarray(M, dtype=np.float64)

    # 1) snap near-zeros
    M_snapped = M.copy()
    M_snapped[np.abs(M_snapped) < zero_tol] = 0.0

    # 2) sign-change indices
    sc = np.where(np.sign(M_snapped[:-1]) * np.sign(M_snapped[1:]) < 0)[0]

    # 3) gather candidate roots + significance metrics
    cand_x = []
    bracket_amp = []  # max(|y1|, |y2|)
    slope_mag = []    # |dM/dx| across the segment
    src_idx = []      # (k) index to segment start

    for k in sc:
        x1, x2 = X[k], X[k+1]
        y1, y2 = M_snapped[k], M_snapped[k+1]
        if y2 == y1:
            continue
        # linear interpolation root
        xr = float(x1 - y1 * (x2 - x1) / (y2 - y1))
        cand_x.append(xr)
        bracket_amp.append(float(max(abs(y1), abs(y2))))
        slope_mag.append(float(abs((y2 - y1) / (x2 - x1))))
        src_idx.append(k)

    # 4) include exact on-grid zeros as candidates (estimate metrics from neighbors)
    z_idx = np.where(M_snapped == 0.0)[0]
    for i in z_idx:
        xi = float(X[i])
        # finite-diff estimate around the point
        if 0 < i < len(M_snapped) - 1:
            yL, yR = M_snapped[i-1], M_snapped[i+1]
            xL, xR = X[i-1], X[i+1]
            amp = float(max(abs(yL), abs(yR)))
            slp = float(abs((yR - yL) / (xR - xL))) if (xR != xL) else 0.0
        else:
            # end points: approximate with one-sided difference if possible
            if i == 0 and len(M_snapped) > 1:
                amp = float(abs(M_snapped[1]))
                slp = float(abs((M_snapped[1] - M_snapped[0]) / (X[1] - X[0])))
            elif i == len(M_snapped) - 1 and len(M_snapped) > 1:
                amp = float(abs(M_snapped[-2]))
                slp = float(abs((M_snapped[-1] - M_snapped[-2]) / (X[-1] - X[-2])))
            else:
                amp, slp = 0.0, 0.0
        cand_x.append(xi)
        bracket_amp.append(amp)
        slope_mag.append(slp)
        src_idx.append(i)

    # Nothing found
    if not cand_x:
        return ([], {}) if return_details else []

    # Filter by minimum amplitude (drop insignificant near-zero noise)
    cand = np.array(list(zip(cand_x, bracket_amp, slope_mag, src_idx)), dtype=float)
    if min_amp > 0.0:
        cand = cand[cand[:,1] >= min_amp]  # keep where bracket_amp >= min_amp
        if cand.size == 0:
            return ([], {}) if return_details else []

    # 5) sort by x and cluster by proximity (min_gap)
    order = np.argsort(cand[:,0])
    cand = cand[order]

    clusters = []
    cur = [cand[0]]
    for row in cand[1:]:
        if row[0] - cur[-1][0] < min_gap:
            cur.append(row)
        else:
            clusters.append(np.array(cur))
            cur = [row]
    clusters.append(np.array(cur))

    # 6) per cluster, choose the "most significant" root
    roots = []
    chosen_idx = []  # indices in X/M related to chosen roots (best-effort)
    for cl in clusters:
        if prefer == "slope":
            # choose by highest slope magnitude
            j = int(np.argmax(cl[:,2]))  # slope_mag
        else:
            # default: choose by largest bracket amplitude
            j = int(np.argmax(cl[:,1]))  # bracket_amp
        roots.append(float(cl[j,0]))
        chosen_idx.append(int(cl[j,3]))

    if return_details:
        details = {
            "roots": roots,
            "chosen_segment_index": chosen_idx,
            "bracket_amp": [float(x) for x in cand[:,1]],
            "slope_mag": [float(x) for x in cand[:,2]],
            "all_candidates": [float(x) for x in cand[:,0]],
        }
        return roots, details
    return roots

class OverhangBeamApp:
    def __init__(self):
        # UI state
        self.span = 0.0
        self.A = 0.0
        self.B = 0.0
        self.E = 0
        self.I = 0.0

        # Raw text inputs
        self.pointloads_txt = "[[]]"
        self.pointmoments_txt = "[[]]"
        self.distributedloads_txt = "[[]]"
        self.linearloads_txt = "[[]]"
        self.trapezoidalloads_txt = "[[]]"

        # Parsed arrays
        self.pointloads = None
        self.pointmoments = None
        self.distributedloads = None
        self.linearloads = None
        self.trapezoidalloads = None

        # Computation state
        self.delta = 5e-5
        self.X = None
        self.reactions = np.array([0.0, 0.0])  # Va, Vb
        self.shearForce = None
        self.bendingMoment = None

        # Records
        self.PL_record = None  # (Va, Vb) per load
        self.PM_record = None
        self.UDL_record = None
        self.LDL_record = None
        self.TDL_record = None

    # ---------- UI ----------
    @staticmethod
    def _convert_to_list_of_lists(input_str: str):
        try:
            result = ast.literal_eval(input_str)
            if isinstance(result, list) and all(isinstance(i, list) for i in result):
                return result
            st.error("Input is not a valid list of lists.")
            return None
        except (ValueError, SyntaxError):
            st.error("Invalid input format.")
            return None

    def render_inputs(self):
        st.title("Overhang Beam Analysis Tool")

        st.markdown("#### Input Parameters")
        self.span = st.number_input("📏 Total length of the beam (in meters):", min_value=0.0, format="%.2f")
        st.markdown("*Units: meters (m) for lengths, kilonewtons (kN) for forces*")

        self.A = st.number_input("🔧 Location of First Support (A) from Left End (in meters):",
                                 min_value=0.0, max_value=self.span, format="%.2f")
        self.B = st.number_input("🔧 Location of Second Support (B) from Left End (in meters):",
                                 min_value=self.A, max_value=self.span, format="%.2f")

        st.markdown("---")
        st.markdown("#### 🔹 Load Inputs")
        self.pointloads_txt = st.text_input(
            "📌 Point Loads",
            help="[[x, Fy], ...]  ➡️ e.g., [[2, -10], [4.5, -15]]",
            value=self.pointloads_txt,
        )
        self.pointmoments_txt = st.text_input(
            "♻️ Point Moments",
            help="[[x, M], ...]  ➡️ e.g., [[1.5, 5], [3, -8]]",
            value=self.pointmoments_txt,
        )
        self.distributedloads_txt = st.text_input(
            "🌫️ Uniform Distributed Loads (UDL)",
            help="[[x_start, x_end, w], ...]  ➡️ e.g., [[0, 2, -5]]",
            value=self.distributedloads_txt,
        )
        self.linearloads_txt = st.text_input(
            "📈 Linear Varying Loads (Triangular)",
            help="[[x_start, x_end, w_start, w_end], ...]  ➡️ e.g., [[1, 3, 0, -10]]",
            value=self.linearloads_txt,
        )
        self.trapezoidalloads_txt = st.text_input(
            "🧱 Trapezoidal Loads",
            help="[[x_start, x_end, w_start, w_end], ...]  ➡️ e.g., [[2, 5, -4, -8]]",
            value=self.trapezoidalloads_txt,
        )

        st.markdown("---")
        st.markdown("#### ⚙️ Material & Beam Properties")
        self.E = st.number_input(
            "🧪 Young's Modulus (E) in Pascals (N/m²):",
            help="Example (steel): 200000000 (200 GPa).",
            min_value=0, max_value=1_000_000_000, value=200_000_000, step=1_000_000,
        )
        self.I = st.number_input(
            "📏 Moment of Inertia (I) in m⁴:",
            help="Example: 0.0020833 for a rectangular section (b=0.3m, h=0.5m).",
            min_value=0.0, max_value=1.0, value=0.0020833, step=0.0000001, format="%.7f",
        )

    def _parse_inputs(self) -> bool:
        self.pointloads = np.array(self._convert_to_list_of_lists(self.pointloads_txt))
        self.pointmoments = np.array(self._convert_to_list_of_lists(self.pointmoments_txt))
        self.distributedloads = np.array(self._convert_to_list_of_lists(self.distributedloads_txt))
        self.linearloads = np.array(self._convert_to_list_of_lists(self.linearloads_txt))
        self.trapezoidalloads = np.array(self._convert_to_list_of_lists(self.trapezoidalloads_txt))

        if any(x is None for x in [
            self.pointloads, self.pointmoments, self.distributedloads, self.linearloads, self.trapezoidalloads
        ]):
            return False

        self.X = np.arange(0, self.span + self.delta, self.delta)
        self.shearForce = np.empty([0, len(self.X)])
        self.bendingMoment = np.empty([0, len(self.X)])
        return True

    # ---------- Reactions ----------
    def _reactions_PL(self, n):
        xp = self.pointloads[n, 0]
        fy = self.pointloads[n, 1]
        la_p = self.A - xp
        mp = fy * la_p
        la_vb = self.B - self.A
        Vb = mp / la_vb
        Va = -fy - Vb
        return Va, Vb

    def _reactions_PM(self, n):
        m = self.pointmoments[n, 1]
        la_vb = self.B - self.A
        Vb = m / la_vb
        Va = -Vb
        return Va, Vb

    def _reactions_UDL(self, n):
        xStart, xEnd, fy = self.distributedloads[n]
        fy_Res = fy * (xEnd - xStart)
        x_Res = xStart + 0.5 * (xEnd - xStart)
        la_p = self.A - x_Res
        mp = fy_Res * la_p
        la_vb = self.B - self.A
        Vb = mp / la_vb
        Va = -fy_Res - Vb
        return Va, Vb

    def _reactions_LDL(self, n):
        xStart, xEnd, fy_start, fy_end = self.linearloads[n]
        if abs(fy_start) > 0:
            fy_Res = 0.5 * fy_start * (xEnd - xStart)
            x_Res = xStart + (1 / 3) * (xEnd - xStart)
        else:
            fy_Res = 0.5 * fy_end * (xEnd - xStart)
            x_Res = xStart + (2 / 3) * (xEnd - xStart)
        la_p = self.A - x_Res
        mp = fy_Res * la_p
        la_vb = self.B - self.A
        Vb = mp / la_vb
        Va = -fy_Res - Vb
        return Va, Vb

    def _reactions_TDL(self, n):
        xStart, xEnd, fy_start, fy_end = self.trapezoidalloads[n]
        fy_Res = 0.5 * (xEnd - xStart) * (fy_start + fy_end)
        x_Res = xStart + (1 / 3) * (((xEnd - xStart) * (fy_start + 2 * fy_end)) / (fy_start + fy_end))
        la_p = self.A - x_Res
        mp = fy_Res * la_p
        la_vb = self.B - self.A
        Vb = mp / la_vb
        Va = -fy_Res - Vb
        return Va, Vb

    def _accumulate_reactions(self):
        nPL = len(self.pointloads[0])
        nPM = len(self.pointmoments[0])
        nUDL = len(self.distributedloads[0])
        nLDL = len(self.linearloads[0])
        nTDL = len(self.trapezoidalloads[0])

        self.PL_record = np.empty([0, 2])
        self.PM_record = np.empty([0, 2])
        self.UDL_record = np.empty([0, 2])
        self.LDL_record = np.empty([0, 2])
        self.TDL_record = np.empty([0, 2])

        if nPL > 0:
            for n, _ in enumerate(self.pointloads):
                va, vb = self._reactions_PL(n)
                self.PL_record = np.append(self.PL_record, [[va, vb]], axis=0)
                self.reactions[0] += va
                self.reactions[1] += vb

        if nPM > 0:
            for n, _ in enumerate(self.pointmoments):
                va, vb = self._reactions_PM(n)
                self.PM_record = np.append(self.PM_record, [[va, vb]], axis=0)
                self.reactions[0] += va
                self.reactions[1] += vb

        if nUDL > 0:
            for n, _ in enumerate(self.distributedloads):
                va, vb = self._reactions_UDL(n)
                self.UDL_record = np.append(self.UDL_record, [[va, vb]], axis=0)
                self.reactions[0] += va
                self.reactions[1] += vb

        if nLDL > 0:
            for n, _ in enumerate(self.linearloads):
                va, vb = self._reactions_LDL(n)
                self.LDL_record = np.append(self.LDL_record, [[va, vb]], axis=0)
                self.reactions[0] += va
                self.reactions[1] += vb

        if nTDL > 0:
            for n, _ in enumerate(self.trapezoidalloads):
                va, vb = self._reactions_TDL(n)
                self.TDL_record = np.append(self.TDL_record, [[va, vb]], axis=0)
                self.reactions[0] += va
                self.reactions[1] += vb

        st.write(f'The vertical reaction at A is {round(self.reactions[0], 2)} kN')
        st.write(f'The vertical reaction at B is {round(self.reactions[1], 2)} kN')

    # ---------- Shear/Moment contributors ----------
    def _sm_PL(self, n):
        xp, fy = self.pointloads[n, 0], self.pointloads[n, 1]
        Va, Vb = self.PL_record[n]
        S = np.zeros(len(self.X))
        M = np.zeros(len(self.X))
        for i, x in enumerate(self.X):
            shear = 0.0
            moment = 0.0
            if x > self.A:
                shear += Va
                moment -= Va * (x - self.A)
            if x >= self.B:  # overhang: reaction at B applies to the right of B
                shear += Vb
                moment -= Vb * (x - self.B)
            if x > xp:
                shear += fy
                moment -= fy * (x - xp)
            S[i] = shear
            M[i] = moment
        return S, M

    def _sm_PM(self, n):
        xm, m = self.pointmoments[n, 0], self.pointmoments[n, 1]
        Va, Vb = self.PM_record[n]
        S = np.zeros(len(self.X))
        M = np.zeros(len(self.X))
        for i, x in enumerate(self.X):
            shear = 0.0
            moment = 0.0
            if x > self.A:
                shear += Va
                moment -= Va * (x - self.A)
            if x >= self.B:
                shear += Vb
                moment -= Vb * (x - self.B)
            if x > xm:
                moment -= m
            S[i] = shear
            M[i] = moment
        return S, M

    def _sm_UDL(self, n):
        xStart, xEnd, fy = self.distributedloads[n]
        Va, Vb = self.UDL_record[n]
        S = np.zeros(len(self.X))
        M = np.zeros(len(self.X))
        for i, x in enumerate(self.X):
            shear = 0.0
            moment = 0.0
            if x > self.A:
                shear += Va
                moment -= Va * (x - self.A)
            if x >= self.B:
                shear += Vb
                moment -= Vb * (x - self.B)
            if x > xStart and x <= xEnd:
                dx = x - xStart
                shear += fy * dx
                moment -= fy * dx * 0.5 * dx
            elif x > xEnd:
                L = xEnd - xStart
                shear += fy * L
                moment -= fy * L * (x - xStart - 0.5 * L)
            S[i] = shear
            M[i] = moment
        return S, M

    def _sm_LDL(self, n):
        xStart, xEnd, fy_start, fy_end = self.linearloads[n]
        Va, Vb = self.LDL_record[n]
        S = np.zeros(len(self.X))
        M = np.zeros(len(self.X))
        for i, x in enumerate(self.X):
            shear = 0.0
            moment = 0.0
            if x > self.A:
                shear += Va
                moment -= Va * (x - self.A)
            if x >= self.B:
                shear += Vb
                moment -= Vb * (x - self.B)

            if x > xStart and x <= xEnd:
                x_base = x - xStart
                if abs(fy_start) > 0:
                    f_cut = fy_start - x_base * (fy_start / (xEnd - xStart))
                    R1 = 0.5 * x_base * (fy_start - f_cut)
                    R2 = x_base * f_cut
                    shear += R1 + R2
                    moment -= R1 * (2 / 3) * x_base - R2 * (x_base / 2)
                else:
                    f_cut = fy_end * (x_base / (xEnd - xStart))
                    R = 0.5 * x_base * f_cut
                    shear += R
                    moment -= R * (x_base / 3)
            elif x > xEnd:
                if abs(fy_start) > 0:
                    R = 0.5 * fy_start * (xEnd - xStart)
                    xr = xStart + (1 / 3) * (xEnd - xStart)
                else:
                    R = 0.5 * fy_end * (xEnd - xStart)
                    xr = xStart + (2 / 3) * (xEnd - xStart)
                shear += R
                moment -= R * (x - xr)

            S[i] = shear
            M[i] = moment
        return S, M

    def _sm_TDL(self, n):
        xStart, xEnd, fy_start, fy_end = self.trapezoidalloads[n]
        Va, Vb = self.TDL_record[n]
        S = np.zeros(len(self.X))
        M = np.zeros(len(self.X))
        for i, x in enumerate(self.X):
            shear = 0.0
            moment = 0.0
            if x > self.A:
                shear += Va
                moment -= Va * (x - self.A)
            if x >= self.B:
                shear += Vb
                moment -= Vb * (x - self.B)

            if x > xStart and x <= xEnd:
                x_base = x - xStart
                if abs(fy_start) > abs(fy_end):
                    f_cut = fy_start - x_base * ((fy_start - fy_end) / (xEnd - xStart))
                    R1 = 0.5 * x_base * (fy_start - f_cut)
                    la_R1 = (2 / 3) * x_base
                    R2 = x_base * f_cut
                    la_R2 = 0.5 * x_base
                else:
                    f_cut = fy_start + x_base * ((fy_end - fy_start) / (xEnd - xStart))
                    R1 = 0.5 * x_base * (f_cut - fy_start)
                    la_R1 = (1 / 3) * x_base
                    R2 = x_base * fy_start
                    la_R2 = 0.5 * x_base
                shear += R1 + R2
                moment -= R1 * la_R1 - R2 * la_R2
            elif x > xEnd:
                R = 0.5 * (xEnd - xStart) * (fy_start + fy_end)
                xr = xStart + (1 / 3) * (((xEnd - xStart) * (fy_start + 2 * fy_end)) / (fy_start + fy_end))
                shear += R
                moment -= R * (x - xr)

            S[i] = shear
            M[i] = moment
        return S, M

    def _assemble_shear_moment(self):
        nPL = len(self.pointloads[0])
        nPM = len(self.pointmoments[0])
        nUDL = len(self.distributedloads[0])
        nLDL = len(self.linearloads[0])
        nTDL = len(self.trapezoidalloads[0])

        if nPL > 0:
            for n, _ in enumerate(self.pointloads):
                S, M = self._sm_PL(n)
                self.shearForce = np.append(self.shearForce, [S], axis=0)
                self.bendingMoment = np.append(self.bendingMoment, [M], axis=0)

        if nPM > 0:
            for n, _ in enumerate(self.pointmoments):
                S, M = self._sm_PM(n)
                self.shearForce = np.append(self.shearForce, [S], axis=0)
                self.bendingMoment = np.append(self.bendingMoment, [M], axis=0)

        if nUDL > 0:
            for n, _ in enumerate(self.distributedloads):
                S, M = self._sm_UDL(n)
                self.shearForce = np.append(self.shearForce, [S], axis=0)
                self.bendingMoment = np.append(self.bendingMoment, [M], axis=0)

        if nLDL > 0:
            for n, _ in enumerate(self.linearloads):
                S, M = self._sm_LDL(n)
                self.shearForce = np.append(self.shearForce, [S], axis=0)
                self.bendingMoment = np.append(self.bendingMoment, [M], axis=0)

        if nTDL > 0:
            for n, _ in enumerate(self.trapezoidalloads):
                S, M = self._sm_TDL(n)
                self.shearForce = np.append(self.shearForce, [S], axis=0)
                self.bendingMoment = np.append(self.bendingMoment, [M], axis=0)

    # ---------- Plotters ----------
    def _plot_shear(self):
        layout = go.Layout(
            title={'text': "Shear Force Diagram", 'y': 0.85, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'},
            titlefont=dict(size=15),
            yaxis=dict(title='Shear Force (kN)'),
            xaxis=dict(title='Distance (m)', range=[-1, self.span + 1]),
            showlegend=False,
        )
        line = go.Scatter(
            x=self.X,
            y=sum(self.shearForce) if len(self.shearForce) else np.zeros_like(self.X),
            mode='lines',
            name='Shear Force',
            fill='tonexty',
            line_color='green',
            fillcolor='rgba(0, 255, 0, 0.1)',
        )
        axis = go.Scatter(x=[0, self.span], y=[0, 0], mode='lines', line_color='black')
        st.plotly_chart(go.Figure(data=[line, axis], layout=layout), use_container_width=True)

    def _plot_moment(self):
        layout = go.Layout(
            title={'text': "Bending Moment Diagram", 'y': 0.85, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'},
            titlefont=dict(size=15),
            yaxis=dict(title='Bending Moment (kNm)', autorange="reversed"),
            xaxis=dict(title='Distance (m)', range=[-1, self.span + 1]),
            showlegend=False,
        )
        line = go.Scatter(
            x=self.X,
            y=-(sum(self.bendingMoment) if len(self.bendingMoment) else np.zeros_like(self.X)),
            mode='lines',
            name='Bending Moment',
            fill='tonexty',
            line_color='red',
            fillcolor='rgba(255, 0, 0, 0.1)',
        )
        axis = go.Scatter(x=[0, self.span], y=[0, 0], mode='lines', line_color='black')
        st.plotly_chart(go.Figure(data=[line, axis], layout=layout), use_container_width=True)

    # ---------- Deflection ----------
    def _calc_deflection(self, M, EI, delx, theta_0, v_0, supportIndexA):
        """Forward integrate from support A using trapezoidal rule."""
        theta_im1 = theta_0
        v_im1 = v_0
        Rotation = np.zeros(len(self.X))
        Deflection = np.zeros(len(self.X))
        Rotation[supportIndexA] = theta_im1
        Deflection[supportIndexA] = v_im1

        for i, _ in enumerate(M[supportIndexA:]):
            ind = i + supportIndexA
            if i > 0:
                M_im1 = M[ind - 1]
                M_i = M[ind]
                M_avg = 0.5 * (M_i + M_im1)
                theta_i = theta_im1 + (M_avg / EI) * delx
                v_i = v_im1 + 0.5 * (theta_i + theta_im1) * delx
                Rotation[ind] = theta_i
                Deflection[ind] = v_i
                theta_im1 = theta_i
                v_im1 = v_i

        return Rotation, Deflection

    def _zero_crossing(self, Deflection, guessStep, initRot, initDef, M, EI, delx, supportIndexB, supportIndexA):
        """Find init rotation that crosses zero deflection at support B."""
        if Deflection[supportIndexB] > 0:
            error_is_positive = True
            while error_is_positive:
                initRot = initRot + guessStep
                _, Deflection = self._calc_deflection(M, EI, delx, initRot, initDef, supportIndexA)
                if Deflection[supportIndexB] <= 0:
                    error_is_positive = False
                    solved = initRot
        else:
            error_is_positive = False
            while not error_is_positive:
                initRot = initRot + guessStep
                _, Deflection = self._calc_deflection(M, EI, delx, initRot, initDef, supportIndexA)
                if Deflection[supportIndexB] >= 0:
                    error_is_positive = True
                    solved = initRot
        return solved

    def _compute_and_plot_deflection(self):
        deltaRot = 5e-4
        initRot = -2.1e-4
        M = -(sum(self.bendingMoment) if len(self.bendingMoment) else np.zeros_like(self.X))
        delx = self.X[1] - self.X[0] if len(self.X) > 1 else 0.0
        EI = self.E * self.I
        initDef = 0.0
        supportIndexA = int(np.where(self.X == self.A)[0].item()) if len(self.X) else 0
        supportIndexB = int(np.where(self.X == self.B)[0].item()) if len(self.X) else 0

        # Test direction for zero crossing
        testDef = np.zeros(3)
        for i, r in enumerate([initRot - deltaRot, initRot, initRot + deltaRot]):
            _, Defl = self._calc_deflection(M, EI, delx, r, initDef, supportIndexA)
            testDef[i] = Defl[supportIndexB]

        if abs(testDef[0]) < abs(testDef[1]):
            solvedInitRotation = self._zero_crossing(Defl, -deltaRot, initRot, initDef, M, EI, delx, supportIndexB, supportIndexA)
        elif abs(testDef[2]) < abs(testDef[1]):
            solvedInitRotation = self._zero_crossing(Defl, +deltaRot, initRot, initDef, M, EI, delx, supportIndexB, supportIndexA)
        else:
            solvedInitRotation = initRot

        Rotation, Deflection = self._calc_deflection(M, EI, delx, solvedInitRotation, initDef, supportIndexA)

        # Left overhang: integrate in reverse from A to x=0
        if self.A != 0:
            theta_im1 = -solvedInitRotation  # rotation continuity
            v_im1 = 0.0                      # deflection at A
            for i in np.arange(supportIndexA - 1, -1, -1):
                M_im1 = M[i + 1]
                M_i = M[i]
                M_avg = 0.5 * (M_i + M_im1)
                theta_i = theta_im1 + (M_avg / EI) * delx
                v_i = v_im1 + 0.5 * (theta_i + theta_im1) * delx
                Rotation[i] = theta_i
                Deflection[i] = v_i
                theta_im1 = theta_i
                v_im1 = v_i

        # Plot
        layout = go.Layout(
            title={'text': "Deflection", 'y': 0.85, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'},
            titlefont=dict(size=15),
            yaxis=dict(title='Deflection'),
            xaxis=dict(title='Distance (m)', range=[-1, self.span + 1]),
            showlegend=False,
        )
        line = go.Scatter(x=self.X, y=Deflection, mode='lines', name='Deflection',
                          line_color='orange', fillcolor='rgba(255, 255, 0, 0.1)')
        axis = go.Scatter(x=[0, self.span], y=[0, 0], mode='lines', line_color='black')
        st.plotly_chart(go.Figure(data=[line, axis], layout=layout), use_container_width=True)

        # ==== Extra results (overhang) ====
        S_total = sum(self.shearForce) if len(self.shearForce) else np.zeros_like(self.X)
        M_total = -(sum(self.bendingMoment) if len(self.bendingMoment) else np.zeros_like(self.X))

        i_max_down = int(np.argmin(Deflection))
        i_max_up   = int(np.argmax(Deflection))

        st.write("Maximum vertical deflection (down):", round(abs(Deflection[i_max_down]), 5), "m at x =", round(self.X[i_max_down], 3), "m")
        if Deflection[i_max_up] > 0:
            st.write("Maximum vertical deflection (up):", round(Deflection[i_max_up], 5), "m at x =", round(self.X[i_max_up], 3), "m")

        # deflection at notable locations
        def _interp_at(xq):
            j = np.searchsorted(self.X, xq)
            if j == 0: return float(Deflection[0])
            if j >= len(self.X): return float(Deflection[-1])
            x1, x2 = self.X[j-1], self.X[j]
            y1, y2 = Deflection[j-1], Deflection[j]
            return float(y1 + (y2 - y1) * (xq - x1) / (x2 - x1))

        st.write("Deflection at left end (x=0):", round(_interp_at(0.0), 5), "m")
        st.write("Deflection at right end (x=L):", round(_interp_at(self.span), 5), "m")
        st.write("Deflection at support A (x=A):", round(_interp_at(self.A), 5), "m")
        st.write("Deflection at support B (x=B):", round(_interp_at(self.B), 5), "m")

        # slopes at supports
        supportIndexA = int(np.where(self.X == self.A)[0].item()) if len(self.X) else 0
        supportIndexB = int(np.where(self.X == self.B)[0].item()) if len(self.X) else 0
        st.write("Rotation at A (slope):", round(float(Rotation[supportIndexA]), 7), "rad")
        st.write("Rotation at B (slope):", round(float(Rotation[supportIndexB]), 7), "rad")

        # shear/moment extrema
        iV = int(np.argmax(np.abs(S_total)))
        iM = int(np.argmax(np.abs(M_total)))
        st.write("Max |Shear|:", round(float(S_total[iV]), 4), "kN at x =", round(self.X[iV], 3), "m")
        st.write("Max |Moment|:", round(float(M_total[iM]), 4), "kNm at x =", round(self.X[iM], 3), "m")

        # contraflexure points
        # sign_changes = np.where(np.sign(M_total[:-1]) * np.sign(M_total[1:]) < 0)[0]
        # if len(sign_changes):
        #     xs = []
        #     for k in sign_changes:
        #         x1, x2 = self.X[k], self.X[k+1]
        #         y1, y2 = M_total[k], M_total[k+1]
        #         xs.append(float(x1 - y1 * (x2 - x1) / (y2 - y1)))
        #     st.write("Points of contraflexure (x where M=0):", [round(x, 3) for x in xs])

        # enforce exact zeros at supports (pins) to reduce jitter
        supportIndexA = int(np.where(self.X == self.A)[0].item()) if len(self.X) else 0
        supportIndexB = int(np.where(self.X == self.B)[0].item()) if len(self.X) else 0
        M_total[supportIndexA] = 0.0
        M_total[supportIndexB] = 0.0

        # contraflexure_xs = find_contraflexure_points(self.X, M_total, zero_tol=1e-2, min_gap=0.1)
        contraflexure_xs = find_contraflexure_points(self.X, M_total, zero_tol=1e-6, min_gap=0.01,
                                             prefer="bracket", min_amp=0.5)

        if contraflexure_xs:
            st.write("Points of contraflexure (M=0):", [round(x, 3) for x in contraflexure_xs])


        # max_deflection = abs(float(np.min(Deflection))) if len(Deflection) else 0.0
        # st.write("Maximum vertical deflection in downward direction:", round(max_deflection, 5), " m")

    # ---------- Public entry ----------
    def render(self):
        self.render_inputs()
        if st.button("Calculate"):
            if not self._parse_inputs():
                return
            self._accumulate_reactions()
            self._assemble_shear_moment()
            self._plot_shear()
            self._plot_moment()
            self._compute_and_plot_deflection()
