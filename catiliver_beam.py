# cantilever.py
import ast
import numpy as np
import plotly.graph_objects as go
import streamlit as st


class CantileverBeamApp:
    def __init__(self):
        # UI state (filled in by render_inputs)
        self.span = 0.0
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

        # Computation grids/records
        self.delta = 5e-5
        self.X = None
        self.reactions = np.array([0.0, 0.0])
        self.shearForce = None
        self.bendingMoment = None

        # Per-load bookkeeping
        self.PL_record = None
        self.PM_record = None
        self.UDL_record = None
        self.LDL_record = None
        self.TDL_record = None

    # ---------- UI helpers ----------
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
        st.title("Cantilever Beam Analysis Tool")

        st.markdown("#### Input Parameters")
        self.span = st.number_input("📏 Span of the beam (in meters):", min_value=0.0, format="%.2f")
        self.B = self.span
        st.markdown("*Units: meters (m) for lengths, kilonewtons (kN) for forces*")

        st.markdown("---")
        st.markdown("#### 🔹 Load Inputs")

        self.pointloads_txt = st.text_input(
            "📌 Point Loads",
            help="Enter point loads as: [[location₁, magnitude₁], [location₂, magnitude₂], ...]\n"
                 "➡️ Example: [[2, -10], [4.5, -15]] means 10kN downward load at 2m and 15kN downward load at 4.5m.",
            value=self.pointloads_txt,
        )
        self.pointmoments_txt = st.text_input(
            "♻️ Point Moments",
            help="Enter point moments as: [[location₁, magnitude₁], [location₂, magnitude₂], ...]\n"
                 "➡️ Example: [[1.5, 5], [3, -8]] means 5kNm clockwise at 1.5m and 8kNm counterclockwise at 3m.",
            value=self.pointmoments_txt,
        )
        self.distributedloads_txt = st.text_input(
            "🌫️ Uniform Distributed Loads (UDL)",
            help="Enter UDLs as: [[start₁, end₁, magnitude₁], [start₂, end₂, magnitude₂], ...]\n"
                 "➡️ Example: [[0, 2, -5], [3, 5, -8]] applies 5kN/m downward load from 0–2m and 8kN/m downward load from 3–5m.",
            value=self.distributedloads_txt,
        )
        self.linearloads_txt = st.text_input(
            "📈 Linear Varying Loads (Triangular)",
            help="Enter triangular loads as: [[start₁, end₁, start_mag₁, end_mag₁], ...]\n"
                 "➡️ Example: [[1, 3, 0, -10]] means load increases from 0kN/m to 10kN/m downward between 1m and 3m.",
            value=self.linearloads_txt,
        )
        self.trapezoidalloads_txt = st.text_input(
            "🧱 Trapezoidal Loads",
            help="Enter trapezoidal loads as: [[start₁, end₁, start_mag₁, end_mag₁], ...]\n"
                 "➡️ Example: [[2, 5, -4, -8]] means load varies from 4kN/m to 8kN/m downward between 2m and 5m.",
            value=self.trapezoidalloads_txt,
        )

        st.markdown("---")
        st.markdown("#### ⚙️ Material & Beam Properties")

        self.E = st.number_input(
            "🧪 Young's Modulus (E) in Pascals (N/m²):",
            help="Elastic property of material. Example for steel: 200000000 (200 GPa).",
            min_value=0, max_value=1_000_000_000, value=200_000_000, step=1_000_000,
        )
        self.I = st.number_input(
            "📏 Moment of Inertia (I) in m⁴:",
            help="Geometric property of the beam cross-section.\n➡️ Example: 0.0020833 for a rectangular section (b=0.3m, h=0.5m).",
            min_value=0.0, max_value=1.0, value=0.0020833, step=0.0000001, format="%.7f",
        )

    def _parse_inputs(self) -> bool:
        """Parse text → arrays; set grid and allocate containers."""
        self.pointloads = np.array(self._convert_to_list_of_lists(self.pointloads_txt))
        self.pointmoments = np.array(self._convert_to_list_of_lists(self.pointmoments_txt))
        self.distributedloads = np.array(self._convert_to_list_of_lists(self.distributedloads_txt))
        self.linearloads = np.array(self._convert_to_list_of_lists(self.linearloads_txt))
        self.trapezoidalloads = np.array(self._convert_to_list_of_lists(self.trapezoidalloads_txt))

        if any(x is None for x in [self.pointloads, self.pointmoments, self.distributedloads,
                                   self.linearloads, self.trapezoidalloads]):
            return False

        self.X = np.arange(0, self.span + self.delta, self.delta)
        self.shearForce = np.empty([0, len(self.X)])
        self.bendingMoment = np.empty([0, len(self.X)])
        return True

    # ---------- Reactions ----------
    def _reactions_PL(self, n):
        xp = self.pointloads[n, 0]
        fy = self.pointloads[n, 1]
        Vb = -fy
        return Vb

    def _reactions_PM(self, n):
        # Point moment does not add vertical reaction in cantilever (only end moment)
        return 0.0

    def _reactions_UDL(self, n):
        xStart, xEnd, fy = self.distributedloads[n, 0], self.distributedloads[n, 1], self.distributedloads[n, 2]
        fy_Res = fy * (xEnd - xStart)
        Vb = -fy_Res
        return Vb

    def _reactions_LDL(self, n):
        xStart, xEnd, fy_start, fy_end = self.linearloads[n]
        if abs(fy_start) > 0:
            fy_Res = 0.5 * fy_start * (xEnd - xStart)
        else:
            fy_Res = 0.5 * fy_end * (xEnd - xStart)
        Vb = -fy_Res
        return Vb

    def _reactions_TDL(self, n):
        xStart, xEnd, fy_start, fy_end = self.trapezoidalloads[n]
        fy_Res = 0.5 * (xEnd - xStart) * (fy_start + fy_end)
        Vb = -fy_Res
        return Vb

    def _accumulate_reactions(self):
        nPL = len(self.pointloads[0])
        nPM = len(self.pointmoments[0])
        nUDL = len(self.distributedloads[0])
        nLDL = len(self.linearloads[0])
        nTDL = len(self.trapezoidalloads[0])

        self.PL_record = np.empty([0, 1])
        self.PM_record = np.empty([0, 1])
        self.UDL_record = np.empty([0, 1])
        self.LDL_record = np.empty([0, 1])
        self.TDL_record = np.empty([0, 1])

        # point loads
        if nPL > 0:
            for n, _ in enumerate(self.pointloads):
                vb = self._reactions_PL(n)
                self.PL_record = np.append(self.PL_record, [[vb]], axis=0)
                self.reactions[0] += vb

        # point moments
        if nPM > 0:
            for n, _ in enumerate(self.pointmoments):
                vb = self._reactions_PM(n)
                self.PM_record = np.append(self.PM_record, [[vb]], axis=0)
                self.reactions[0] += vb

        # UDL
        if nUDL > 0:
            for n, _ in enumerate(self.distributedloads):
                vb = self._reactions_UDL(n)
                self.UDL_record = np.append(self.UDL_record, [[vb]], axis=0)
                self.reactions[0] += vb

        # LDL
        if nLDL > 0:
            for n, _ in enumerate(self.linearloads):
                vb = self._reactions_LDL(n)
                self.LDL_record = np.append(self.LDL_record, [[vb]], axis=0)
                self.reactions[0] += vb

        # TDL
        if nTDL > 0:
            for n, _ in enumerate(self.trapezoidalloads):
                vb = self._reactions_TDL(n)
                self.TDL_record = np.append(self.TDL_record, [[vb]], axis=0)
                self.reactions[0] += vb

        st.write(f'The vertical reaction at fixed end is {round(self.reactions[0], 2)} kN')

    # ---------- Shear/Moment contributions ----------
    def _sm_PL(self, n):
        xp, fy = self.pointloads[n, 0], self.pointloads[n, 1]
        Vb = self.PL_record[n, 0]
        Shear = np.zeros(len(self.X))
        Moment = np.zeros(len(self.X))
        for i, x in enumerate(self.X):
            shear = 0.0
            moment = 0.0
            if x >= self.B:
                shear += Vb
                moment -= Vb * (x - self.B)
            if x > xp:
                shear += fy
                moment -= fy * (x - xp)
            Shear[i] = shear
            Moment[i] = moment
        return Shear, Moment

    def _sm_PM(self, n):
        xm, m = self.pointmoments[n, 0], self.pointmoments[n, 1]
        Vb = self.PM_record[n, 0]
        Shear = np.zeros(len(self.X))
        Moment = np.zeros(len(self.X))
        for i, x in enumerate(self.X):
            shear = 0.0
            moment = 0.0
            if x >= self.B:
                shear += Vb
                moment -= Vb * (x - self.B)
            if x > xm:
                moment -= m
            Shear[i] = shear
            Moment[i] = moment
        return Shear, Moment

    def _sm_UDL(self, n):
        xStart, xEnd, fy = self.distributedloads[n]
        Vb = self.UDL_record[n, 0]
        Shear = np.zeros(len(self.X))
        Moment = np.zeros(len(self.X))
        for i, x in enumerate(self.X):
            shear = 0.0
            moment = 0.0
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
            Shear[i] = shear
            Moment[i] = moment
        return Shear, Moment

    def _sm_LDL(self, n):
        xStart, xEnd, fy_start, fy_end = self.linearloads[n]
        Vb = self.LDL_record[n, 0]
        Shear = np.zeros(len(self.X))
        Moment = np.zeros(len(self.X))
        for i, x in enumerate(self.X):
            shear = 0.0
            moment = 0.0
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

            Shear[i] = shear
            Moment[i] = moment
        return Shear, Moment

    def _sm_TDL(self, n):
        xStart, xEnd, fy_start, fy_end = self.trapezoidalloads[n]
        Vb = self.TDL_record[n, 0]
        Shear = np.zeros(len(self.X))
        Moment = np.zeros(len(self.X))
        for i, x in enumerate(self.X):
            shear = 0.0
            moment = 0.0
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

            Shear[i] = shear
            Moment[i] = moment
        return Shear, Moment

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
            yaxis=dict(title='Bending Moment (kNm)'),
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

    def _compute_and_plot_deflection(self):
        # Reverse integrate from fixed end B to free end
        M = -(sum(self.bendingMoment) if len(self.bendingMoment) else np.zeros_like(self.X))
        delx = self.X[1] - self.X[0] if len(self.X) > 1 else 0.0
        EI = self.E * self.I
        supportIndexB = int(np.where(self.X == self.B)[0].item()) if EI != 0 and delx > 0 else 0

        theta_im1 = 0.0
        v_im1 = 0.0
        Rotation = np.zeros(len(self.X))
        Deflection = np.zeros(len(self.X))
        Rotation[supportIndexB] = theta_im1
        Deflection[supportIndexB] = v_im1

        for i in np.arange(supportIndexB - 1, -1, -1):
            M_im1 = M[i + 1]
            M_i = M[i]
            M_avg = 0.5 * (M_i + M_im1)
            theta_i = theta_im1 + (M_avg / EI) * delx if EI != 0 else 0.0
            v_i = v_im1 + 0.5 * (theta_i + theta_im1) * delx
            Rotation[i] = theta_i
            Deflection[i] = v_i
            theta_im1 = theta_i
            v_im1 = v_i

        layout = go.Layout(
            title={'text': "Deflection", 'y': 0.85, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'},
            titlefont=dict(size=15),
            yaxis=dict(title='Deflection'),
            xaxis=dict(title='Distance (m)', range=[-1, self.span + 1]),
            showlegend=False,
        )
        line = go.Scatter(x=self.X, y=Deflection, mode='lines', name='Deflection', line_color='orange',
                          fillcolor='rgba(255, 255, 0, 0.1)')
        axis = go.Scatter(x=[0, self.span], y=[0, 0], mode='lines', line_color='black')
        st.plotly_chart(go.Figure(data=[line, axis], layout=layout), use_container_width=True)

        # ==== Extra results (cantilever) ====
        S_total = sum(self.shearForce) if len(self.shearForce) else np.zeros_like(self.X)
        M_total = -(sum(self.bendingMoment) if len(self.bendingMoment) else np.zeros_like(self.X))

        i_max_down = int(np.argmin(Deflection))
        i_max_up   = int(np.argmax(Deflection))

        st.write("Maximum vertical deflection (down):", round(abs(Deflection[i_max_down]), 5), "m at x =", round(self.X[i_max_down], 3), "m")
        if Deflection[i_max_up] > 0:
            st.write("Maximum vertical deflection (up):", round(Deflection[i_max_up], 5), "m at x =", round(self.X[i_max_up], 3), "m")

        # deflection at key points
        def _interp_at(xq):
            j = np.searchsorted(self.X, xq)
            if j == 0: return float(Deflection[0])
            if j >= len(self.X): return float(Deflection[-1])
            x1, x2 = self.X[j-1], self.X[j]
            y1, y2 = Deflection[j-1], Deflection[j]
            return float(y1 + (y2 - y1) * (xq - x1) / (x2 - x1))

        mid = self.span * 0.5
        st.write("Deflection at midspan (x=L/2):", round(_interp_at(mid), 5), "m")
        st.write("Deflection at free end (x=0):", round(float(Deflection[0]), 5), "m")
        st.write("Rotation at free end (slope, x=0):", round(float(Rotation[0]), 7), "rad")

        # shear/moment extrema
        iV = int(np.argmax(np.abs(S_total)))
        iM = int(np.argmax(np.abs(M_total)))
        st.write("Max |Shear|:", round(float(S_total[iV]), 4), "kN at x =", round(self.X[iV], 3), "m")
        st.write("Max |Moment|:", round(float(M_total[iM]), 4), "kNm at x =", round(self.X[iM], 3), "m")

        # points of contraflexure (M = 0)
        sign_changes = np.where(np.sign(M_total[:-1]) * np.sign(M_total[1:]) < 0)[0]
        if len(sign_changes):
            xs = []
            for k in sign_changes:
                x1, x2 = self.X[k], self.X[k+1]
                y1, y2 = M_total[k], M_total[k+1]
                xs.append(float(x1 - y1 * (x2 - x1) / (y2 - y1)))
            st.write("Points of contraflexure (x where M=0):", [round(x, 3) for x in xs])


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
