# simply_supported.py
import ast
import numpy as np
import plotly.graph_objects as go
import streamlit as st


class SimplySupportedBeamApp:
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
        self.PL_record = None  # (Va, Vb) per point load
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
        st.title("Simply Supported Beam Analysis Tool")

        st.markdown("#### Input Parameters")
        self.span = st.number_input("📏 Span of the beam (in meters):", min_value=0.0, format="%.2f")
        self.A = 0.0
        self.B = self.span
        st.markdown("*Units: meters (m) for lengths, kilonewtons (kN) for forces*")

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
        xm = self.pointmoments[n, 0]
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
            if x >= self.B:
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
        theta = theta_0
        v = v_0
        Rotation = np.zeros(len(self.X))
        Deflection = np.zeros(len(self.X))
        Rotation[supportIndexA] = theta
        Deflection[supportIndexA] = v

        for i in range(supportIndexA + 1, len(M)):
            M_avg = 0.5 * (M[i] + M[i - 1])
            theta += (M_avg / EI) * delx
            v += 0.5 * (theta + Rotation[i - 1]) * delx
            Rotation[i] = theta
            Deflection[i] = v

        return Rotation, Deflection

    def _zero_crossing(self, M, EI, delx, init_rot, init_def, step, supportIndexB, supportIndexA):
        """Adjust initial rotation to get near-zero deflection at support B."""
        prev_sign = np.sign(self._calc_deflection(M, EI, delx, init_rot, init_def, supportIndexA)[1][supportIndexB])
        while True:
            init_rot += step
            _, defl = self._calc_deflection(M, EI, delx, init_rot, init_def, supportIndexA)
            current = defl[supportIndexB]
            current_sign = np.sign(current)
            if current_sign != prev_sign:
                return init_rot
            prev_sign = current_sign

    def _compute_and_plot_deflection(self):
        deltaRot = 5e-5
        initRot = -2.1e-4
        M = -(sum(self.bendingMoment) if len(self.bendingMoment) else np.zeros_like(self.X))
        delx = self.X[1] - self.X[0] if len(self.X) > 1 else 0.0
        EI = self.E * self.I
        supportIndexA = int(np.where(self.X == self.A)[0].item()) if len(self.X) else 0
        supportIndexB = int(np.where(self.X == self.B)[0].item()) if len(self.X) else 0
        initDef = 0.0

        # pick direction based on reduced deflection
        test_rotations = [initRot - deltaRot, initRot, initRot + deltaRot]
        test_deflections = [
            self._calc_deflection(M, EI, delx, r, initDef, supportIndexA)[1][supportIndexB]
            for r in test_rotations
        ]

        if abs(test_deflections[0]) < abs(test_deflections[1]):
            solved_init_rotation = self._zero_crossing(M, EI, delx, initRot, initDef, -deltaRot, supportIndexB, supportIndexA)
        elif abs(test_deflections[2]) < abs(test_deflections[1]):
            solved_init_rotation = self._zero_crossing(M, EI, delx, initRot, initDef, +deltaRot, supportIndexB, supportIndexA)
        else:
            solved_init_rotation = initRot

        Rotation, Deflection = self._calc_deflection(M, EI, delx, solved_init_rotation, initDef, supportIndexA)

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

        # ==== Extra results (simply supported) ====
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

        for label, xq in [("L/4", self.span*0.25), ("L/2", self.span*0.5), ("3L/4", self.span*0.75)]:
            st.write(f"Deflection at {label}:", round(_interp_at(xq), 5), "m")

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
