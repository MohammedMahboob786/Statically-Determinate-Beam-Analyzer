import streamlit as st
import numpy as np
import plotly.graph_objects as go
import ast
from PIL import Image

# Set page config
st.set_page_config(page_title="Beam Analysis Tool", page_icon="🔧")


# Function for Cantilever Beam Analysis
def cantilever_beam():
    st.title("Cantiliver Beam Analysis Tool")

    # Input parameters
    span = st.number_input("Enter the span of the beam in meters:")
    B = span   
    st.subheader("Ensure units are in 'm' and 'kN'.")
    st.subheader("Point Loads")
    pointloads = st.text_input("Format for Point Loads: [[location 1,Magnitude 1], [location 2,Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Point Moments")
    pointmoments = st.text_input("Format for Point Moments: [[location 1,Magnitude 1], [location 2,Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Distributed Loads")
    distributedloads = st.text_input("Format for Distributed Loads: [[start location 1, end location 1,Magnitude 1], [start location 2, end location 2,Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Linear Loads")
    linearloads = st.text_input("Format for Linear Loads: [[start location 1, end location 1, start Magnitude 1, end Magnitude 1], [start location 2, end location 2, start Magnitude 2, end Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Trapezoidal Loads")
    trapezoidalloads = st.text_input("Format for Trapezoidal Loads: [[start location 1, end location 1, start Magnitude 1, end Magnitude 1], [start location 2, end location 2, start Magnitude 2, end Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Young's Modulus")
    E = st.number_input("Young's Modulus of the material in Pascals (N/m²):", min_value=0, max_value=1000000000, value=200000000)
    
    st.subheader("Moment of Inertia")
    I = st.number_input("Moment of Inertia of the beam in meters (m⁴):", min_value=0.0, max_value=1000000000.0, value=0.0020833, step=0.0000001, format="%.7f")


    # Convert input text to numpy arrays
    def convert_to_list_of_lists(input_str):
        try:
            # Use ast.literal_eval to safely evaluate the string
            result = ast.literal_eval(input_str)
            # Check if the result is a list of lists
            if isinstance(result, list) and all(isinstance(i, list) for i in result):
                return result
            else:
                st.error("Input is not a valid list of lists.")
                return None
        except (ValueError, SyntaxError):
            st.error("Invalid input format.")
            return None
        
    if st.button("Calculate"):
        pointloads = np.array(convert_to_list_of_lists(pointloads))
        pointmoments = np.array(convert_to_list_of_lists(pointmoments))
        distributedloads = np.array(convert_to_list_of_lists(distributedloads))
        linearloads = np.array(convert_to_list_of_lists(linearloads))
        trapezoidalloads = np.array(convert_to_list_of_lists(trapezoidalloads))
        
        delta = 0.00005
        X = np.arange(0, span+delta, delta)
        nPL = len(pointloads[0])
        nPM = len(pointmoments[0])
        nUDL = len(distributedloads[0])
        nLDL = len(linearloads[0])
        nTDL = len(trapezoidalloads[0])
            
        reactions = np.array([0.0, 0.0])
        shearForce = np.empty([0, len(X)])
        bendingMoment = np.empty([0, len(X)])
        
        
        def reactions_PL(n):
            xp = pointloads[n, 0]
            fy = pointloads[n, 1]

            Vb = -fy
            return Vb

        def reactions_PM(n):
            xm = pointmoments[n, 0]
            m = pointmoments[n, 1]
            
            Vb = 0
            return Vb

        def reactions_UDL(n):
            xStart = distributedloads[n,0]
            xEnd = distributedloads[n,1]
            fy = distributedloads[n,2]
            
            fy_Res = fy*(xEnd - xStart)
            x_Res = xStart + 0.5*(xEnd - xStart)
            
            Vb = -fy_Res
            return Vb
        
        def reactions_LDL(n):
            xStart = linearloads[n,0]
            xEnd = linearloads[n,1]
            fy_start = linearloads[n,2]
            fy_end = linearloads[n,3]
            
            #Determine location and magnitude of force resultant
            if abs(fy_start)>0:
                fy_Res = 0.5*fy_start*(xEnd-xStart)
            # x_Res = xStart + (1/3)*(xEnd-xStart)
            else:
                fy_Res = 0.5*fy_end*(xEnd-xStart)
            # x_Res = xStart + (2/3)*(xEnd-xStart)
            
            Vb = -fy_Res
            return Vb
        
        def reactions_TDL(n):
            xStart = trapezoidalloads[n,0]
            xEnd = trapezoidalloads[n,1]
            fy_start = trapezoidalloads[n,2]
            fy_end = trapezoidalloads[n,3]
            
            #determine loacation and magnitude of force resultant
            if abs(fy_start) > abs(fy_end):
                fy_Res = 0.5*(xEnd-xStart)*(fy_start+fy_end)
            # x_Res = xStart + ((1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end)))
            else:
                fy_Res = 0.5*(xEnd-xStart)*(fy_start+fy_end)
            # x_Res = xStart + ((1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end)))
            
            Vb = -fy_Res  
            return Vb

        PL_record = np.empty([0,1])
        if(nPL>0):
            for n, p in enumerate(pointloads):
                vb =reactions_PL(n)
                PL_record = np.append(PL_record, [np.array([vb])], axis =0)

                reactions[0] = reactions[0] + vb

        PM_record = np.empty([0,1])
        if(nPM>0):
            for n, p in enumerate(pointmoments):
                vb = reactions_PM(n)
                PM_record = np.append(PM_record, [np.array([vb])], axis =0)
                
                reactions[0] = reactions[0] + vb
                
        UDL_record = np.empty([0,1])
        if(nUDL>0):
            for n, p in enumerate(distributedloads):
                vb = reactions_UDL(n)
                UDL_record = np.append(UDL_record, [np.array([vb])], axis =0)
                
                reactions[0] = reactions[0] + vb
        
        LDL_record = np.empty([0,1])
        if(nLDL>0):
            for n, p in enumerate(linearloads):
                vb = reactions_LDL(n)
                LDL_record = np.append(LDL_record, [np.array([vb])], axis =0)
                
                reactions[0] = reactions[0] + vb
                
        TDL_record = np.empty([0,1])
        if(nTDL>0):
            for n, p in enumerate(trapezoidalloads):
                vb = reactions_TDL(n)
                TDL_record = np.append(TDL_record, [np.array([vb])], axis =0)
                
                reactions[0] = reactions[0] + vb        

        #Plotting and printing

        st.write('The vertical reaction at fixed end is {} kN'.format(round(reactions[0], 2)))        
                ##          SHEAR AND MOMENT CALCULATIONS
        #Define function to calculate shear forces and bending moments due to point loads

        #7
        def shear_moment_PL(n):
            xp = pointloads[n, 0]
            fy = pointloads[n, 1]
            Vb = PL_record[n, 0]
            
            #cycle through the structure and calculate the shearforce and bending moment at each point
            Shear = np.zeros(len(X))
            Moment  = np.zeros(len(X))

            for i, x in enumerate(X):
                shear = 0
                moment = 0
                
                if x >= B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                    
                if x > xp:
                    # calculate shear and moment due to reaction at point load
                    shear = shear + fy
                    moment = moment - fy*(x - xp)
                    
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
        #Define function to calculate shear forces and bending moments due to point moments
        def shear_moment_PM(n):
            xm = pointmoments[n, 0]
            m = pointmoments[n, 1]
            Vb = PM_record[n, 0]
            
            #cycle through the structure and calculate the shearforce and bending moment at each point
            Shear = np.zeros(len(X))
            Moment  = np.zeros(len(X))

            for i, x in enumerate(X):
                shear = 0
                moment = 0
                
                if x >= B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                
                if x > xm:
                    # calculate shear and moment due to reaction at point load
                    moment = moment - m
                    
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
        def shear_moment_UDL(n):
            xStart = distributedloads[n,0]
            xEnd = distributedloads[n,1]
            fy = distributedloads[n,2]
            Vb = UDL_record[n,0]
            
            #cycle through the structure and calculate the shearforce and bending moment at each point
            Shear = np.zeros(len(X))
            Moment  = np.zeros(len(X))

            for i, x in enumerate(X):
                shear = 0
                moment = 0
                
                if x >= B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)

                if x>xStart and x<=xEnd:
                    # calculate shear and moment due to reaction at point load
                    shear = shear + fy*(x-xStart)
                    moment = moment - fy*(x-xStart)*0.5*(x-xStart)
                elif (x>xEnd):
                    shear = shear + fy*(xEnd-xStart)
                    moment = moment - fy*(xEnd-xStart)*(x-xStart-0.5*(xEnd-xStart))
                    
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
                #Define function to calculate shear forces and bending moments due to LDLs
        def shear_moment_LDL(n):
            xStart = linearloads[n,0]
            xEnd = linearloads[n,1]
            fy_start = linearloads[n,2]
            fy_end = linearloads[n,3]
            Vb = LDL_record[n,0]
            
            
            Shear = np.zeros(len(X))
            Moment = np.zeros(len(X))
            
            for i, x in enumerate(X):
                shear = 0
                moment = 0
                
                if x >= B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                
                if x>xStart and x<=xEnd:
                    if abs(fy_start)>0:
                        x_base = x-xStart
                        f_cut = fy_start - x_base*(fy_start/(xEnd-xStart))
                        R1 = 0.5*x_base*(fy_start-f_cut)
                        R2 = x_base*f_cut
                        shear = shear + R1 + R2
                        moment = moment - R1*(2/3)*x_base - R2*(x_base/2)       
                    else:
                        x_base  = x-xStart
                        f_cut = fy_end*(x_base/(xEnd-xStart))
                        R = 0.5*x_base*f_cut
                        shear = shear + R
                        moment = moment -R*(x_base/3)
                        
                elif x>xEnd:
                    if abs(fy_start)>0:
                        R = 0.5*fy_start*(xEnd-xStart)
                        xr = xStart + (1/3)*(xEnd-xStart)
                        shear = shear + R
                        moment = moment - R*(x-xr)
                    else:
                        R=0.5*fy_end*(xEnd-xStart)
                        xr = xStart + (2/3)*(xEnd-xStart)
                        shear = shear + R
                        moment = moment - R*(x-xr)
                        
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
                #Define function to calculate shear forces and bending moments due to TDLs
        def shear_moment_TDL(n):
            xStart = trapezoidalloads[n,0]
            xEnd = trapezoidalloads[n,1]
            fy_start = trapezoidalloads[n,2]
            fy_end = trapezoidalloads[n,3]
            Vb = TDL_record[n,0]
            
            Shear = np.zeros(len(X))
            Moment = np.zeros(len(X))
            
            for i, x in enumerate(X):
                shear = 0
                moment = 0
                
                if x >= B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                
                if x>xStart and x<=xEnd:
                    if abs(fy_start)>abs(fy_end):
                        x_base = x-xStart
                        f_cut=fy_start-x_base*((fy_start-fy_end)/(xEnd-xStart))
                        
                        R1= 0.5*(x_base)*(fy_start-f_cut)
                        la_R1=(2/3)*(x_base)
                        
                        R2= x_base*f_cut
                        la_R2=0.5*(x_base)
                        
                        shear = shear + R1 + R2
                        moment = moment - R1*la_R1 - R2*la_R2
                        
                    else:
                        x_base = x-xStart
                        f_cut=fy_start+x_base*((fy_end-fy_start)/(xEnd-xStart))
                        
                        R1= 0.5*(x_base)*(f_cut-fy_start)
                        la_R1=(1/3)*(x_base)
                        
                        R2= (x_base)*(fy_start)
                        la_R2=0.5*(x_base)
                        
                        shear = shear + R1 + R2
                        moment = moment - R1*la_R1 - R2*la_R2
                    
                elif x>xEnd:
                    if abs(fy_start)>abs(fy_end):
                        R = 0.5*(xEnd-xStart)*(fy_start+fy_end)
                        xr = xStart + (1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end))
                        
                        shear = shear + R
                        moment = moment - R*(x-xr)
                    
                    else:
                        R=0.5*(xEnd-xStart)*(fy_end+fy_start)
                        xr=xStart + (1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end))
                        
                        shear = shear + R
                        moment = moment - R*(x-xr)
                
                Shear[i] = shear
                Moment[i] = moment

            return Shear, Moment
        
        #Cycle through all the point loads and determine shear and moment

        #8
        if(nPL>0):
            for n, p in enumerate(pointloads):
                Shear, Moment = shear_moment_PL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each point load
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each point load        
        if(nPM>0):
            for n, p in enumerate(pointmoments):
                Shear, Moment = shear_moment_PM(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each point moments
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each point moments
        if(nUDL>0):
            for n, p in enumerate(distributedloads):
                Shear, Moment = shear_moment_UDL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each point moments
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each point moments      
        if(nLDL>0):
            for n, p in enumerate(linearloads):
                Shear, Moment = shear_moment_LDL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each LDL
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each LDL  
        if(nTDL>0):
            for n, p in enumerate(trapezoidalloads):
                Shear, Moment = shear_moment_TDL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each TDL
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each TDL
        #Define the layout object
        layout = go.Layout(
            title={
                'text': "Shear Force Diagram",
                'y':0.85,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'},
            titlefont=dict(size=15),
            yaxis = dict(
                title='Shear Force (kN)'
            ),
            xaxis = dict(
                title='Distance (m)',       
                range=[-1, span+1]
            ),
            showlegend=False,        
        )

        #Define the shear force trace
        line = go.Scatter(
            x = X,
            y = sum(shearForce),
            mode='lines',
            name='Shear Force',
            fill='tonexty',
            line_color='green',
            fillcolor='rgba(0, 255, 0, 0.1)'
        )

        #Define a horizontal line to represent the structure
        axis = go.Scatter(
            x = [0, span],
            y = [0,0],
            mode='lines',
            line_color='black'
        )

        #Generate and view the figure
        fig = go.Figure(data=[line, axis], layout=layout)
        st.plotly_chart(fig)
      
      
        
        layout = go.Layout(
        title={
            'text': "Bending Moment Diagram",
            'y':0.85,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        titlefont=dict(size=15),
        yaxis = dict(
            title='Bending Moment (kNm)'
        #    autorange="reversed",   
        ),
        xaxis = dict(
            title='Distance (m)',       
            range=[-1, span+1]
        ),
        showlegend=False 
    )

    #Define the Bending Moment trace
        line = go.Scatter(
            x = X,
            y = -sum(bendingMoment),
            mode='lines',
            name='Bending Moment',
            fill='tonexty',
            line_color='red',
            fillcolor='rgba(255, 0, 0, 0.1)'
        )

        #Define a horizontal line to represent the structure
        axis = go.Scatter(
            x = [0, span],
            y = [0,0],
            mode='lines',
            line_color='black'
        )

        #Generate and view the figure
        fig = go.Figure(data=[line, axis], layout=layout)
        st.plotly_chart(fig)
        
        
        deltaRot = 0.000005
        initRot = -0.00021       

        M = -sum(bendingMoment)
        delx = X[1]-X[0]
        EI = E*I
        initDef = 0
        supportIndexB = np.where(X==B)[0].item()
            
        print("solve for deflection by integrating in reverse direction")
            
        theta_im1 = 0 #Rotation on other side of support A
        v_im1 = 0 #Vertical deflection at support A
        Rotation = np.zeros(len(X))
        Deflection = np.zeros(len(X))
        Rotation[supportIndexB] = theta_im1
        Deflection[supportIndexB] = v_im1

        #Generate a range of indices in reverse direction from support A to left endge of beam
        reverseRange = np.arange(supportIndexB-1,-1,-1) 

        #Loop through data and integrate (Trapezoidal rule) - REVERSE DIRECTION
        for i in reverseRange:                        
            M_im1 = M[i+1] #(300) - Assign previous value of M (reverse direction)
            M_i = M[i] #(299) - Assign current value of M (reverse direction)
            M_avg = 0.5*(M_i + M_im1)   
                
            theta_i = theta_im1 + (M_avg/EI)*delx #Integrate moment values to get rotations                     
            v_i = v_im1 + 0.5*(theta_i+theta_im1)*delx #Integrate rotation values to get displacements
                
            #Store data
            Rotation[i] = theta_i        
            Deflection[i] = v_i

            #Update values for next loop iteration
            theta_im1 = theta_i     
            v_im1 = v_i 

        #Define the shear force trace
        layout = go.Layout(
            title={
                'text': "Deflection",
                'y':0.85,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'},
            titlefont=dict(size=15),
            yaxis = dict(
                title='Deflection'
            ),
            xaxis = dict(
                title='Distance (m)',       
                range=[-1, span+1]
            ),
            showlegend=False,        
        )

        #Define the shear force trace
        line = go.Scatter(
            x = X,
            y = Deflection,
            mode='lines',
            name='Deflection',
            line_color='orange',
        # fill='tonexty',
            fillcolor='rgba(255, 255, 0, 0.1)'
        )

        #Define a horizontal line to represent the structure
        axis = go.Scatter(
            x = [0, span],
            y = [0,0],
            mode='lines',
            line_color='black'
        )

        #Generate and view the figure
        fig = go.Figure(data=[line,axis], layout=layout)
        st.plotly_chart(fig)   

        max_deflection=abs(min(Deflection))
        st.write("Maximum vertical deflection in downward direction:",round(max_deflection,5), " m")
    

# Function for Simply Supported Beam Analysis
def simply_supported_beam():
    st.title("Simply Supported Beam Analysis Tool")

    # Input parameters
    span = st.number_input("Enter the span of the beam in meters:")
    A = 0
    B = span
    st.subheader("Ensure units are in 'm' and 'kN'.")
    st.subheader("Point Loads")
    pointloads = st.text_input("Format for Point Loads: [[location 1,Magnitude 1], [location 2,Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Point Moments")
    pointmoments = st.text_input("Format for Point Moments: [[location 1,Magnitude 1], [location 2,Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Distributed Loads")
    distributedloads = st.text_input("Format for Distributed Loads: [[start location 1, end location 1,Magnitude 1], [start location 2, end location 2,Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Linear Loads")
    linearloads = st.text_input("Format for Linear Loads: [[start location 1, end location 1, start Magnitude 1, end Magnitude 1], [start location 2, end location 2, start Magnitude 2, end Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Trapezoidal Loads")
    trapezoidalloads = st.text_input("Format for Trapezoidal Loads: [[start location 1, end location 1, start Magnitude 1, end Magnitude 1], [start location 2, end location 2, start Magnitude 2, end Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Young's Modulus")
    E = st.number_input("Young's Modulus of the material in Pascals (N/m²):", min_value=0, max_value=1000000000, value=200000000)
    
    st.subheader("Moment of Inertia")
    I = st.number_input("Moment of Inertia of the beam in meters (m⁴):", min_value=0.0, max_value=1000000000.0, value=0.0020833, step=0.0000001, format="%.7f")


    # Convert input text to numpy arrays
    def convert_to_list_of_lists(input_str):
        try:
            # Use ast.literal_eval to safely evaluate the string
            result = ast.literal_eval(input_str)
            # Check if the result is a list of lists
            if isinstance(result, list) and all(isinstance(i, list) for i in result):
                return result
            else:
                st.error("Input is not a valid list of lists.")
                return None
        except (ValueError, SyntaxError):
            st.error("Invalid input format.")
            return None
        
    if st.button("Calculate"):
        pointloads = np.array(convert_to_list_of_lists(pointloads))
        pointmoments = np.array(convert_to_list_of_lists(pointmoments))
        distributedloads = np.array(convert_to_list_of_lists(distributedloads))
        linearloads = np.array(convert_to_list_of_lists(linearloads))
        trapezoidalloads = np.array(convert_to_list_of_lists(trapezoidalloads))
        
        delta = 0.00005
        X = np.arange(0, span+delta, delta)
        nPL = len(pointloads[0])
        nPM = len(pointmoments[0])
        nUDL = len(distributedloads[0])
        nLDL = len(linearloads[0])
        nTDL = len(trapezoidalloads[0])
            
        reactions = np.array([0.0, 0.0])
        shearForce = np.empty([0, len(X)])
        bendingMoment = np.empty([0, len(X)])
        
        
        def reactions_PL(n):
            xp = pointloads[n, 0]
            fy = pointloads[n, 1]

            la_p = A-xp
            mp = fy*la_p
            la_vb = B-A
            Vb = mp/la_vb
            Va = -fy-Vb

            return Va, Vb

        def reactions_PM(n):
            xm = pointmoments[n, 0]
            m = pointmoments[n, 1]
            la_vb = B-A
            
            Vb = m/la_vb
            Va = -Vb
            
            return Va, Vb

        def reactions_UDL(n):
            xStart = distributedloads[n,0]
            xEnd = distributedloads[n,1]
            fy = distributedloads[n,2]
            
            fy_Res = fy*(xEnd - xStart)
            x_Res = xStart + 0.5*(xEnd - xStart)
            
            la_p = A-x_Res
            mp = fy_Res*la_p
            la_vb = B-A
            
            Vb = mp/la_vb
            Va = -fy_Res -Vb
            
            return Va, Vb
        
        def reactions_LDL(n):
            xStart = linearloads[n,0]
            xEnd = linearloads[n,1]
            fy_start = linearloads[n,2]
            fy_end = linearloads[n,3]
            
            #Determine location and magnitude of force resultant
            if abs(fy_start)>0:
                fy_Res = 0.5*fy_start*(xEnd-xStart)
                x_Res = xStart + (1/3)*(xEnd-xStart)
            else:
                fy_Res = 0.5*fy_end*(xEnd-xStart)
                x_Res = xStart + (2/3)*(xEnd-xStart)
            
            la_p = A-x_Res
            mp = fy_Res*la_p
            la_vb = B-A
            
            Vb = mp/la_vb
            Va = -fy_Res -Vb
            
            return Va, Vb
        
        def reactions_TDL(n):
            xStart = trapezoidalloads[n,0]
            xEnd = trapezoidalloads[n,1]
            fy_start = trapezoidalloads[n,2]
            fy_end = trapezoidalloads[n,3]
            
            #determine loacation and magnitude of force resultant
            if abs(fy_start) > abs(fy_end):
                fy_Res = 0.5*(xEnd-xStart)*(fy_start+fy_end)
                x_Res = xStart + ((1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end)))
            else:
                fy_Res = 0.5*(xEnd-xStart)*(fy_start+fy_end)
                x_Res = xStart + ((1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end)))
                
            la_p = A-x_Res
            mp = fy_Res*la_p
            la_vb = B-A
            Vb = mp/la_vb
            Va = -fy_Res -Vb
            
            return Va, Vb

        PL_record = np.empty([0,2])
        if(nPL>0):
            for n, p in enumerate(pointloads):
                va, vb =reactions_PL(n)
                PL_record = np.append(PL_record, [np.array([va, vb])], axis =0)

                reactions[0] = reactions[0] + va
                reactions[1] = reactions[1] + vb

        PM_record = np.empty([0,2])
        if(nPM>0):
            for n, p in enumerate(pointmoments):
                va, vb = reactions_PM(n)
                PM_record = np.append(PM_record, [np.array([va, vb])], axis =0)
                
                reactions[0] = reactions[0] + va
                reactions[1] = reactions[1] + vb
                
        UDL_record = np.empty([0,2])
        if(nUDL>0):
            for n, p in enumerate(distributedloads):
                va, vb = reactions_UDL(n)
                UDL_record = np.append(UDL_record, [np.array([va, vb])], axis =0)
                
                reactions[0] = reactions[0] + va
                reactions[1] = reactions[1] + vb
        
        LDL_record = np.empty([0,2])
        if(nLDL>0):
            for n, p in enumerate(linearloads):
                va, vb = reactions_LDL(n)
                LDL_record = np.append(LDL_record, [np.array([va, vb])], axis =0)
                
                reactions[0] = reactions[0] + va
                reactions[1] = reactions[1] + vb
                
        TDL_record = np.empty([0,2])
        if(nTDL>0):
            for n, p in enumerate(trapezoidalloads):
                va, vb = reactions_TDL(n)
                TDL_record = np.append(TDL_record, [np.array([va, vb])], axis =0)
                
                reactions[0] = reactions[0] + va
                reactions[1] = reactions[1] + vb        

        #Plotting and printing

        st.write('The vertical reaction at A is {} kN'.format(round(reactions[0], 2)))
        st.write('The vertical reaction at B is {} kN'.format(round(reactions[1], 2)))
        
                ##          SHEAR AND MOMENT CALCULATIONS
        #Define function to calculate shear forces and bending moments due to point loads

        #7
        def shear_moment_PL(n):
            xp = pointloads[n, 0]
            fy = pointloads[n, 1]
            Va = PL_record[n, 0]
            Vb = PL_record[n, 1]

            #cycle through the structure and calculate the shearforce and bending moment at each point
            Shear = np.zeros(len(X))
            Moment  = np.zeros(len(X))

            for i, x in enumerate(X):
                shear = 0
                moment = 0

                if x > A:
                    #calculate shear and moment due to reaction at A
                    shear = shear + Va
                    moment = moment - Va*(x-A)
                
                if x >= B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                
                if x > xp:
                    # calculate shear and moment due to reaction at point load
                    shear = shear + fy
                    moment = moment - fy*(x - xp)
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
        #Define function to calculate shear forces and bending moments due to point moments
        def shear_moment_PM(n):
            xm = pointmoments[n, 0]
            m = pointmoments[n, 1]
            Va = PM_record[n, 0]
            Vb = PM_record[n, 1]
            
            #cycle through the structure and calculate the shearforce and bending moment at each point
            Shear = np.zeros(len(X))
            Moment  = np.zeros(len(X))

            for i, x in enumerate(X):
                shear = 0
                moment = 0

                if x > A:
                    #calculate shear and moment due to reaction at A
                    shear = shear + Va
                    moment = moment - Va*(x - A)
                    
                if x >= B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                
                if x > xm:
                    # calculate shear and moment due to reaction at point load
                    moment = moment - m
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
        def shear_moment_UDL(n):
            xStart = distributedloads[n,0]
            xEnd = distributedloads[n,1]
            fy = distributedloads[n,2]
            Va = UDL_record[n,0]
            Vb = UDL_record[n,1]
            
            #cycle through the structure and calculate the shearforce and bending moment at each point
            Shear = np.zeros(len(X))
            Moment  = np.zeros(len(X))

            for i, x in enumerate(X):
                shear = 0
                moment = 0

                if x > A:
                    #calculate shear and moment due to reaction at A
                    shear = shear + Va
                    moment = moment - Va*(x-A)
                    
                if x >= B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
            
                if x>xStart and x<=xEnd:
                    # calculate shear and moment due to reaction at point load
                    shear = shear + fy*(x-xStart)
                    moment = moment - fy*(x-xStart)*0.5*(x-xStart)
                elif (x>xEnd):
                    shear = shear + fy*(xEnd-xStart)
                    moment = moment - fy*(xEnd-xStart)*(x-xStart-0.5*(xEnd-xStart))
                    
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
                #Define function to calculate shear forces and bending moments due to LDLs
        def shear_moment_LDL(n):
            xStart = linearloads[n,0]
            xEnd = linearloads[n,1]
            fy_start = linearloads[n,2]
            fy_end = linearloads[n,3]
            Va = LDL_record[n,0]
            Vb = LDL_record[n,1]
            
            Shear = np.zeros(len(X))
            Moment = np.zeros(len(X))
            
            for i, x in enumerate(X):
                shear = 0
                moment = 0
                
                if x > A:
                    #calculate shear and moment due to reaction at A
                    shear = shear + Va
                    moment = moment - Va*(x-A)
                    
                if x >= B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                
                if x>xStart and x<=xEnd:
                    if abs(fy_start)>0:
                        x_base = x-xStart
                        f_cut = fy_start - x_base*(fy_start/(xEnd-xStart))
                        R1 = 0.5*x_base*(fy_start-f_cut)
                        R2 = x_base*f_cut
                        shear = shear + R1 + R2
                        moment = moment - R1*(2/3)*x_base - R2*(x_base/2)      
                    else:
                        x_base  = x-xStart
                        f_cut = fy_end*(x_base/(xEnd-xStart))
                        R = 0.5*x_base*f_cut
                        shear = shear + R
                        moment = moment -R*(x_base/3)
                        
                elif x>xEnd:
                    if abs(fy_start)>0:
                        R = 0.5*fy_start*(xEnd-xStart)
                        xr = xStart + (1/3)*(xEnd-xStart)
                        shear = shear + R
                        moment = moment - R*(x-xr)
                    else:
                        R=0.5*fy_end*(xEnd-xStart)
                        xr = xStart + (2/3)*(xEnd-xStart)
                        shear = shear + R
                        moment = moment - R*(x-xr)
                        
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
                #Define function to calculate shear forces and bending moments due to TDLs
        def shear_moment_TDL(n):
            xStart = trapezoidalloads[n,0]
            xEnd = trapezoidalloads[n,1]
            fy_start = trapezoidalloads[n,2]
            fy_end = trapezoidalloads[n,3]
            Va = TDL_record[n,0]
            Vb = TDL_record[n,1]
            
            Shear = np.zeros(len(X))
            Moment = np.zeros(len(X))
            
            for i, x in enumerate(X):
                shear = 0
                moment = 0
                
                if x > A:
                    #calculate shear and moment due to reaction at A
                    shear = shear + Va
                    moment = moment - Va*(x-A)
                if x >= B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                    
                if x>xStart and x<=xEnd:
                    if abs(fy_start)>abs(fy_end):
                        x_base = x-xStart
                        f_cut=fy_start-x_base*((fy_start-fy_end)/(xEnd-xStart))
                        
                        R1= 0.5*(x_base)*(fy_start-f_cut)
                        la_R1=(2/3)*(x_base)
                        
                        R2= x_base*f_cut
                        la_R2=0.5*(x_base)
                        
                        shear = shear + R1 + R2
                        moment = moment - R1*la_R1 - R2*la_R2
                        
                    else:
                        x_base = x-xStart
                        f_cut=fy_start+x_base*((fy_end-fy_start)/(xEnd-xStart))
                        
                        R1= 0.5*(x_base)*(f_cut-fy_start)
                        la_R1=(1/3)*(x_base)
                        
                        R2= (x_base)*(fy_start)
                        la_R2=0.5*(x_base)
                        
                        shear = shear + R1 + R2
                        moment = moment - R1*la_R1 - R2*la_R2
                    
                elif x>xEnd:
                    if abs(fy_start)>abs(fy_end):
                        R = 0.5*(xEnd-xStart)*(fy_start+fy_end)
                        xr = xStart + (1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end))
                        
                        shear = shear + R
                        moment = moment - R*(x-xr)
                    
                    else:
                        R=0.5*(xEnd-xStart)*(fy_end+fy_start)
                        xr=xStart + (1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end))
                        
                        shear = shear + R
                        moment = moment - R*(x-xr)
                
                Shear[i] = shear
                Moment[i] = moment

            return Shear, Moment
        
        #Cycle through all the point loads and determine shear and moment

        #8
        if(nPL>0):
            for n, p in enumerate(pointloads):
                Shear, Moment = shear_moment_PL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each point load
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each point load
        
        if(nPM>0):
            for n, p in enumerate(pointmoments):
                Shear, Moment = shear_moment_PM(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each point moments
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each point moments        

        if(nUDL>0):
            for n, p in enumerate(distributedloads):
                Shear, Moment = shear_moment_UDL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each UDL
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each UDL        

        if(nLDL>0):
            for n, p in enumerate(linearloads):
                Shear, Moment = shear_moment_LDL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each LDL
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each LDL        

        if(nTDL>0):
            for n, p in enumerate(trapezoidalloads):
                Shear, Moment = shear_moment_TDL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each TDL
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each TDL

        #Define the layout object
        layout = go.Layout(
            title={
                'text': "Shear Force Diagram",
                'y':0.85,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'},
            titlefont=dict(size=15),
            yaxis = dict(
                title='Shear Force (kN)'
            ),
            xaxis = dict(
                title='Distance (m)',       
                range=[-1, span+1]
            ),
            showlegend=False,        
        )

        #Define the shear force trace
        line = go.Scatter(
            x = X,
            y = sum(shearForce),
            mode='lines',
            name='Shear Force',
            fill='tonexty',
            line_color='green',
            fillcolor='rgba(0, 255, 0, 0.1)'
        )

        #Define a horizontal line to represent the structure
        axis = go.Scatter(
            x = [0, span],
            y = [0,0],
            mode='lines',
            line_color='black'
        )

        #Generate and view the figure
        fig = go.Figure(data=[line, axis], layout=layout)
        st.plotly_chart(fig)
      
      
        
        layout = go.Layout(
        title={
            'text': "Bending Moment Diagram",
            'y':0.85,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        titlefont=dict(size=15),
        yaxis = dict(
            title='Bending Moment (kNm)',
            autorange="reversed",   
        ),
        xaxis = dict(
            title='Distance (m)',       
            range=[-1, span+1]
        ),
        showlegend=False 
    )

        #Define the Bending Moment trace
        line = go.Scatter(
            x = X,
            y = -sum(bendingMoment),
            mode='lines',
            name='Bending Moment',
            fill='tonexty',
            line_color='red',
            fillcolor='rgba(255, 0, 0, 0.1)'
        )

        #Define a horizontal line to represent the structure
        axis = go.Scatter(
            x = [0, span],
            y = [0,0],
            mode='lines',
            line_color='black'
        )

        #Generate and view the figure
        fig = go.Figure(data=[line, axis], layout=layout)
        st.plotly_chart(fig)
        
        deltaRot = 0.000005
        initRot = -0.00021       

        M = -sum(bendingMoment)
        delx = X[1]-X[0]
        EI = E*I
        initDef = 0
        supportIndexA = np.where(X==A)[0].item()
        supportIndexB = np.where(X==B)[0].item()
            
        # Define a function to calculate deflection by forward integrating (from left support) using the Trapezoidal Rule

        def calcDeflection(M, EI, delx, theta_0, v_0):
                
            theta_im1 = theta_0
            v_im1 = v_0
                
            Rotation = np.zeros(len(X))
            Deflection = np.zeros(len(X))
            Rotation[supportIndexA] = theta_im1
            Deflection[supportIndexA] = v_im1
                
            #loop through the data and integrate (Trapezoidal Rule)
            for i, m in enumerate(M[supportIndexA::]):
                ind = i+supportIndexA
                if i>0:
                    M_im1 = M[ind-1]
                    M_i = M[ind]
                    M_avg = 0.5*(M_i+M_im1)
                        
                    theta_i = theta_im1 + (M_avg/EI)*delx
                    v_i = v_im1 + 0.5*(theta_i+theta_im1)*delx
                        
                    #store data
                    Rotation[ind] = theta_i
                    Deflection[ind] = v_i
                        
                    #update values for next loop iteration
                    theta_im1 = theta_i
                    v_im1 = v_i
                        
            return Rotation, Deflection
            
        def zeroCrossing(Deflection, guessStep, initRot, initDef):
            """
            Find the value of initial rotation that minimised deflection at right side 
            support by identifying where error crosses zero.
            """
            
            #If the deflection error is positive
            if Deflection[supportIndexB]>0:
                errorIsPositive = True #Set flag for error sign
                
                #Keep testing lower initial rotation values until error turns NEGATIVE
                while errorIsPositive:
                    initRot = initRot + guessStep
                    Rotation, Deflection = calcDeflection(M, EI, delx, initRot, initDef)
                    
                    #If error has turned NEGATIVE, switch the flag to allow loop to stop
                    if Deflection[supportIndexB]<=0:
                        errorIsPositive = False
                        solvedInitRotation = initRot #Save the 'solved' value that minimised the error
            
            #Else if deflection error is negative
            elif (Deflection[supportIndexB]<0):
                errorIsPositive = False #Set flag for error sign
                
                #Keep testing lower initial rotation values until error turns POSITIVE
                while not errorIsPositive:
                    initRot = initRot + guessStep
                    Rotation, Deflection = calcDeflection(M, EI, delx, initRot, initDef)
                    
                    #If error has turned POSITIVE, switch the flag to allow loop to stop
                    if Deflection[supportIndexB]>=0:
                        errorIsPositive = True
                        solvedInitRotation = initRot #Save the 'solved' value that minimised the error
            
            return solvedInitRotation
            
        #Test whether reducing or increasing initial rotation leads to reduction in disp error at other support
        testDef = np.zeros(3)
        for i, r in enumerate([initRot-deltaRot, initRot, initRot+deltaRot]):
            Rotation, Deflection = calcDeflection(M, EI, delx, r, initDef)
            testDef[i] = Deflection[supportIndexB]
                
        if(abs(testDef[0])<abs(testDef[1])):
            #Need to test in the negative rotation direction by reducing the initial rotation guess
            print('Need to test in the negative direction')    
            solvedInitRotation = zeroCrossing(Deflection, -deltaRot, initRot, initDef)            

        elif(abs(testDef[2])<abs(testDef[1])):
            #Need to test in the positive rotation direction by incresing the initial rotation guess 
            print('Need to test in the positive direction')    
            solvedInitRotation = zeroCrossing(Deflection, +deltaRot, initRot, initDef)      

        #Run the deflection calculation with the solved value of initial rotation
        Rotation, Deflection = calcDeflection(M, EI, delx, solvedInitRotation, initDef)    
        print('Solved initial rotation is {one}'.format(one=solvedInitRotation))
        print('The error in displacement at support is {one}'.format(one=Deflection[supportIndexB])) 
            
        #Define the layout object
        layout = go.Layout(
            title={
                'text': "Deflection",
                'y':0.85,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'},
            titlefont=dict(size=15),
            yaxis = dict(
                title='Deflection'
            ),
            xaxis = dict(
                title='Distance (m)',       
                range=[-1, span+1]
            ),
            showlegend=False,        
        )

        #Define the shear force trace
        line = go.Scatter(
            x = X,
            y = Deflection,
            mode='lines',
            name='Deflection',
            line_color='orange',
            # fill='tonexty',
            fillcolor='rgba(255, 255, 0, 0.1)'
        )

        #Define a horizontal line to represent the structure
        axis = go.Scatter(
            x = [0, span],
            y = [0,0],
            mode='lines',
            line_color='black'
        )
        #Generate and view the figure
        fig = go.Figure(data=[line,axis], layout=layout)
        st.plotly_chart(fig)   

        max_deflection=abs(min(Deflection))
        st.write("Maximum vertical deflection in downward direction:",round(max_deflection,5), " m")

# Function for Overhang Beam Analysis
def overhang_beam():
    st.title("Overhang Beam Analysis Tool")

    # Input parameters
    span = st.number_input("Enter the span of the beam in meters:")
    st.subheader("Ensure units are in 'm' and 'kN'.")
    A = st.number_input("Enter the span of the beam from left side till the 1st support (A) in meters:")
    B = st.number_input("Enter the span of the beam from left side till the 2nd support (B) in meters:")   
    st.subheader("Point Loads")
    pointloads = st.text_input("Format for Point Loads: [[location 1,Magnitude 1], [location 2,Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Point Moments")
    pointmoments = st.text_input("Format for Point Moments: [[location 1,Magnitude 1], [location 2,Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Distributed Loads")
    distributedloads = st.text_input("Format for Distributed Loads: [[start location 1, end location 1,Magnitude 1], [start location 2, end location 2,Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Linear Loads")
    linearloads = st.text_input("Format for Linear Loads: [[start location 1, end location 1, start Magnitude 1, end Magnitude 1], [start location 2, end location 2, start Magnitude 2, end Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Trapezoidal Loads")
    trapezoidalloads = st.text_input("Format for Trapezoidal Loads: [[start location 1, end location 1, start Magnitude 1, end Magnitude 1], [start location 2, end location 2, start Magnitude 2, end Magnitude 2], ...]:", value="[[]]")
    
    st.subheader("Young's Modulus")
    E = st.number_input("Young's Modulus of the material in Pascals (N/m²):", min_value=0, max_value=1000000000, value=200000000)
    
    st.subheader("Moment of Inertia")
    I = st.number_input("Moment of Inertia of the beam in meters (m⁴):", min_value=0.0, max_value=1000000000.0, value=0.0020833, step=0.0000001, format="%.7f")


    # Convert input text to numpy arrays
    def convert_to_list_of_lists(input_str):
        try:
            # Use ast.literal_eval to safely evaluate the string
            result = ast.literal_eval(input_str)
            # Check if the result is a list of lists
            if isinstance(result, list) and all(isinstance(i, list) for i in result):
                return result
            else:
                st.error("Input is not a valid list of lists.")
                return None
        except (ValueError, SyntaxError):
            st.error("Invalid input format.")
            return None
        
    if st.button("Calculate"):
        pointloads = np.array(convert_to_list_of_lists(pointloads))
        pointmoments = np.array(convert_to_list_of_lists(pointmoments))
        distributedloads = np.array(convert_to_list_of_lists(distributedloads))
        linearloads = np.array(convert_to_list_of_lists(linearloads))
        trapezoidalloads = np.array(convert_to_list_of_lists(trapezoidalloads))
        
        delta = 0.00005
        X = np.arange(0, span+delta, delta)
        nPL = len(pointloads[0])
        nPM = len(pointmoments[0])
        nUDL = len(distributedloads[0])
        nLDL = len(linearloads[0])
        nTDL = len(trapezoidalloads[0])
            
        reactions = np.array([0.0, 0.0])
        shearForce = np.empty([0, len(X)])
        bendingMoment = np.empty([0, len(X)])
        
        
        def reactions_PL(n):
            xp = pointloads[n, 0]
            fy = pointloads[n, 1]

            la_p = A-xp
            mp = fy*la_p
            la_vb = B-A
            Vb = mp/la_vb
            Va = -fy-Vb

            return Va, Vb

        def reactions_PM(n):
            xm = pointmoments[n, 0]
            m = pointmoments[n, 1]
            la_vb = B-A
            
            Vb = m/la_vb
            Va = -Vb
            
            return Va, Vb

        def reactions_UDL(n):
            xStart = distributedloads[n,0]
            xEnd = distributedloads[n,1]
            fy = distributedloads[n,2]
            
            fy_Res = fy*(xEnd - xStart)
            x_Res = xStart + 0.5*(xEnd - xStart)
            
            la_p = A-x_Res
            mp = fy_Res*la_p
            la_vb = B-A
            
            Vb = mp/la_vb
            Va = -fy_Res -Vb
            
            return Va, Vb
        
        def reactions_LDL(n):
            xStart = linearloads[n,0]
            xEnd = linearloads[n,1]
            fy_start = linearloads[n,2]
            fy_end = linearloads[n,3]
            
            #Determine location and magnitude of force resultant
            if abs(fy_start)>0:
                fy_Res = 0.5*fy_start*(xEnd-xStart)
                x_Res = xStart + (1/3)*(xEnd-xStart)
            else:
                fy_Res = 0.5*fy_end*(xEnd-xStart)
                x_Res = xStart + (2/3)*(xEnd-xStart)
            
            la_p = A-x_Res
            mp = fy_Res*la_p
            la_vb = B-A
            
            Vb = mp/la_vb
            Va = -fy_Res -Vb
            
            return Va, Vb
        
        def reactions_TDL(n):
            xStart = trapezoidalloads[n,0]
            xEnd = trapezoidalloads[n,1]
            fy_start = trapezoidalloads[n,2]
            fy_end = trapezoidalloads[n,3]
            
            #determine loacation and magnitude of force resultant
            if abs(fy_start) > abs(fy_end):
                fy_Res = 0.5*(xEnd-xStart)*(fy_start+fy_end)
                x_Res = xStart + ((1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end)))
            else:
                fy_Res = 0.5*(xEnd-xStart)*(fy_start+fy_end)
                x_Res = xStart + ((1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end)))
                
            la_p = A-x_Res
            mp = fy_Res*la_p
            la_vb = B-A
            Vb = mp/la_vb
            Va = -fy_Res -Vb
            
            return Va, Vb

        PL_record = np.empty([0,2])
        if(nPL>0):
            for n, p in enumerate(pointloads):
                va, vb =reactions_PL(n)
                PL_record = np.append(PL_record, [np.array([va, vb])], axis =0)

                reactions[0] = reactions[0] + va
                reactions[1] = reactions[1] + vb

        PM_record = np.empty([0,2])
        if(nPM>0):
            for n, p in enumerate(pointmoments):
                va, vb = reactions_PM(n)
                PM_record = np.append(PM_record, [np.array([va, vb])], axis =0)
                
                reactions[0] = reactions[0] + va
                reactions[1] = reactions[1] + vb
                
        UDL_record = np.empty([0,2])
        if(nUDL>0):
            for n, p in enumerate(distributedloads):
                va, vb = reactions_UDL(n)
                UDL_record = np.append(UDL_record, [np.array([va, vb])], axis =0)
                
                reactions[0] = reactions[0] + va
                reactions[1] = reactions[1] + vb
        
        LDL_record = np.empty([0,2])
        if(nLDL>0):
            for n, p in enumerate(linearloads):
                va, vb = reactions_LDL(n)
                LDL_record = np.append(LDL_record, [np.array([va, vb])], axis =0)
                
                reactions[0] = reactions[0] + va
                reactions[1] = reactions[1] + vb
                
        TDL_record = np.empty([0,2])
        if(nTDL>0):
            for n, p in enumerate(trapezoidalloads):
                va, vb = reactions_TDL(n)
                TDL_record = np.append(TDL_record, [np.array([va, vb])], axis =0)
                
                reactions[0] = reactions[0] + va
                reactions[1] = reactions[1] + vb        

        #Plotting and printing

        st.write('The vertical reaction at A is {} kN'.format(round(reactions[0], 2)))
        st.write('The vertical reaction at B is {} kN'.format(round(reactions[1], 2)))
        
                ##          SHEAR AND MOMENT CALCULATIONS
        #Define function to calculate shear forces and bending moments due to point loads

        #7
        def shear_moment_PL(n):
            xp = pointloads[n, 0]
            fy = pointloads[n, 1]
            Va = PL_record[n, 0]
            Vb = PL_record[n, 1]

            #cycle through the structure and calculate the shearforce and bending moment at each point
            Shear = np.zeros(len(X))
            Moment  = np.zeros(len(X))

            for i, x in enumerate(X):
                shear = 0
                moment = 0

                if x > A:
                    #calculate shear and moment due to reaction at A
                    shear = shear + Va
                    moment = moment - Va*(x-A)
                if x > B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                if x > xp:
                    # calculate shear and moment due to reaction at point load
                    shear = shear + fy
                    moment = moment - fy*(x - xp)
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
                #Define function to calculate shear forces and bending moments due to point moments
        def shear_moment_PM(n):
            xm = pointmoments[n, 0]
            m = pointmoments[n, 1]
            Va = PM_record[n, 0]
            Vb = PM_record[n, 1]
            
            #cycle through the structure and calculate the shearforce and bending moment at each point
            Shear = np.zeros(len(X))
            Moment  = np.zeros(len(X))

            for i, x in enumerate(X):
                shear = 0
                moment = 0

                if x > A:
                    #calculate shear and moment due to reaction at A
                    shear = shear + Va
                    moment = moment - Va*(x - A)
                if x > B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                if x > xm:
                    # calculate shear and moment due to reaction at point load
                    moment = moment - m
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
        def shear_moment_UDL(n):
            xStart = distributedloads[n,0]
            xEnd = distributedloads[n,1]
            fy = distributedloads[n,2]
            Va = UDL_record[n,0]
            Vb = UDL_record[n,1]
            
            #cycle through the structure and calculate the shearforce and bending moment at each point
            Shear = np.zeros(len(X))
            Moment  = np.zeros(len(X))

            for i, x in enumerate(X):
                shear = 0
                moment = 0

                if x > A:
                    #calculate shear and moment due to reaction at A
                    shear = shear + Va
                    moment = moment - Va*(x-A)
                if x > B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x - B)
                if x>xStart and x<=xEnd:
                    # calculate shear and moment due to reaction at point load
                    shear = shear + fy*(x-xStart)
                    moment = moment - fy*(x-xStart)*0.5*(x-xStart)
                elif (x>xEnd):
                    shear = shear + fy*(xEnd-xStart)
                    moment = moment - fy*(xEnd-xStart)*(x-xStart-0.5*(xEnd-xStart))
                    
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
                #Define function to calculate shear forces and bending moments due to LDLs
        def shear_moment_LDL(n):
            xStart = linearloads[n,0]
            xEnd = linearloads[n,1]
            fy_start = linearloads[n,2]
            fy_end = linearloads[n,3]
            Va = LDL_record[n,0]
            Vb = LDL_record[n,1]
            
            Shear = np.zeros(len(X))
            Moment = np.zeros(len(X))
            
            for i, x in enumerate(X):
                shear = 0
                moment = 0
                
                if x > A:
                    #calculate shear and moment due to reaction at A
                    shear = shear + Va
                    moment = moment - Va*(x-A)
                if x > B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x-B)
                if x>xStart and x<=xEnd:
                    if abs(fy_start)>0:
                        x_base = x-xStart
                        f_cut = fy_start - x_base*(fy_start/(xEnd-xStart))
                        R1 = 0.5*x_base*(fy_start-f_cut)
                        R2 = x_base*f_cut
                        shear = shear + R1 + R2
                        moment = moment - R1*(2/3)*x_base - R2*(x_base/2)
                        
                    else:
                        x_base  = x-xStart
                        f_cut = fy_end*(x_base/(xEnd-xStart))
                        R = 0.5*x_base*f_cut
                        shear = shear + R
                        moment = moment -R*(x_base/3)
                elif x>xEnd:
                    if abs(fy_start)>0:
                        R = 0.5*fy_start*(xEnd-xStart)
                        xr = xStart + (1/3)*(xEnd-xStart)
                        shear = shear + R
                        moment = moment - R*(x-xr)
                    else:
                        R=0.5*fy_end*(xEnd-xStart)
                        xr = xStart + (2/3)*(xEnd-xStart)
                        shear = shear + R
                        moment = moment - R*(x-xr)
                        
                # store shear and moment for this location
                Shear[i] = shear
                Moment[i] = moment
            return Shear, Moment
        
                #Define function to calculate shear forces and bending moments due to TDLs
        def shear_moment_TDL(n):
            xStart = trapezoidalloads[n,0]
            xEnd = trapezoidalloads[n,1]
            fy_start = trapezoidalloads[n,2]
            fy_end = trapezoidalloads[n,3]
            Va = TDL_record[n,0]
            Vb = TDL_record[n,1]
            
            Shear = np.zeros(len(X))
            Moment = np.zeros(len(X))
            
            for i, x in enumerate(X):
                shear = 0
                moment = 0
                
                if x > A:
                    #calculate shear and moment due to reaction at A
                    shear = shear + Va
                    moment = moment - Va*(x-A)
                if x > B:
                    # calculate shear and moment due to reaction at B
                    shear = shear + Vb
                    moment = moment - Vb*(x-B)
                if x>xStart and x<=xEnd:
                    if abs(fy_start)>abs(fy_end):
                        x_base = x-xStart
                        f_cut=fy_start-x_base*((fy_start-fy_end)/(xEnd-xStart))
                        
                        R1= 0.5*(x_base)*(fy_start-f_cut)
                        la_R1=(2/3)*(x_base)
                        
                        R2= x_base*f_cut
                        la_R2=0.5*(x_base)
                        
                        shear = shear + R1 + R2
                        moment = moment - R1*la_R1 - R2*la_R2
                        
                    else:
                        x_base = x-xStart
                        f_cut=fy_start+x_base*((fy_end-fy_start)/(xEnd-xStart))
                        
                        R1= 0.5*(x_base)*(f_cut-fy_start)
                        la_R1=(1/3)*(x_base)
                        
                        R2= (x_base)*(fy_start)
                        la_R2=0.5*(x_base)
                        
                        shear = shear + R1 + R2
                        moment = moment - R1*la_R1 - R2*la_R2
                    
                elif x>xEnd:
                    if abs(fy_start)>abs(fy_end):
                        R = 0.5*(xEnd-xStart)*(fy_start+fy_end)
                        xr = xStart + (1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end))
                        
                        shear = shear + R
                        moment = moment - R*(x-xr)
                    
                    else:
                        R=0.5*(xEnd-xStart)*(fy_end+fy_start)
                        xr=xStart + (1/3)*(((xEnd-xStart)*(fy_start+2*fy_end))/(fy_start+fy_end))
                        
                        shear = shear + R
                        moment = moment - R*(x-xr)
                
                Shear[i] = shear
                Moment[i] = moment

            return Shear, Moment
        
        #Cycle through all the point loads and determine shear and moment

        #8
        if(nPL>0):
            for n, p in enumerate(pointloads):
                Shear, Moment = shear_moment_PL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each point load
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each point load
        
        if(nPM>0):
            for n, p in enumerate(pointmoments):
                Shear, Moment = shear_moment_PM(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each point moments
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each point moments        

        if(nUDL>0):
            for n, p in enumerate(distributedloads):
                Shear, Moment = shear_moment_UDL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each UDL
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each UDL        

        if(nLDL>0):
            for n, p in enumerate(linearloads):
                Shear, Moment = shear_moment_LDL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each LDL
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each LDL        

        if(nTDL>0):
            for n, p in enumerate(trapezoidalloads):
                Shear, Moment = shear_moment_TDL(n)
                shearForce = np.append(shearForce, [Shear], axis=0)        #Store shear force record for each TDL
                bendingMoment = np.append(bendingMoment, [Moment], axis=0) #Store bending moment record for each TDL

        #Define the layout object
        layout = go.Layout(
            title={
                'text': "Shear Force Diagram",
                'y':0.85,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'},
            titlefont=dict(size=15),
            yaxis = dict(
                title='Shear Force (kN)'
            ),
            xaxis = dict(
                title='Distance (m)',       
                range=[-1, span+1]
            ),
            showlegend=False,        
        )

        #Define the shear force trace
        line = go.Scatter(
            x = X,
            y = sum(shearForce),
            mode='lines',
            name='Shear Force',
            fill='tonexty',
            line_color='green',
            fillcolor='rgba(0, 255, 0, 0.1)'
        )

        #Define a horizontal line to represent the structure
        axis = go.Scatter(
            x = [0, span],
            y = [0,0],
            mode='lines',
            line_color='black'
        )

        #Generate and view the figure
        fig = go.Figure(data=[line, axis], layout=layout)
        st.plotly_chart(fig)
      
      
        
        layout = go.Layout(
        title={
            'text': "Bending Moment Diagram",
            'y':0.85,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        titlefont=dict(size=15),
        yaxis = dict(
            title='Bending Moment (kNm)',
            autorange="reversed",   
        ),
        xaxis = dict(
            title='Distance (m)',       
            range=[-1, span+1]
        ),
        showlegend=False 
    )

        #Define the Bending Moment trace
        line = go.Scatter(
            x = X,
            y = -sum(bendingMoment),
            mode='lines',
            name='Bending Moment',
            fill='tonexty',
            line_color='red',
            fillcolor='rgba(255, 0, 0, 0.1)'
        )

        #Define a horizontal line to represent the structure
        axis = go.Scatter(
            x = [0, span],
            y = [0,0],
            mode='lines',
            line_color='black'
        )

        #Generate and view the figure
        fig = go.Figure(data=[line, axis], layout=layout)
        st.plotly_chart(fig)
        
        deltaRot = 0.000005
        initRot = -0.00021       

        M = -sum(bendingMoment)
        delx = X[1]-X[0]
        EI = E*I
        initDef = 0
        supportIndexA = np.where(X==A)[0].item()
        supportIndexB = np.where(X==B)[0].item()
            
        # Define a function to calculate deflection by forward integrating (from left support) using the Trapezoidal Rule

        def calcDeflection(M, EI, delx, theta_0, v_0):
                
            theta_im1 = theta_0
            v_im1 = v_0
                
            Rotation = np.zeros(len(X))
            Deflection = np.zeros(len(X))
            Rotation[supportIndexA] = theta_im1
            Deflection[supportIndexA] = v_im1
                
            #loop through the data and integrate (Trapezoidal Rule)
            for i, m in enumerate(M[supportIndexA::]):
                ind = i+supportIndexA
                if i>0:
                    M_im1 = M[ind-1]
                    M_i = M[ind]
                    M_avg = 0.5*(M_i+M_im1)
                        
                    theta_i = theta_im1 + (M_avg/EI)*delx
                    v_i = v_im1 + 0.5*(theta_i+theta_im1)*delx
                        
                    #store data
                    Rotation[ind] = theta_i
                    Deflection[ind] = v_i
                        
                    #update values for next loop iteration
                    theta_im1 = theta_i
                    v_im1 = v_i
                        
            return Rotation, Deflection
            
        def zeroCrossing(Deflection, guessStep, initRot, initDef):
            """
            Find the value of initial rotation that minimised deflection at right side 
            support by identifying where error crosses zero.
            """
            
            #If the deflection error is positive
            if Deflection[supportIndexB]>0:
                errorIsPositive = True #Set flag for error sign
                
                #Keep testing lower initial rotation values until error turns NEGATIVE
                while errorIsPositive:
                    initRot = initRot + guessStep
                    Rotation, Deflection = calcDeflection(M, EI, delx, initRot, initDef)
                    
                    #If error has turned NEGATIVE, switch the flag to allow loop to stop
                    if Deflection[supportIndexB]<=0:
                        errorIsPositive = False
                        solvedInitRotation = initRot #Save the 'solved' value that minimised the error
            
            #Else if deflection error is negative
            elif (Deflection[supportIndexB]<0):
                errorIsPositive = False #Set flag for error sign
                
                #Keep testing lower initial rotation values until error turns POSITIVE
                while not errorIsPositive:
                    initRot = initRot + guessStep
                    Rotation, Deflection = calcDeflection(M, EI, delx, initRot, initDef)
                    
                    #If error has turned POSITIVE, switch the flag to allow loop to stop
                    if Deflection[supportIndexB]>=0:
                        errorIsPositive = True
                        solvedInitRotation = initRot #Save the 'solved' value that minimised the error
            
            return solvedInitRotation
            
        #Test whether reducing or increasing initial rotation leads to reduction in disp error at other support
        testDef = np.zeros(3)
        for i, r in enumerate([initRot-deltaRot, initRot, initRot+deltaRot]):
            Rotation, Deflection = calcDeflection(M, EI, delx, r, initDef)
            testDef[i] = Deflection[supportIndexB]
                
        if(abs(testDef[0])<abs(testDef[1])):
            #Need to test in the negative rotation direction by reducing the initial rotation guess
            print('Need to test in the negative direction')    
            solvedInitRotation = zeroCrossing(Deflection, -deltaRot, initRot, initDef)            

        elif(abs(testDef[2])<abs(testDef[1])):
            #Need to test in the positive rotation direction by incresing the initial rotation guess 
            print('Need to test in the positive direction')    
            solvedInitRotation = zeroCrossing(Deflection, +deltaRot, initRot, initDef)      

        #Run the deflection calculation with the solved value of initial rotation
        Rotation, Deflection = calcDeflection(M, EI, delx, solvedInitRotation, initDef)    
        print('Solved initial rotation is {one}'.format(one=solvedInitRotation))
        print('The error in displacement at support is {one}'.format(one=Deflection[supportIndexB]))
            
        if A!=0:
            print("There is an overhand on the left side - solve for deflection by integrating in reverse direction")
            
            theta_im1 = -solvedInitRotation #Rotation on other side of support A
            v_im1 = 0 #Vertical deflection at support A

            #Generate a range of indices in reverse direction from support A to left endge of beam
            reverseRange = np.arange(supportIndexA-1,-1,-1) 

            #Loop through data and integrate (Trapezoidal rule) - REVERSE DIRECTION
            for i in reverseRange:                        
                M_im1 = M[i+1] #(300) - Assign previous value of M (reverse direction)
                M_i = M[i] #(299) - Assign current value of M (reverse direction)
                M_avg = 0.5*(M_i + M_im1)   
                
                theta_i = theta_im1 + (M_avg/EI)*delx #Integrate moment values to get rotations                     
                v_i = v_im1 + 0.5*(theta_i+theta_im1)*delx #Integrate rotation values to get displacements
                
                #Store data
                Rotation[i] = theta_i        
                Deflection[i] = v_i

                #Update values for next loop iteration
                theta_im1 = theta_i     
                v_im1 = v_i   
            
        #Define the layout object
        layout = go.Layout(
            title={
                'text': "Deflection",
                'y':0.85,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'},
            titlefont=dict(size=15),
            yaxis = dict(
                title='Deflection'
            ),
            xaxis = dict(
                title='Distance (m)',       
                range=[-1, span+1]
            ),
            showlegend=False,        
        )

        #Define the shear force trace
        line = go.Scatter(
            x = X,
            y = Deflection,
            mode='lines',
            name='Deflection',
            line_color='orange',
          #  fill='tonexty',
            fillcolor='rgba(255, 255, 0, 0.1)'
        )

        #Define a horizontal line to represent the structure
        axis = go.Scatter(
            x = [0, span],
            y = [0,0],
            mode='lines',
            line_color='black'
        )
        #Generate and view the figure
        fig = go.Figure(data=[line,axis], layout=layout)
        st.plotly_chart(fig)   

        max_deflection=abs(min(Deflection))
        st.write("Maximum vertical deflection in downward direction:",round(max_deflection,5), " m")

# Main function to handle page navigation
def main():
    st.title("Statically Determinate Beam Analyzer")
    st.write("This tool makes interactable plots of Shear Force Diagram, Bending Moment Diagram, and Deflection for Cantilever Beam, Simply Supported Beam, and Overhang Beam")
    st.write(' ')

    # Initialize session state if it doesn't exist
    if "page" not in st.session_state:
        st.session_state["page"] = "home"

    # Display beam selection buttons or analysis based on the current page
    if st.session_state["page"] == "home":
        display_home()
    else:
        display_analysis()

def display_home():
    col1,_, col2,_, col3 = st.columns(5)
    cantiliver = Image.open("images/cantiliver.png")
    simply_supported = Image.open("images/simply_supported.png")
    overhang = Image.open("images/overhang.png")
    with col1:
        st.image(cantiliver, use_column_width=True, channels="RGB")
        st.markdown("<div style='text-align: center; margin-top: 10px;'>", unsafe_allow_html=True)
        # st.markdown("<div style='height: 32px; padding: 10rem 20 rem;'></div>", unsafe_allow_html=True)
        st.button("Cantilever Beam", on_click=set_page, args=("cantilever",))
    with col2:
        st.image(simply_supported, use_column_width=True, channels="RGB")
        st.markdown("<div style='text-align: center; margin-top: 10px;'>", unsafe_allow_html=True)
        # st.markdown("<div style='height: 35px;padding-right: 200rem'></div>", unsafe_allow_html=True)
        st.button("Simply Supported Beam", on_click=set_page, args=("simply_supported",))
    with col3:
        st.image(overhang, use_column_width=True, channels="RGB")
        st.markdown("<div style='text-align: center; margin-top: 10px;'>", unsafe_allow_html=True)
        # st.markdown("<div style='height: 30px;'></div>", unsafe_allow_html=True)
        st.button("Overhang Beam", on_click=set_page, args=("overhang",))
        
    
def set_page(page):
    st.session_state["page"] = page

def display_analysis():
    selected_page = st.session_state["page"]
    if selected_page == "cantilever":
        cantilever_beam()
    elif selected_page == "simply_supported":
        simply_supported_beam()
    elif selected_page == "overhang":
        overhang_beam()

    if st.button("Return to Beam Selection", on_click=set_page, args=("home",)):
        st.session_state["page"] = "home"


if __name__ == '__main__':
    main()
