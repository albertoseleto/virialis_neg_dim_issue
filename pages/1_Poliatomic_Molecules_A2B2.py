# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:00:27 2023

@author: Alberto
"""
#########################
#LIBRARY IMPORTS
import streamlit as st
import numpy as np
import pandas as pd
import vegas
import time
from scipy.optimize import curve_fit
from statistics import mean
from numpy import arange
from pandas import read_csv
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from constant import constants

###############################

# TIME MEASURE

# DATA ENTRY
st.title('Second Virial Coefficient Calculator for A2B2 Molecule')

uploaded_file = st.file_uploader("upload a file")


potential = st.selectbox(
    'What potential energy do you want to use?',
    ('Rydberg Potential', 'Improved Leonard-Jonnes Potential'))


step = st.selectbox(
    'What is the Temperature step that you want to use?',
    (100, 50, 25, 200, 300))


st.write('The calculus will be made using ', potential,
          'using a temperature step of', step)


def ILJ():
    st.write('LEONARD-JONNES IS IN CONSTRUCTION...')

    data = pd.read_csv(uploaded_file, sep="\s+", header=None)
    data.columns = ["alpha", "beta", "mp", "De", "Req" ]
    st.subheader('DataFrame')
    st.write(data)


    #dataframes de input
    h_alpha = data.loc[0, 'alpha']
    h_beta = data.loc[0, 'beta']
    h_mp = data.loc[0, 'mp']
    h_De = data.loc[0, 'De'] 
    h_Req = data.loc[0, 'Req']
        

    x_alpha = data.loc[1, 'alpha']
    x_beta = data.loc[1, 'beta']
    x_mp = data.loc[1, 'mp']
    x_De = data.loc[1, 'De'] 
    x_Req = data.loc[1, 'Req']
        

    z_alpha = data.loc[2, 'alpha']
    z_beta = data.loc[2, 'beta']
    z_mp = data.loc[2, 'mp']
    z_De = data.loc[2, 'De'] 
    z_Req = data.loc[2, 'Req']
        

    ta_alpha = data.loc[3, 'alpha']
    ta_beta = data.loc[3, 'beta']
    ta_mp = data.loc[3, 'mp']
    ta_De = data.loc[3, 'De'] 
    ta_Req = data.loc[3, 'Req']
    

    tb_alpha = data.loc[4, 'alpha']
    tb_beta = data.loc[4, 'beta']
    tb_mp = data.loc[4, 'mp']
    tb_De = data.loc[4, 'De'] 
    tb_Req = data.loc[4, 'Req']
        

    l_alpha = data.loc[5, 'alpha']
    l_beta = data.loc[5, 'beta']
    l_mp = data.loc[5, 'mp']
    l_De = data.loc[5, 'De'] 
    l_Req = data.loc[5, 'Req']
        


    #def de funcoes
    class LC1: 
        def UH(r):
            n = h_beta + h_alpha * (r / h_Req) ** 2
            return  h_De * ((h_mp/(n - h_mp) * (h_Req/r) ** n) - (n/(n - h_mp) * (h_Req/r) ** h_mp))

        def UX(r):

            n = x_beta + x_alpha * (r / x_Req) ** 2
            return x_De * ((x_mp/(n - x_mp) * (x_Req/r) ** n) - (n/(n - x_mp) * (x_Req/r) ** x_mp))

        def UZ(r):

            n = z_beta + z_alpha * (r / z_Req) ** 2
            return z_De * ((z_mp/(n - z_mp) * (z_Req/r) ** n) - (n/(n - z_mp) * (z_Req/r) ** z_mp))

        def UTa(r):
            n = ta_beta + ta_alpha * (r / ta_Req) ** 2
            return ta_De * ((ta_mp/(n - ta_mp) * (ta_Req/r) ** n) - (n/(n - ta_mp) * (ta_Req/r) ** ta_mp))

        def UTb(r):
            n = tb_beta + tb_alpha * (r / tb_Req) ** 2
            return tb_De * ((tb_mp/(n - tb_mp) * (tb_Req/r) ** n) - (n/(n - tb_mp) * (tb_Req/r) ** tb_mp))

        def UL(r):
            n = l_beta + l_alpha * (r / l_Req) ** 2
            return l_De * ((l_mp/(n - l_mp) * (l_Req/r) ** n) - (n/(n - l_mp) * (l_Req/r) ** l_mp))
        '''
        def USa(r):
            n = Sa.beta + Sa.alpha * (r / Sa.Req) ** 2
            return Sa.De * ((Sa.mp/(n - Sa.mp) * (Sa.Req/r) ** n) - (n/(n - Sa.mp) * (Sa.Req/r) ** Sa.mp))

        def USb(r):
            n = Sb.beta + Sb.alpha * (r / Sb.Req) ** 2
            return Sb.De * ((Sb.mp/(n - Sb.mp) * (Sb.Req/r) ** n) - (n/(n - Sb.mp) * (Sb.Req/r) ** Sb.mp))
        '''

    class complexa1:
        def UM_000(r):
            UM_000 = (2*LC1.UH(r)+LC1.UL(r) + 2 *
                        (LC1.UTa(r)+LC1.UTb(r)+LC1.UX(r)))/9
            return UM_000

        def UM_202(r):
            UM_202 = 2*(LC1.UH(r)-LC1.UL(r)+LC1.UTa(r) -
                        2*LC1.UTb(r)+LC1.UX(r))/(9*(5**(1/2)))
            return UM_202

        def UM_022(r):
            UM_022 = 2*(LC1.UH(r)-LC1.UL(r)-2*LC1.UTa(r) +
                        LC1.UTb(r)+LC1.UX(r))/(9*(5**(1/2)))
            return UM_022

        def UM_220(r):
            UM_220 = 2*(4*LC1.UH(r)-LC1.UL(r)-5*(LC1.UTa(r) + 
                        LC1.UTb(r) + LC1.UX(r))+12*LC1.UZ(r))/(45*(5**(1/2)))
            return UM_220

        def UM_222(r):
            UM_222 = ((2/7)**(1/2))*(13*LC1.UH(r)-LC1.UL(r)+7 *
                                        (LC1.UTa(r)+LC1.UTb(r)-2*LC1.UX(r))-12*LC1.UZ(r))/45
            return UM_222

        def UM_224(r):
            UM_224 = ((2/35)**(1/2)*8 * LC1.UH(r)+LC1.UL(r)+2*LC1.UZ(r))/15
            
            return UM_224

    def UM_FINAL1(r, th_a, th_b, phi):  # potencial que vai ser usado na equação 3
        UM_FINAL = (complexa1.UM_000(r)+(5**(1/2))/4 * (3*(np.cos(2*th_a))+1) * complexa1.UM_202(r) 
            + (5**(1/2))/4 * (3*(np.cos(2*th_b))+1) * complexa1.UM_022(r) 
            + (5**(1/2))/16 * complexa1.UM_220(r) * ((3*(np.cos(2*th_a))+1) * (3*(np.cos(2*th_b))+1) 
            + 12*np.sin(2*th_a)* np.sin(2*th_b)*np.cos(phi)
            + 3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)) 
            - 14**(1/2)*5/112 * complexa1.UM_222(r) * ((3*np.cos(2*th_a))+1) * (3*(np.cos(2*th_b))+1) 
            + 6*np.sin(2*th_a)*np.sin(2*th_b) * np.cos(phi) 
            - 3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b))*np.cos(2*phi)
            + (3*(70)**(1/2))/112 * complexa1.UM_224(r) * ((3*(np.cos(2*th_a))+1)*(3*(np.cos(2*th_b))+1) 
            - 8* np.sin(2*th_a) * np.sin(2*th_b)*np.cos(phi)+((1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi))/2))
        #print('r:',r ,'000:',UM_000(r),'202:',UM_202(r),'022:',UM_022(r),'220:',UM_220(r),'222:',UM_222(r),'224:',UM_224(r))
        #st.write('um',UM_FINAL)
        return UM_FINAL
    

    #funcao final a ser integrada
    def integrand_vegas(x):
        r = x[0]
        th_a = x[1]
        th_b = x[2]
        phi = x[3]
        F = UM_FINAL1(r, th_a, th_b, phi)*constants.rk/T

        return constants.N_A/4 * np.sin(th_a)  * np.sin(th_b)*(r**2)*(1-np.exp(-F))

    
    #integrador e def de limites de integracao
    integ = vegas.Integrator(
        [[0, 10], [0, np.pi], [0, np.pi], [0, 2*np.pi]])

    B_clas = []
    Tstr = []

    B_main = []
    
    for T in range(50, 1000, step):

        integ(integrand_vegas, nitn=10, neval=10000)

        

        result = integ(integrand_vegas, nitn=10, neval=10000)

        # st.write(result.summary())

        st.write('result of classic virial = ', result, 'for T = ', T)
        
        B_clas.append(result.mean)
        Tstr.append(T)



    st.write('result of final virial =', result.mean, 'for T = ', T)

    B_main.append(result.mean )



    r = np.linspace(0, 10, 100)
    th_a = np.linspace(0, np.pi, 100)
    th_b = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2*np.pi, 100)

    fig1 = go.Figure()

    fig1.add_trace(
        go.Scatter(x=r, y=LC1.UH(r), name='UH')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC1.UX(r), name='UX')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC1.UZ(r), name='UZ')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC1.UTa(r), name='UTa')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC1.UTb(r), name='UTb')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC1.UL(r), name='UL')
    )


    fig2 = go.Figure()
    fig2.add_trace(
        go.Scatter(x=r, y=complexa1.UM_000(r), name='UM_000')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=complexa1.UM_202(r), name='UM_202')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=complexa1.UM_220(r), name='UM_220')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=complexa1.UM_022(r), name='UM_022')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=complexa1.UM_222(r), name='UM_222')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=complexa1.UM_224(r), name='UM_224')
    )


    fig3 = go.Figure()



    fig3.add_trace(
        go.Scatter(x=Tstr, y=B_clas, name='B classic')
    )

    fig1.update_xaxes(range=[0, 10])
    fig1.update_yaxes(range=[-3*h_De, 3*h_De])

    fig2.update_xaxes(range=[0, 10])
    fig2.update_yaxes(range=[-3*h_De, 5*h_De])



    fig1.update_layout(height=600, width=800,
                    title_text="Graph of Leading Configurations per Distance", xaxis_title=r"$r[A]$",
                    yaxis_title='U[kcal/mol]'
                        )
    st.plotly_chart(fig1, use_container_width=True)

    fig2.update_layout(height=600, width=800,
                    title_text="Graph of Complex Energies per Distance", xaxis_title=r"$r[A]$",
                    yaxis_title='U[kcal/mol]')
    st.plotly_chart(fig2, use_container_width=True)

    fig3.update_layout(height=600, width=800,
                    title_text="Graph of the Second Virial Coefficient B per Temperature T", yaxis_title='B(T)[cm^3/mol]',
                        xaxis_title="T[K]")
    st.plotly_chart(fig3, use_container_width=True)


    st.write('evaluation time {}'.format(time.strftime(
        "%H:%M:%S", time.gmtime(time.time()-start))))
    st.write('number of temperatures:', len(Tstr))

    st.success('Calculations finished!', icon="✅")


def ryd():
    data = pd.read_csv(uploaded_file, sep="\s+", header=None)

    data.columns = ["a1", "a2", "a3", "a4", 'a5',
                    'De', 'Req', 'Eref']
    st.subheader('DataFrame')
    st.write(data)
    #st.subheader('Descriptive Statistics')
    # st.write(data.describe())

    class H:

        a1 = data.loc[0, 'a1']
        a2 = data.loc[0, 'a2']
        a3 = data.loc[0, 'a3']
        a4 = data.loc[0, 'a4']
        a5 = data.loc[0, 'a5'] 
        De = data.loc[0, 'De'] 
        Req = data.loc[0, 'Req']
        Eref = data.loc[0, 'Eref']

    class X:

        a1 = data.loc[1, 'a1']
        a2 = data.loc[1, 'a2']
        a3 = data.loc[1, 'a3']
        a4 = data.loc[1, 'a4']
        a5 = data.loc[1, 'a5']
        De = data.loc[1, 'De']  # eV
        Req = data.loc[1, 'Req']
        Eref = data.loc[1, 'Eref']

    class Z:

        a1 = data.loc[2, 'a1']
        a2 = data.loc[2, 'a2']
        a3 = data.loc[2, 'a3']
        a4 = data.loc[2, 'a4']
        a5 = data.loc[2, 'a5']
        De = data.loc[2, 'De']
        Req = data.loc[2, 'Req']
        Eref = data.loc[2, 'Eref']

    class Ta:

        a1 = data.loc[3, 'a1']
        a2 = data.loc[3, 'a2']
        a3 = data.loc[3, 'a3']
        a4 = data.loc[3, 'a4']
        a5 = data.loc[3, 'a5']
        De = data.loc[3, 'De']
        Req = data.loc[3, 'Req']
        Eref = data.loc[3, 'Eref']

    class Tb:

        a1 = data.loc[4, 'a1']
        a2 = data.loc[4, 'a2']
        a3 = data.loc[4, 'a3']
        a4 = data.loc[4, 'a4']
        a5 = data.loc[4, 'a5']
        De = data.loc[4, 'De']
        Req = data.loc[4, 'Req']
        Eref = data.loc[4, 'Eref']

    class L:

        a1 = data.loc[5, 'a1']
        a2 = data.loc[5, 'a2']
        a3 = data.loc[5, 'a3']
        a4 = data.loc[5, 'a4']
        a5 = data.loc[5, 'a5']
        De = data.loc[5, 'De']
        Req = data.loc[5, 'Req']
        Eref = data.loc[5, 'Eref']

    class LC:
        def UH(r):
            y = (r-H.Req)/H.Req
            UH = -H.De*(1 + H.a1*y + H.a2*y**2 + H.a3*y**3 +
                        H.a4*y**4 + H.a5*y**5) * np.exp(-H.a1*y) + H.Eref
            return UH

        def UX(r):

            y = (r-X.Req)/X.Req
            UX = -X.De*(1 + X.a1*y + X.a2*y**2 + X.a3*y**3 +
                        X.a4*y**4 + X.a5*y**5) * np.exp(-X.a1*y) + X.Eref
            return UX

        def UZ(r):

            y = (r-Z.Req)/Z.Req
            UZ = -Z.De*(1 + Z.a1*y + Z.a2*y**2 + Z.a3*y**3 +
                        Z.a4*y**4 + Z.a5*y**5) * np.exp(-Z.a1*y) + Z.Eref
            return UZ

        def UTa(r):
            y = (r-Ta.Req)/Ta.Req
            UTa = -Ta.De*(1 + Ta.a1*y + Ta.a2*y**2 + Ta.a3*y**3 +
                            Ta.a4*y**4 + Ta.a5*y**5) * np.exp(-Ta.a1*y) + Ta.Eref
            return UTa

        def UTb(r):
            y = (r-Tb.Req)/Tb.Req
            UTb = -Tb.De*(1 + Tb.a1*y + Tb.a2*y**2 + Tb.a3*y**3 +
                            Tb.a4*y**4 + Tb.a5*y**5) * np.exp(-Tb.a1*y) + Tb.Eref
            return UTb

        def UL(r):
            y = (r-L.Req)/L.Req
            UL = -L.De*(1 + L.a1*y + L.a2*y**2 + L.a3*y**3 +
                        L.a4*y**4 + L.a5*y**5) * np.exp(-L.a1*y) + L.Eref
            return UL

    class complexa:
        def UM_000(r):
            UM_000 = (2*LC.UH(r)+LC.UL(r) + 2 *
                        (LC.UTa(r)+LC.UTb(r)+LC.UX(r)))/9
            return UM_000

        def UM_202(r):
            UM_202 = 2*(LC.UH(r)-LC.UL(r)+LC.UTa(r) -
                        2*LC.UTb(r)+LC.UX(r))/(9*(5**(1/2)))
            return UM_202

        def UM_022(r):
            UM_022 = 2*(LC.UH(r)-LC.UL(r)-2*LC.UTa(r) +
                        LC.UTb(r)+LC.UX(r))/(9*(5**(1/2)))
            return UM_022

        def UM_220(r):
            UM_220 = 2*(4*LC.UH(r)-LC.UL(r)-5*(LC.UTa(r) + 
                        LC.UTb(r) + LC.UX(r))+12*LC.UZ(r))/(45*(5**(1/2)))
            return UM_220

        def UM_222(r):
            UM_222 = ((2/7)**(1/2))*(13*LC.UH(r)-LC.UL(r)+7 *
                                        (LC.UTa(r)+LC.UTb(r)-2*LC.UX(r))-12*LC.UZ(r))/45
            return UM_222

        def UM_224(r):
            UM_224 = ((2/35)**(1/2)*8 * LC.UH(r)+LC.UL(r)+2*LC.UZ(r))/15
            return UM_224

    def UM_FINAL(r, th_a, th_b, phi):  # potencial que vai ser usado na equação 3
        UM_FINAL = (complexa.UM_000(r)+(5**(1/2))/4 * (3*(np.cos(2*th_a))+1) * complexa.UM_202(r) 
            + (5**(1/2))/4 * (3*(np.cos(2*th_b))+1) * complexa.UM_022(r) 
            + (5**(1/2))/16 * complexa.UM_220(r) * ((3*(np.cos(2*th_a))+1) * (3*(np.cos(2*th_b))+1) 
            + 12*np.sin(2*th_a)* np.sin(2*th_b)*np.cos(phi)
            + 3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)) 
            - 14**(1/2)*5/112 * complexa.UM_222(r) * ((3*np.cos(2*th_a))+1) * (3*(np.cos(2*th_b))+1) 
            + 6*np.sin(2*th_a)*np.sin(2*th_b) * np.cos(phi) 
            - 3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b))*np.cos(2*phi)
            + (3*(70)**(1/2))/112 * complexa.UM_224(r) * ((3*(np.cos(2*th_a))+1)*(3*(np.cos(2*th_b))+1) 
            - 8* np.sin(2*th_a) * np.sin(2*th_b)*np.cos(phi)+((1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi))/2))
        #print('r:',r ,'000:',UM_000(r),'202:',UM_202(r),'022:',UM_022(r),'220:',UM_220(r),'222:',UM_222(r),'224:',UM_224(r))
        # print('um',UM_FINAL)
        return UM_FINAL
    
    #funcao final a ser integrada
    def integrand_vegas(x):
        r = x[0]
        th_a = x[1]
        th_b = x[2]
        phi = x[3]
        F = UM_FINAL(r, th_a, th_b, phi)*constants.rk/T

        return constants.N_A/4 * np.sin(th_a)  * np.sin(th_b)*(r**2)*(1-np.exp(-F))

    
    #integrador e def de limites de integracao
    integ = vegas.Integrator(
        [[0, 10], [0, np.pi], [0, np.pi], [0, 2*np.pi]])

    B_clas = []
    Tstr = []

    B_main = []
    

    for T in range(50, 1000, step):

        integ(integrand_vegas, nitn=10, neval=10000)

        

        result = integ(integrand_vegas, nitn=10, neval=10000)

        # st.write(result.summary())

        st.write('result of classic virial = ', result, 'for T = ', T)
        
        B_clas.append(result.mean)
        Tstr.append(T)



    st.write('result of final virial =', result.mean, 'for T = ', T)

    B_main.append(result.mean )



    r = np.linspace(0, 10, 100)
    th_a = np.linspace(0, np.pi, 100)
    th_b = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2*np.pi, 100)

    fig1 = go.Figure()

    fig1.add_trace(
        go.Scatter(x=r, y=LC.UH(r), name='UH')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC.UX(r), name='UX')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC.UZ(r), name='UZ')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC.UTa(r), name='UTa')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC.UTb(r), name='UTb')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC.UL(r), name='UL')
    )


    fig2 = go.Figure()
    fig2.add_trace(
        go.Scatter(x=r, y=complexa.UM_000(r), name='UM_000')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=complexa.UM_202(r), name='UM_202')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=complexa.UM_220(r), name='UM_220')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=complexa.UM_022(r), name='UM_022')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=complexa.UM_222(r), name='UM_222')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=complexa.UM_224(r), name='UM_224')
    )


    fig3 = go.Figure()



    fig3.add_trace(
        go.Scatter(x=Tstr, y=B_clas, name='B classic')
    )

    fig1.update_xaxes(range=[0, 10])
    fig1.update_yaxes(range=[-3*H.De, 3*H.De])

    fig2.update_xaxes(range=[0, 10])
    fig2.update_yaxes(range=[-3*H.De, 5*H.De])



    fig1.update_layout(height=600, width=800,
                    title_text="Graph of Leading Configurations per Distance", xaxis_title=r"$r[A]$",
                    yaxis_title='U[kcal/mol]'
                        )
    st.plotly_chart(fig1, use_container_width=True)

    fig2.update_layout(height=600, width=800,
                    title_text="Graph of Complex Energies per Distance", xaxis_title=r"$r[A]$",
                    yaxis_title='U[kcal/mol]')
    st.plotly_chart(fig2, use_container_width=True)

    fig3.update_layout(height=600, width=800,
                    title_text="Graph of the Second Virial Coefficient B per Temperature T", yaxis_title='B(T)[cm^3/mol]',
                        xaxis_title="T[K]")
    st.plotly_chart(fig3, use_container_width=True)


    st.write('evaluation time {}'.format(time.strftime(
        "%H:%M:%S", time.gmtime(time.time()-start))))
    st.write('number of temperatures:', len(Tstr))

    st.success('Calculations finished!', icon="✅")

if st.button('Calculate'):
    start = time.time()
    if uploaded_file is not None:
        

        if potential == 'Improved Leonard-Jonnes Potential':
            ILJ()
        elif potential == 'Rydberg Potential':
            ryd()
            
       
    else:
        st.info('☝️ Upload a .dat file')
