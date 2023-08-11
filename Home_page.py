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
from PIL import Image


from plotly.subplots import make_subplots
import plotly.graph_objects as go


st.title('Welcome to Virialis!')


st.write("Below, you can see a diagram on how to use Virialis")
image = Image.open('Fitting (2).png')

st.image(image)
