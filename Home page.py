import streamlit as st
import numpy as np
import altair as alt
import pandas as pd

import vegas
import random
import time
import math
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from statistics import mean
from numpy import arange
from pandas import read_csv
from matplotlib import pyplot
from PIL import Image


from plotly.subplots import make_subplots
import plotly.graph_objects as go


st.title('Welcome to Virialis!')


st.write("Below, you can see a diagram on how to use Virialis")
image = Image.open('Fitting (2).png')

st.image(image)