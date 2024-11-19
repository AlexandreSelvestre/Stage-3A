# print(1)
# import numpy as np
# import datetime 
# print(2, datetime.datetime.now())
# array = np.array([1, 2, 3, 4, 5])
# other = np.array([0,2,3,5,6])
# dic = {0: array}
# array = other
# dic[1] = array
# print(dic, datetime.datetime.now())


#from __future__ import print_function
print("1")
from datetime import datetime
from scipy.interpolate import PchipInterpolator
from scipy.integrate import quad, simpson
import os
from radiomics import featureextractor, getFeatureClasses
print("2", datetime.now())
import radiomics
import SimpleITK as sitk
import glob
import pandas as pd
import matplotlib.pyplot as plt
print("3", datetime.now())
import tqdm #dejafe
import numpy as np
import itertools
print("4", datetime.now())
from utils import mask_superpose_simple, resample_image_to_reference
from utils_advanced_area import *
if __name__ == '__main__':
    current_dir = os.getcwd()
    import matplotlib
    if current_dir != "/gpfs/users/selvestra/python_codes":
        matplotlib.use('TKAgg')
from utils_recaler import eliminate_temps
from natsort import natsorted
print(5, datetime.now())