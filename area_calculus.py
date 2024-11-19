
from __future__ import print_function
import sys
import os
import logging
import six
from radiomics import featureextractor, getFeatureClasses
import radiomics
import SimpleITK as sitk
import nibabel as nib
import glob
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
import numpy as np
import logging
from utils import mask_superpose_simple, resample_image_to_reference
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('TKAgg')
    
def get_profile(image, mask, spacing):
    #donne la fonction de répatition de la surface en fonction de la profondeur: dictionnaire avec la profondeur en clé et la surface en valeur
    settings = {}
    settings['binWidth'] = 25
    settings['resampledPixelSpacing'] = [2,2]
    settings['interpolator'] = 'sitkBSpline'
    settings['verbose'] = True
    settings['force2D'] = True
    settings['force2Ddimension'] = 2
    extractor_area = featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor_area.disableAllFeatures()
    extractor_area.enableFeaturesByName(**{"shape2D": ["MeshSurface"]})
    new_image, resampler = resample_image_to_reference(image, image, spacing, force_size = None)
    new_mask = resampler.Execute(mask)
    dic_areas = {}
    current_sum = 0
    for i in range(new_image.GetSize()[2]):
        depth = new_image.TransformContinuousIndexToPhysicalPoint([0,0,i])[2]
        try:
            feature = extractor_area.execute(new_image[:,:,i], new_mask[:,:,i])
            area = feature['original_shape2D_MeshSurface']
            current_sum = area*spacing[2] + current_sum
            print(area*spacing[2])
        except Exception as e:
            #print("erreur", e)
            pass
        dic_areas[depth] = current_sum
    return dic_areas

def from_profile_get_depths(dic_areas, n_slices):
    values = np.array(list(dic_areas.values()))
    area_tot = values[-1]
    li_depths = []
    for i in range(n_slices):
        area_searched = ((i + 0.5)/n_slices)*area_tot
        to_minimize = abs(values - area_searched)
        index = np.argmin(to_minimize)
        good_depth = list(dic_areas.keys())[index]
        li_depths.append(good_depth)
    return li_depths

def from_depths_get_slices(li_depths, image):
    #on suposera que l'image et le mask sont bien à l'échelle à laquelle on veut les échantillonner
    li_slice_num = []
    origin = image.GetOrigin()
    for i, depth in enumerate(li_depths):
        slice_num = int(image.TransformPhysicalPointToContinuousIndex([origin[0],origin[1],depth])[2])
        li_slice_num.append(slice_num)
    return li_slice_num
    


        
        
    
    
    
    
if __name__ == '__main__':    
    ls_image =sorted(glob.glob('./*_NAT*/*.nii'))
    image = sitk.ReadImage(ls_image[0])
    mask = sitk.ReadImage(ls_image[0].replace('_NAT','').replace('.nii','_masked.nii')) != 0
    size_phys  =(np.array(image.GetSize())-1)*np.array(image.GetSpacing())
    origin = np.array(image.GetOrigin())
    new_spacing = np.array(image.GetSpacing())/2
    new_image, resampler  =resample_image_to_reference(image, image, new_spacing, force_size = None)
    new_mask = resampler.Execute(mask)
    dic_areas = get_profile(image, mask, new_spacing/2)
    li_depths = from_profile_get_depths(dic_areas, n_slices = 10)
    li_slices_num = from_depths_get_slices(li_depths, new_image)
    print(li_slices_num)
    for i in range(10):
        slice_num = li_slices_num[i]
        slice_image = new_image[:,:,slice_num]
        slice_mask = new_mask[:,:,slice_num]
        mask_superpose_simple(slice_image,slice_mask, num = i + 1, max_num = 5,line = 2)
    plt.show()

    #pour l'instant slices [55, 60, 64, 66, 69, 71, 73, 75, 77, 79]. Bien avec new_spacing mais faux avec spacing initial

    # point_z = image.TransformContinuousIndexToPhysicalPoint([0,0,36])[2]
    # slice_num = int(image.TransformPhysicalPointToContinuousIndex([0,0,point_z])[2])
    # slice_image = image[:,:,slice_num]
    # slice_mask = mask[:,:,slice_num]
    # mask_superpose_simple(slice_image,slice_mask, num = 1, max_num = 2)
    # print(origin, new_image.GetOrigin())    
    # slice_num = int(new_image.TransformPhysicalPointToContinuousIndex([origin[0],origin[1],point_z])[2])
    # print(new_image.TransformPhysicalPointToContinuousIndex([origin[0],origin[1],point_z]))
    # slice_image = new_image[:,:,slice_num]
    # slice_mask = new_mask[:,:,slice_num]
    # mask_superpose_simple(slice_image,slice_mask, num = 2, max_num = 2)
    # plt.show()

