
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
#print("3", datetime.now())
import tqdm #dejafe
import numpy as np
import itertools
#print("4", datetime.now())
from utils import mask_superpose_simple, resample_image_to_reference
from utils_advanced_area import *
if __name__ == '__main__':
    current_dir = os.getcwd()
    import matplotlib
    if current_dir != "/gpfs/users/selvestra/python_codes":
        matplotlib.use('TKAgg')
from utils_recaler import eliminate_temps
from natsort import natsorted
    

#print("5", datetime.now())
#arrive après interpolate_all_reciprocals
def kill_depassement(dic_interpolators, li_interpolations, best_li_decals):
    #li_interpolations va de l'espace vers l'aire non normalisée
    #ici area est normalisé!!!
    min_depth_global = - np.inf
    max_depth_global = np.inf
    for image_num in dic_interpolators.keys():
        dic_interpolator = dic_interpolators[image_num]
        min_depth = dic_interpolator['min_depth']
        max_depth = dic_interpolator['max_depth']
        min_depth_global = max(min_depth_global, min_depth)
        max_depth_global = min(max_depth_global, max_depth)
    for image_num in dic_interpolators.keys():
        renorm_constant = dic_interpolators[image_num]["renorm_constant"]
        decal = best_li_decals[image_num]
        area_to_kill_min = li_interpolations[image_num]["interpolator"].__call__(min_depth_global - decal)/renorm_constant
        area_to_kill_max = li_interpolations[image_num]["interpolator"].__call__(max_depth_global - decal)/renorm_constant
        dic_interpolators[image_num]["area_to_kill_min"] = area_to_kill_min
        dic_interpolators[image_num]["area_to_kill_max"] = area_to_kill_max
    return dic_interpolators

def get_final_interpolator(dic_interpolators, image_num = None):
    #ici area est normalisé!!!
    if image_num is None:
        def mean_interpolator(area):
            depth= 0
            for image_num in dic_interpolators.keys():
                interpolator = dic_interpolators[image_num]["interpolator"]
                area_to_kill_min = dic_interpolators[image_num]["area_to_kill_min"]
                area_to_kill_max = dic_interpolators[image_num]["area_to_kill_max"]
                depth += interpolator.__call__((1 - area)*area_to_kill_min + area*area_to_kill_max) #on commence à l'aire qui garantit 
            depth = depth/len(dic_interpolators.keys())
            return depth
        mean_interpolator = np.vectorize(mean_interpolator)
        return mean_interpolator
    else:
        def local_interpolator(area):
            interpolator = dic_interpolators[image_num]["interpolator"]
            area_to_kill_min = dic_interpolators[image_num]["area_to_kill_min"]
            area_to_kill_max = dic_interpolators[image_num]["area_to_kill_max"]
            depth = interpolator.__call__((1 - area)*area_to_kill_min + area*area_to_kill_max) 
            return depth
        local_interpolator = np.vectorize(local_interpolator)
        return local_interpolator
        
            

def get_n_depths_equalized(dic_interpolators, best_li_decals,n_slices, do_average = True,show = False):
    dic_depths_to_extract = {}
    dic_securite = {}
    if do_average:
        mean_interpolator = get_final_interpolator(dic_interpolators)
    else:
        dic_final_interpolators = {image_num:get_final_interpolator(dic_interpolators, image_num) for image_num in dic_interpolators.keys()}
    areas = np.linspace(0, ((n_slices-1)/n_slices), n_slices)
    areas = areas + (0.5/n_slices)
    for image_num in dic_interpolators.keys():
        if do_average:
            depths_to_extract = mean_interpolator(areas)
        else:
            local_interpolator = dic_final_interpolators[image_num]
            depths_to_extract = local_interpolator(areas)
        securite_high = dic_interpolators[image_num]["securite_high"]
        securite_low = dic_interpolators[image_num]["securite_low"]
        dic_securite[image_num] = {"securite_high": securite_high, "securite_low": securite_low}
        new_depths_to_extract = np.copy(depths_to_extract) - best_li_decals[image_num]
        #print(best_li_decals[image_num])
        dic_depths_to_extract[image_num] = np.copy(new_depths_to_extract)
        if show:
            area_show = np.linspace(0, 1, 100)
            if do_average:
                depth_show = mean_interpolator(area_show)
            else:
                depth_show = local_interpolator(area_show)
            if (do_average and image_num == list(dic_interpolators.keys())[0]) or (not do_average):
                plt.plot(area_show, depth_show)
                plt.scatter(areas, depths_to_extract)
    if show:
        plt.show()
    return dic_depths_to_extract, dic_securite 





if __name__ == '__main__':
    #print("6", datetime.now())
    path_data_brut = "../data/radio_brut"
    current_dir = os.getcwd()
    show = True
    if current_dir == "/gpfs/users/selvestra/python_codes":
        path_data_brut = "/gpfs/workdir/selvestra/data/radio_brut"
        show = False
    path_radios = './*_NAT*/*.nii'
    path_brut = os.path.join(path_data_brut, path_radios)
    path_data = os.path.dirname(path_data_brut)
    path_save = os.path.join(
        path_data, "multislice_excel_with_shape_2D_entre_deux.xlsx")
    ls_image = sorted(glob.glob(path_brut))
    time_inj = ["_ART", "_PORT","_TARD"]
    #time_inj = ["_ART", "_PORT","_VEIN","_TARD"]
    ls_image_prime = ls_image
    ls_image_no_time = natsorted(list(set([eliminate_temps(x) for x in ls_image_prime])))
    ls_image_full_times_surconfiance = [[name +time + ".nii" for time in time_inj] for name in ls_image_no_time]
    ls_image_full_time = []
    dict_image_full_time = {}
    for li_names in ls_image_full_times_surconfiance:
        li_true_names = []
        for i,name in enumerate(li_names):
            if name in ls_image_prime:
                li_true_names.append(name)
        ls_image_full_time.append(li_true_names)
        classe_name = li_names[0].split('/')[-2].split('_')[1]
        patient_num = li_names[0].split('/')[-1].split('_')[0]
        dict_image_full_time[(patient_num, classe_name)] = li_true_names
    
    #print("7", datetime.now())
    num = 5
    #print(ls_image_full_time[num])
    li_images = [sitk.ReadImage(image_name) for image_name in ls_image_full_time[num]]
    li_masks = [sitk.ReadImage(image_name.replace('_NAT','').replace('.nii','_masked.nii')) != 0 for image_name in ls_image_full_time[num]]
    # print(li_images[0].GetOrigin(),li_images[0].TransformContinuousIndexToPhysicalPoint([0,0,0]) )
    li_interpolations, best_li_decals = recaler_interpolation(li_images, li_masks, discretisation = 3)
    #print("best decals:",best_li_decals)
    dic_points, dic_securite = generate_points_initial_space(li_interpolations, best_li_decals, show = False) 
    dic_interpolators = interpolate_all_reciprocals(dic_points, dic_securite, renorm  = True) 
    if show:
        print("ggo")
        plot_reciprocals(dic_interpolators)
    print("8", datetime.datetime.now())
    dic_interpolators = kill_depassement(dic_interpolators, li_interpolations, best_li_decals) 
    depths_to_extract, dic_securite = get_n_depths_equalized(dic_interpolators, best_li_decals, 10, do_average = False,show = show) #change
    #print("depths", depths_to_extract)
    dic_slices_to_extract = get_n_slices(li_images, li_masks, depths_to_extract, dic_securite) 
    #print(dic_slices_to_extract)
    print("9", datetime.datetime.now())