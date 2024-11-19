
#from __future__ import print_function
import datetime
from scipy.interpolate import PchipInterpolator
from scipy.integrate import quad, simpson
import os
from radiomics import featureextractor, getFeatureClasses
import radiomics
import SimpleITK as sitk
import glob
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
import numpy as np
import logging
import itertools
from utils import mask_superpose_simple, resample_image_to_reference


if __name__ == '__main__':
    current_dir = os.getcwd()
    import matplotlib
    if current_dir != "/gpfs/users/selvestra/python_codes":
        matplotlib.use('TKAgg')
    from utils_recaler import eliminate_temps
    from natsort import natsorted
    
    
# Le parametre de securite assure qu'on extraira jamais une slice dans une zone sans masque (avec interpolation neighboor)
# Les paramètres true_min et true_max disent où prendre les intégrales (en esquivant la zone nulle)
# La min area assure de ne pas demander une aire plus basse que l'aire minimale (lié à la mauvaise construction de la courbe d'interpolation)
    
def get_profile_brut(image, mask):
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
    dic_areas = {}
    spacing = image.GetSpacing()
    current_sum = 0
    fin_croissance = False
    last_li_null_depths = []
    origin = image.GetOrigin()
    securite_low = origin[2]
    securite_high = origin[2] + spacing[2]*image.GetSize()[2]
    for i in range(image.GetSize()[2]):
        depth = image.TransformContinuousIndexToPhysicalPoint([0,0,i])[2]
        try:
            feature = extractor_area.execute(image[:,:,i], mask[:,:,i])
            area = feature['original_shape2D_MeshSurface']
            #print(area*spacing[2])
        except Exception as e:
            area = 0
        if area > 0 and current_sum == 0:
            securite_low = depth - 0.5*spacing[2]
        if area == 0 and current_sum > 0:
            fin_croissance = True
            last_li_null_depths.append(depth)
            securite_high = depth + 0.5*spacing[2]
        if area >0 and fin_croissance:
            for null_depth in last_li_null_depths:
                current_sum += previous_non_null_area*spacing[2]
                dic_areas[null_depth] = current_sum
            current_sum += area*spacing[2]
            dic_areas[depth] = current_sum
            last_li_null_depths = []
            fin_croissance = False
        else:
            if area > 0:
                previous_non_null_area = area
            current_sum = area*spacing[2] + current_sum
            dic_areas[depth] = current_sum
        #print("slice",i,"area",area, "sum",current_sum)
    print(securite_low)
    return dic_areas, securite_low, securite_high


def extend_out_of_bounds(f, a,b, to_null = False):
    def new_f(x):
        if x < a:
            if to_null:
                return 0
            else:
                return f(a)
        if x > b:
            if to_null:
                return 0
            else:
                return f(b)
        return f(x)
    return new_f

def interpolate_initial_space(dic_area, n_points = None):
    #On interpole dans domaine de départ. On sample dans ce domaine pour inverser. On re-interpole à larrivée. Superposer les deux graphes
    # si n_points est renseigné on retourne les depths associées aux areas. Sinon, on retourne une fonction
    depths = np.array(list(dic_area.keys()))
    areas = np.array(list(dic_area.values()))
    interpolator = PchipInterpolator(depths, areas, extrapolate = False)
    if n_points is None:
        return interpolator
    else:
        depths_x = np.linspace(depths[0], depths[-1], n_points)
        areas = interpolator.__call__(depths_x)
        return depths_x, areas


#CHecker les ordres de grandeur de décalages réels. On faisait avant au plus 10 slices de 1.5 (soit 15 dans chaque sens)

def recaler_interpolation(li_images,li_masks, discretisation = 99, time_ref = 0):
    # la discretisation regle le nombre de tentatives de décalages pour chaque courbe. Note: on autorise des décalages jusqu'à 10% de la taille totale de la radio
    li_interpolations = []
    mean_depth_min = 0
    mean_depth_max = 0
    borne_min = np.inf
    borne_max = - np.inf
    li_derivative_unbounded = [None for _ in range(len(li_images))]
    li_derivative_bounded = [None for _ in range(len(li_images))]
    for image_num, image in enumerate(li_images):
        mask = li_masks[image_num]
        dic_areas, securite_low, securite_high = get_profile_brut(image, mask)
        interpolator = interpolate_initial_space(dic_areas)
        min_depth=  min(list(dic_areas.keys()))
        max_depth = max(list(dic_areas.keys()))
        mean_depth_min += min_depth
        mean_depth_max += max_depth
        borne_min = min(borne_min, min_depth)
        borne_max = max(borne_max, max_depth)
        def create_derivative_unbounded(interpolator, min_depth, max_depth):
            def derivative_unbounded(x):
                if x < min_depth or x > max_depth:
                    return 0
                return interpolator.__call__(x, nu = 1)
            return derivative_unbounded
        li_derivative_unbounded[image_num] = create_derivative_unbounded(interpolator, min_depth, max_depth)
        li_derivative_unbounded[image_num]= np.vectorize(li_derivative_unbounded[image_num])
        x= np.linspace(min_depth, max_depth, discretisation*10)
        y = li_derivative_unbounded[image_num](x)
        # plt.plot(x,y, label = "image_num" + str(image_num))
        true_min = x[max(np.where(y > 0)[0][0] -1,0)]
        true_max = x[min(np.where(y > 0)[0][-1] + 1, len(x)-1)]
        integral_local, _=  quad(li_derivative_unbounded[image_num], a = true_min, b = true_max, limit = 1000, epsrel = 1e-5, epsabs = 20)
        # print("ici erreur", _, "integrale", integral_local)
        def create_divide(f, factor):
            def new_f(x):
                return f(x)/factor
            return new_f
        li_derivative_bounded[image_num] = create_divide(li_derivative_unbounded[image_num], integral_local)
        li_interpolations.append({"interpolator": interpolator, "min_depth": min_depth, "max_depth": max_depth, "true_min": true_min, "true_max": true_max, "securite_low": securite_low, "securite_high": securite_high, "depths": list(dic_areas.keys()), "areas": list(dic_areas.values())})
    # plt.legend()
    # plt.show()
    mean_depth_min = mean_depth_min/len(li_images)
    mean_depth_max = mean_depth_max/len(li_images)
    decalage_autorise = 0.1*(mean_depth_max - mean_depth_min) #en positif et en négatif
    borne_min = borne_min - decalage_autorise
    borne_max = borne_max + decalage_autorise
    if discretisation//2 == 0:
        discretisation += 1
    li_decalages_no_prod = [np.linspace(-decalage_autorise, decalage_autorise, discretisation) for _ in range(len(li_images)-1)]
    
    if time_ref is None:
        iterator = itertools.product(*li_decalages_no_prod) #produit cartésien de décalages
        best_result = -np.inf
        marge_best_int_low = -np.inf
        marge_best_int_high = - np.inf

        for li_decals in tqdm.tqdm(iterator):
            li_decals = list(li_decals)
            li_decals_augmented = [0] + li_decals
            def function_area(x):
                value = 1
                for i,delta_x in enumerate(li_decals_augmented):
                    value *= li_derivative_bounded[i](x - delta_x)
                return value
            
            function_area = np.vectorize(function_area)
            integral, error = quad(function_area,a= borne_min, b = borne_max, limit = 100, epsrel = 0.1e-5, epsabs = 0.1e-5)
            marge_new_int_high = integral + abs(error)
            marge_new_int_low = integral - abs(error)
            incertitude = False
            if integral > best_result:
                if marge_new_int_low <= marge_best_int_high:
                    incertitude = True
                else:
                    incertitude = False
                best_result = integral
                best_li_decals = li_decals_augmented
                marge_best_int_low = integral - abs(error)
                marge_best_int_high = integral + abs(error)
            elif marge_new_int_high >= marge_best_int_low:
                marge_best_int_high = max(marge_best_int_high, marge_new_int_high) #on sait jamais...
                incertitude = True
        if best_result <= 0:
            print("ouille, int nulle")
        if incertitude is True:
            print("incertitude")
    else:
        best_li_decals = [0 for _ in range(len(li_images))]
        for image_num in range(len(li_images)):
            if image_num == time_ref:
                pass
            else:
                best_result = -np.inf
                marge_best_int_low = -np.inf
                marge_best_int_high = - np.inf
                li_errors_relat = []
                for delta_x in li_decalages_no_prod[0]:
                    def function_area(x):
                        value = li_derivative_bounded[image_num](x - delta_x)* li_derivative_bounded[time_ref](x)
                        #print(li_derivative_bounded[time_ref](x), li_derivative_bounded[image_num](x - delta_x))
                        return value
                    function_area = np.vectorize(function_area)
                    x = np.linspace(borne_min, borne_max, discretisation* 10)
                    y = function_area(x)
                    if np.sum(y) > 0:
                        # true_min = x[np.where(y > 0)[0][0]]
                        # true_max = x[np.where(y > 0)[0][-1]]
                        true_min = x[max(np.where(y > 0)[0][0] -1,0)]
                        true_max = x[min(np.where(y > 0)[0][-1] + 1, len(x)-1)]
                    else:
                        true_min = borne_min
                        true_max = borne_max
                    li_x =np.linspace(true_min, true_max, 1000)
                    li_y = function_area(li_x)
                    integral = simpson(li_y, li_x)
                    li_x_autre = np.linspace(borne_min, borne_max, round(1000*0.9))
                    li_y_autre = function_area(li_x_autre)
                    integral_autre = simpson(li_y_autre, li_x_autre)
                    error = abs(integral - integral_autre)
                    if integral >0:
                        error_relat = error/integral
                    elif error == 0:
                        error_relat = 0
                    else:
                        error_relat = 0
                        print("il y a un inf mais on l'ignore pour delta_x =", delta_x)
                    li_errors_relat.append(error_relat)
                    
                    #print(delta_x, integral)
                    # y = function_area(np.linspace(borne_min, borne_max, 1000))
                    # plt.plot(np.linspace(borne_min, borne_max, 1000), y)
                    # plt.show()
                    marge_new_int_low = integral - abs(error)
                    marge_new_int_high = integral + abs(error)
                    if integral > best_result:
                        if marge_new_int_low <= marge_best_int_high:
                            incertitude = True
                        else:
                            incertitude = False
                        best_result = integral
                        best_li_decals[image_num] = delta_x
                        marge_best_int_low = integral - abs(error)
                        marge_best_int_high = integral + abs(error)
                    elif marge_new_int_high >= marge_best_int_low:
                        marge_best_int_high = max(marge_best_int_high, marge_new_int_high)
                        incertitude = True
                    #print(incertitude, integral)
                print("error_relat_moyenne", np.mean(np.array(li_errors_relat)))
                if best_result <= 0:
                    print("ouille, int nulle")
                if incertitude is True:
                    print("incertitude")
    return li_interpolations, best_li_decals


#Construire les courbes décalées: l'abscisse des depths_x est en équivalent dans l'espace de la première image

def generate_points_initial_space(li_interpolations, li_decals, n_points = 1000, show = False):
    dic_points = {}
    min_show = np.inf
    max_show = -np.inf
    dic_securite = {}
    for image_num in range(len(li_interpolations)):
        dic_securite[image_num] = {}
        interpolator = li_interpolations[image_num]["interpolator"]
        true_min = li_interpolations[image_num]["true_min"]
        true_max = li_interpolations[image_num]["true_max"]
        decal = li_decals[image_num]
        #decal = 0
        depths_x = np.linspace(true_min + decal, true_max + decal, n_points)
        if show:
            min_show = min(min_show, true_min)
            max_show = max(max_show, true_max)
        areas = interpolator.__call__(depths_x - decal)
        dic_points[image_num] = {"depths_x": depths_x, "areas": areas}
        dic_securite[image_num]["securite_low"] = li_interpolations[image_num]["securite_low"]
        dic_securite[image_num]["securite_high"] = li_interpolations[image_num]["securite_high"]
    #la sortie est dans l'espace de comparaison commun (premiere image)
    if show:
        depths_x_show = np.linspace(min_show, max_show, n_points)
        #for image_num in range(len(li_interpolations)):
        for image_num in range(1):
            #rien n'est recalé ici
            depths_points = li_interpolations[image_num]["depths"]
            areas_points = li_interpolations[image_num]["areas"]
            area_show = li_interpolations[image_num]["interpolator"].__call__(depths_x_show)
            #plt.plot(depths_x_show, area_show, label = "num_temps " + str(image_num))
            plt.plot(depths_x_show, area_show)
            plt.scatter(depths_points, areas_points)
        #plt.title("Plot of volume travelled in the tumor CCK 1 (time arterial) as a function of depth (along the z axis)")
        plt.xlabel("depths")
        plt.ylabel("volume travelled in the tumor")
        plt.legend()
        plt.show()
    return dic_points, dic_securite

def interpolate_reciprocal(depths_x, areas, renorm = False):
    #on est déjà ,recalés
    """renorm définit si on remet les aires entre 0 et 1 pour comparer bien les courbes"""
    
    index_to_interpolate  = np.where((areas > np.min(areas)) & (areas < np.max(areas)))[0]
    #print(index_to_interpolate)
    index_to_interpolate = index_to_interpolate.tolist()
    index_to_interpolate.append(max(index_to_interpolate)+1)
    index_to_interpolate.insert(0,min(index_to_interpolate)-1)
    index_to_interpolate= np.array(index_to_interpolate)
    depths_to_interpolate = depths_x[index_to_interpolate]
    areas_to_interpolate = areas[index_to_interpolate]
    renorm_constant = np.max(areas_to_interpolate)
    if renorm:
        areas_to_interpolate = areas_to_interpolate/renorm_constant
    #print(areas_to_interpolate)
    interpolator = PchipInterpolator(areas_to_interpolate, depths_to_interpolate, extrapolate = False)
    min_depth = np.min(depths_to_interpolate)
    max_depth = np.max(depths_to_interpolate)
    min_area = np.min(areas_to_interpolate)
    max_area = np.max(areas_to_interpolate)
    #print("basic", min_area, max_area)
    return interpolator, min_depth, max_depth, min_area, max_area, renorm_constant

def interpolate_all_reciprocals(dic_points, dic_securite, renorm = False): #change
    dic_interpolators = {}
    for image_num in dic_points.keys():
        depths_x = dic_points[image_num]["depths_x"]
        areas = dic_points[image_num]["areas"]
        interpolator, min_depth, max_depth, min_area, max_area, renorm_constant = interpolate_reciprocal(depths_x, areas, renorm = renorm)
        securite_low = dic_securite[image_num]["securite_low"]
        securite_high = dic_securite[image_num]["securite_high"]
        dic_interpolators[image_num] = {"interpolator": interpolator, "min_depth": min_depth, "max_depth": max_depth, "min_area": min_area, "max_area": max_area, "securite_low": securite_low, "securite_high": securite_high, "renorm_constant": renorm_constant}
    return dic_interpolators

def plot_reciprocals(dic_interpolators):
    for image_num in dic_interpolators.keys():
        interpolator = dic_interpolators[image_num]["interpolator"]
        min_area = dic_interpolators[image_num]["min_area"]
        max_area = dic_interpolators[image_num]["max_area"]
        areas = np.linspace(min_area, max_area, 1000)
        depths = interpolator.__call__(areas) + 0.001* image_num
        plt.plot(areas, depths, label = "image_num" + str(image_num))
    plt.legend()
    plt.show()
    
def get_n_depths(dic_interpolators, best_li_decals,n_slices):
    dic_depths_to_extract = {}
    dic_securite = {}
    for image_num in dic_interpolators.keys():
        securite_high = dic_interpolators[image_num]["securite_high"]
        securite_low = dic_interpolators[image_num]["securite_low"]
        dic_securite[image_num] = {"securite_high": securite_high, "securite_low": securite_low}
        interpolator = dic_interpolators[image_num]["interpolator"]
        min_area = dic_interpolators[image_num]["min_area"]
        max_area = dic_interpolators[image_num]["max_area"]
        areas = np.linspace(0, ((n_slices-1)/n_slices)*max_area, n_slices)
        areas = areas + (0.5/n_slices)*max_area
        areas = np.maximum(areas, min_area)
        # print("min_area", min_area)
        depths_to_extract = interpolator.__call__(areas)
        depths_to_extract = depths_to_extract - best_li_decals[image_num]
        #print(depths_to_extract)
        dic_depths_to_extract[image_num] = np.copy(depths_to_extract)
        #print(areas)
        if min_area >0:
            print("attention: min_area > 0 et vaut", min_area)
    return dic_depths_to_extract, dic_securite 

def get_n_slices(li_images, li_masks, dic_depths_to_extract, dic_securite):
    dic_slices_to_extract = {}
    for image_num in range(len(li_images)):
        image = li_images[image_num]
        mask = li_masks[image_num]
        depths_to_extract = dic_depths_to_extract[image_num]
        differences = depths_to_extract[1:] - depths_to_extract[:-1]
        min_diff = np.min(differences)
        spacing_z = max(min_diff/10, 0.2) # au cas où l'écrasement des aires aurait joué des tours
        securite_low = dic_securite[image_num]["securite_low"] + spacing_z
        securite_high = dic_securite[image_num]["securite_high"] - spacing_z
        if np.max(depths_to_extract) > securite_high:
            print("attention: max_depth > securite_high:", np.max(depths_to_extract), securite_high)
        if np.min(depths_to_extract) < securite_low:
            print("attention: min_depth < securite_low:", np.min(depths_to_extract), securite_low, "next depths:", depths_to_extract)
        depths_to_extract = np.maximum(depths_to_extract, securite_low)
        depths_to_extract = np.minimum(depths_to_extract, securite_high)
        new_spacing = np.array(image.GetSpacing())
        new_spacing[2] = spacing_z
        new_image = resample_image_to_reference(image, image, new_spacing, force_size = None, interpolator = "bspline")
        li_slices = []
        for depth in depths_to_extract:
            origin = new_image.GetOrigin()
            slice_num = int(new_image.TransformPhysicalPointToContinuousIndex([origin[0],origin[1],depth])[2])
            li_slices.append(slice_num)
        dic_slices_to_extract[image_num] = {"li_slices": li_slices, "spacing_z": spacing_z}
    return dic_slices_to_extract
        
    

    


if __name__ == '__main__':
    path_data_brut = "../data/radio_brut"
    current_dir = os.getcwd()
    if current_dir == "/gpfs/users/selvestra/python_codes":
        path_data_brut = "/gpfs/workdir/selvestra/data/radio_brut"
    path_radios = './*_NAT*/*.nii'
    path_brut = os.path.join(path_data_brut, path_radios)
    path_data = os.path.dirname(path_data_brut)
    path_save = os.path.join(
        path_data, "multislice_excel_with_shape_2D_entre_deux.xlsx")
    ls_image = sorted(glob.glob(path_brut))
    time_inj = ["_ART", "_PORT","_TARD"]
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
    
    num = 6
    li_images = [sitk.ReadImage(image_name) for image_name in ls_image_full_time[num]]
    li_masks = [sitk.ReadImage(image_name.replace('_NAT','').replace('.nii','_masked.nii')) != 0 for image_name in ls_image_full_time[num]]
    print("1",datetime.datetime.now())
    # print(li_images[0].GetOrigin(),li_images[0].TransformContinuousIndexToPhysicalPoint([0,0,0]) )
    li_interpolations, best_li_decals = recaler_interpolation(li_images, li_masks)
    print("best decals:",best_li_decals, "2",datetime.datetime.now())
    dic_points, dic_securite = generate_points_initial_space(li_interpolations, best_li_decals) #change!!
    dic_interpolators = interpolate_all_reciprocals(dic_points, dic_securite) #change!!
    #plot_reciprocals(dic_interpolators)
    depths_to_extract, dic_securite = get_n_depths(dic_interpolators, best_li_decals,10) #change!!
    print(depths_to_extract, "3", datetime.datetime.now())
    #get_n_slices prend un temps fou!!!!
    dic_slices_to_extract = get_n_slices(li_images, li_masks, depths_to_extract, dic_securite) #change
    print(dic_slices_to_extract, "4", datetime.datetime.now())

