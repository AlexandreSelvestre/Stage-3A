
from __future__ import print_function
from utils_advanced_recaler import *
from utils_advanced_area import *
from area_calculus import *
from utils_autour_3D import *
from utils_recaler import *
from utils import *
from utils_mask_sain import *
import numpy as np
import datetime
import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import glob
import SimpleITK as sitk
import os
import logging
from radiomics import featureextractor, getFeatureClasses
import radiomics
from natsort import natsorted


path_data_brut = "../data/radio_brut"
current_dir = os.getcwd()
if current_dir == "/gpfs/users/selvestra/python_codes":
    path_data_brut = "/gpfs/workdir/selvestra/data/radio_brut"
path_radios = './*_NAT*/*.nii'
path_brut = os.path.join(path_data_brut, path_radios)
path_data = os.path.dirname(path_data_brut)
path_save = os.path.join(
    path_data, "multislice_excel_with_shape_2D_5_area_tot.xlsx")
ls_image = sorted(glob.glob(path_brut))


time_inj = ["_ART", "_PORT", "_TARD"]
ls_image_prime = ls_image
ls_image_no_time = natsorted(
    list(set([eliminate_temps(x) for x in ls_image_prime])))
ls_image_full_times_surconfiance = [
    [name + time + ".nii" for time in time_inj] for name in ls_image_no_time]
ls_image_full_time = []
dict_image_full_time = {}
for li_names in ls_image_full_times_surconfiance:
    li_true_names = []
    for i, name in enumerate(li_names):
        if name in ls_image_prime:
            li_true_names.append(name)
    ls_image_full_time.append(li_true_names)
    classe_name = li_names[0].split('/')[-2].split('_')[1]
    patient_num = li_names[0].split('/')[-1].split('_')[0]
    dict_image_full_time[(patient_num, classe_name)] = li_true_names

ls_image_full_time = ls_image_full_time

median_spacing = [0.59375, 0.59375, 1.5]


settings = {}
settings['binWidth'] = 25
settings['resampledPixelSpacing'] = [2, 2]
settings['interpolator'] = 'sitkBSpline'
settings['verbose'] = True
settings['force2D'] = True
settings['force2Ddimension'] = 2


extractor = featureextractor.RadiomicsFeatureExtractor(**settings)
extractor.settings


mode = "area"
fenetre_x = None
fenetre_y = None
fenetre_z = range(-10, 10)
label = extractor.settings.get('label')
label = 1
minDims = 2
extractor.disableAllFeatures()
# utiliser des try et excepts ...
li_featureClass = ['firstorder', 'shape2D', 'glrlm', 'glszm', 'gldm', 'ngtdm']
li_carac_glcm = ['Autocorrelation',
                 'ClusterProminence',
                 'ClusterShade',
                 'ClusterTendency',
                 'Contrast',
                 'Correlation',
                 'DifferenceAverage',
                 'DifferenceEntropy',
                 'DifferenceVariance',
                 'Id',
                 'Idm',
                 'Idmn',
                 'Idn',
                 'Imc1',
                 'Imc2',
                 'InverseVariance',
                 'JointEnergy',
                 'JointEntropy',
                 'MCC',
                 'MaximumProbability',
                 'SumAverage',
                 'SumEntropy',
                 'SumSquares'
                 ]
for featureClass in li_featureClass:
    extractor.enableFeatureClassByName(featureClass)
extractor.enableFeaturesByName(**{"glcm": li_carac_glcm})
logger = logging.getLogger("radiomics.glcm")
logger.setLevel(logging.ERROR)
al_dat_slice = pd.DataFrame()


do_nothing = False  # Si ZERO RESAMPLING
no_decal = False  # SI ZERO DECALAGE
show_plots = False
n_slices = 5
keep_whole_curve = True
if keep_whole_curve is False:
    average_curve = False

for num, li_image_names in tqdm.tqdm(enumerate(ls_image_full_time)):
    print("num:",  num, "datetime:", datetime.datetime.now())
    # li_mask_names = [image_name.replace('_NAT','').replace('.nii','_masked.nii') for image_name in li_image_names]

    classe_name = li_image_names[0].split('/')[-2].split('_')[1]
    patient_num = li_image_names[0].split('/')[-1].split('_')[0]

    sizes = [np.array(sitk.ReadImage(image_name).GetSize())[2]
             for image_name in li_image_names]
    max_size = np.max(sizes)
    difference = max_size - sizes

    li_images = [sitk.ReadImage(image_name) for image_name in li_image_names]
    li_masks = [sitk.ReadImage(image_name.replace('_NAT', '').replace(
        '.nii', '_masked.nii')) != 0 for image_name in li_image_names]

    # même échelle selon x et y pour tout le monde.
    li_num_slices = np.array([image.GetSize()[2] for image in li_images])
    image_num_with_max_slices = np.argmax(li_num_slices)
    for i, image in enumerate(li_images):
        resampled_image = resample_image_to_reference(
            image, image, median_spacing[:2], interpolator="bspline")
        resampled_mask = resample_image_to_reference(
            li_masks[i] != 0, li_masks[i] != 0, median_spacing[:2], interpolator="nearest") != 0
        li_images[i] = resampled_image
        li_masks[i] = resampled_mask

    li_shapes = [list(image.GetSize()) for image in li_images]
    li_shapes_masks = [list(mask.GetSize()) for mask in li_masks]
    li_caled = []

    
    print("debut recalage",datetime.datetime.now())
    if keep_whole_curve:
        discretisation = 3
    else:
        discretisation = 99
    li_interpolations, best_li_decals = recaler_interpolation(
        li_images, li_masks, discretisation=discretisation)
    print("decals", best_li_decals)
    dic_points, dic_securite = generate_points_initial_space(
        li_interpolations, best_li_decals)
    if keep_whole_curve:
        dic_interpolators = interpolate_all_reciprocals(
            dic_points, dic_securite)
        # plot_reciprocals(dic_interpolators)
        depths_to_extract, dic_securite = get_n_depths(
            dic_interpolators, best_li_decals, n_slices)
        # print(depths_to_extract)
        dic_slices_to_extract = get_n_slices(
            li_images, li_masks, depths_to_extract, dic_securite)
    else:
        dic_interpolators = interpolate_all_reciprocals(
            dic_points, dic_securite, renorm=True)
        dic_interpolators = kill_depassement(
            dic_interpolators, li_interpolations, best_li_decals)
        depths_to_extract, dic_securite = get_n_depths_equalized(
            dic_interpolators, best_li_decals, n_slices, do_average=average_curve, show=False)
        dic_slices_to_extract = get_n_slices(
            li_images, li_masks, depths_to_extract, dic_securite)

    print("fin recalage",datetime.datetime.now())
    for image_num_local, image in enumerate(li_images):
        li_caled.append({})  # une liste par radio
        dic_slices_local = dic_slices_to_extract[image_num_local]
        spacing_z = dic_slices_local["spacing_z"]
        li_slices_num = dic_slices_local["li_slices"]
        mask = li_masks[image_num_local]
        spacing = np.copy(median_spacing)
        spacing[2] = spacing_z
        # print(spacing)
        image_sampled = resample_image_to_reference(
            image, image, spacing, interpolator="bspline")
        mask_sampled = resample_image_to_reference(
            mask != 0, mask != 0, spacing, interpolator="nearest") != 0
        for slice_num_ref in range(n_slices):
            slice_num = li_slices_num[slice_num_ref]
            slice_mask = mask_sampled[:, :, slice_num]  # format sitk
            slice_img = image_sampled[:, :, slice_num]  # format sitk
            lsif = sitk.LabelStatisticsImageFilter()
            lsif.Execute(slice_img, slice_mask)
            boundingBox = np.array(lsif.GetBoundingBox(label))
            # UBound - LBound + 1 = Size
            ndims = np.sum((boundingBox[1::2] - boundingBox[0::2] + 1) > 3)
            li_caled[-1][slice_num_ref] = {"img": slice_img, "mask": slice_mask,
                                           "accept": False, "time": time_inj[image_num_local], "slice_num_true": slice_num}
            imageName = li_image_names[image_num_local]
            if (sitk.GetArrayFromImage(slice_mask).sum() > 0) & (ndims >= minDims):
                try:
                    # Comprendre qui est accepté
                    li_caled[-1][slice_num_ref]["accept"] = True
                    featureVector = extractor.execute(slice_img, slice_mask)
                    # Tres important: on ne veut pas le slice_num mais le slice_num_ref. Et pas de slice 0 dans le dataframe
                    featureVector['slice_num'] = slice_num_ref + 1
                    df_idx = pd.DataFrame(
                        dict([(k, pd.Series(v)) for k, v in featureVector.items()])).iloc[0]
                    df_idx['classe_name'] = classe_name
                    temps_inj = imageName.split('_')[-1].split('.')[0]
                    df_idx['temps_inj'] = temps_inj
                    df_idx['patient_num'] = patient_num
                    # al_dat_slice = al_dat_slice.append(df_idx)
                    df_idx = pd.DataFrame(df_idx).T
                    al_dat_slice = pd.concat([al_dat_slice, df_idx])
                    # print(al_dat_slice)
                except Exception as e:
                    print("problem with slice", slice_num,
                          "of image", imageName, "erreur", e)


al_dat_slice.to_excel(
    path_save, index=False)
