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
import cv2
import scipy.ndimage as nd
from utils import get_all_good_slices, give_ls, get_good_slice, rescale_image, rescale_image_float, mask_superpose_simple

def cercle_image(image_name, slice_num = None,show = False):
    image = sitk.ReadImage(image_name)
    name_mask = image_name.replace('_NAT','').replace('.nii','_masked.nii')
    mask = sitk.ReadImage(name_mask)
    if slice_num is None:
        slice_num = get_good_slice(image_name)
    slice_good = sitk.GetArrayFromImage(image[:,:,slice_num])
    mask_slice_good = sitk.GetArrayFromImage(mask[:,:,slice_num]) !=0
    mask_slice_good = mask_slice_good.astype(np.uint8)
    mask_dilation_first = cv2.dilate(mask_slice_good, np.ones((3,3), np.uint8), iterations = 5) #10 et 18
    mask_dilation_second = cv2.dilate(mask_slice_good, np.ones((3,3), np.uint8), iterations = 15)
    # mask_dilation_3 = cv2.dilate(mask_slice_good, np.ones((3,3), np.uint8), iterations = 40)
    # mask_dilation_4 = cv2.dilate(mask_slice_good, np.ones((3,3), np.uint8), iterations = 50)
    mask_cercle = mask_dilation_second - mask_dilation_first # + mask_dilation_4 - mask_dilation_3
    
    #histogramme sur variance locale
    
    filtered = rescale_image_float(slice_good).astype(np.float64)
    size = 5
    mean = nd.uniform_filter(filtered, size)
    mean_sq = nd.uniform_filter(filtered**2, size)
    variance = (mean_sq - mean**2)
    #hist, bins = np.histogram(variance[mask_cercle == 1].flatten(), bins = "doane")
    hist, bins = np.histogram(variance[mask_cercle == 1].flatten(), bins = 100)
    #On a les bins maintenant on classe les points
    shape_init = slice_good.shape
    bin_indices = np.digitize(variance.flatten(), bins) - 1
    bin_indices = bin_indices.reshape(shape_init)
    bin_indices[mask_cercle == 0] = -1
    
    #On choisit les bins qu'on garde de l'histogramme: à côté les unes des autres
    pic = np.argmax(hist)
    pic = 0 ##################### ATTENTION C'est car on cherche le minimum de variance
    deja_pris = [pic]
    taille_max = np.sum(mask_cercle == 1)*0.2
    taille = hist[pic]
    nbins = len(hist)
    while taille < taille_max:
        possible = []
        if max(deja_pris) + 1< nbins :
            possible.append(max(deja_pris) + 1)
        if min(deja_pris) > 0:
            possible.append(min(deja_pris) - 1)
        values = [hist[i] for i in possible]
        new_bar = possible[np.argmax(values)]
        deja_pris.append(new_bar)
        taille = taille + hist[new_bar]
    
    bool_matrix = np.isin(bin_indices, deja_pris)
    
    quantile_low_lum =0.5* np.quantile(slice_good.flatten(),0.9)+ 0.5 *np.min(slice_good.flatten()) #pas être dans le noir complet
    cercle_lumineux = bool_matrix & (slice_good > quantile_low_lum)
    
    
    # quantile_inf_cercle = np.quantile(slice_good[mask_cercle == 1].flatten(),0.5)
    # quantile_sup_tumeur = np.quantile(slice_good[mask_slice_good == 1].flatten(),0.7)
    # if quantile_inf_cercle < np.quantile(slice_good[mask_slice_good == 1].flatten(),0.3):
    #     cercle_lumineux = (slice_good > quantile_inf_cercle) & (mask_cercle == 1) & (slice_good < quantile_sup_tumeur)
        
    # else:
    #     cercle_lumineux = (mask_cercle == 1) & (slice_good > quantile_inf_cercle)
    cercle_lumineux = cercle_lumineux.astype(np.uint8)
    if show:
        plt.figure(figsize=(20, 5))
        plt.subplot(1, 3, 1)
        plt.imshow(bool_matrix)
        plt.subplot(1, 3, 2)
        plt.imshow(variance)
        plt.subplot(1, 3, 3)
        plt.hist(variance[mask_cercle == 1].flatten(), bins=bins, edgecolor='black')
        plt.title("Histogramme des valeurs des pixels")
        plt.xlabel("Valeur des pixels")
        plt.ylabel("Fréquence")
        plt.show()
        mask_superpose_simple(image_name, slice_num = slice_num, other_mask = mask_cercle)
        mask_superpose_simple(image_name, slice_num = slice_num, other_mask = cercle_lumineux)
    return cercle_lumineux, slice_good, mask_slice_good, mask_cercle

def compute_score(i,components, slice_good, mask_slice_good ,mask_cercle_lumineux, intens_standard, coeff = 0):
    #intens_standard doit standardiser l'écart à l'intensité tumorale. Doit être calculé à partir du max des composantes lumineuses de taille significative.
    component_i = (components == i).astype(np.uint8)
    intens_median_tumor = np.quantile(slice_good[mask_slice_good == 1].flatten(), 0.5)
    size_cercle = np.sum(mask_cercle_lumineux)
    size_relat_component_i = np.sum(component_i)/size_cercle
    intens_compnent_i = np.quantile(slice_good[component_i == 1].flatten(), 0.5)
    ecart_intens_relat_component_i = np.abs(intens_compnent_i - intens_median_tumor)/ intens_standard
    score = size_relat_component_i - coeff*ecart_intens_relat_component_i
    #print(size_relat_component_i, coeff*ecart_intens_relat_component_i)
    return score
    
def compute_standard(components, num_comp, slice_good, mask_slice_good):
    max_intens = np.quantile(slice_good[mask_slice_good == 1].flatten(), 0.5)
    for i in range(1, num_comp+1):
        max_intens = max(max_intens, np.quantile(slice_good[components == i].flatten(), 0.5))
    return max_intens

def cercle_image_compos_connexe(image_name, slice_num = None,show = False):
    if slice_num is None:
        slice_num = get_good_slice(image_name)
    mask_cercle_lumineux, slice_good, mask_slice_good, mask_cercle = cercle_image(image_name, slice_num= slice_num,show = False)
    components, num_comp = nd.label(mask_cercle_lumineux) #l'indexation des composantes commence à 1
    intens_standard = compute_standard(components, num_comp, slice_good, mask_slice_good)
    dict_sizes = {str(i) : compute_score(i,components, slice_good, mask_slice_good ,mask_cercle_lumineux, intens_standard) for i in range(1, num_comp + 1)}
    dict_sizes = dict(sorted(dict_sizes.items(), key=lambda item: item[1], reverse = True))
    size_to_cover  = min(np.sum(mask_cercle) * 0.05, np.sum(mask_cercle_lumineux) * 0.5)
    covered_size = 0
    li_keys_compos_connexe = []
    i = 0
    while covered_size < size_to_cover:
        key = list(dict_sizes.keys())[i]
        li_keys_compos_connexe.append(key)
        covered_size += np.sum(components == int(key))
        i += 1
    if(len(li_keys_compos_connexe) >0):
        li_keys_compos_connexe = [li_keys_compos_connexe[0]]
    else:
        print("echec mask vide image", image_name, "slice",slice_num)
        return np.zeros_like(mask_cercle_lumineux)
    mask_cercle_compos_connexe = np.zeros_like(mask_cercle_lumineux)
    for key in li_keys_compos_connexe:
        mask_cercle_compos_connexe += (components == int(key)).astype(np.uint8)
    if show:
        mask_superpose_simple(image_name, slice_num = slice_num, other_mask = mask_cercle_lumineux)
        mask_superpose_simple(image_name, slice_num = slice_num, other_mask = mask_cercle_compos_connexe)
    return mask_cercle_compos_connexe
    
        
    