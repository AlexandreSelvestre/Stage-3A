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

def extrac_circle(image, mask , interieur = 5, exterieur = 15):
    mask = sitk.Cast(mask != 0, sitk.sitkUInt8)
    filter_int  = sitk.BinaryDilateImageFilter()
    filter_ext  = sitk.BinaryDilateImageFilter()
    filter_int.SetForegroundValue(1)
    filter_ext.SetForegroundValue(1)
    filter_int.SetKernelRadius(interieur)
    filter_ext.SetKernelRadius(exterieur)
    mask_exte  = filter_ext.Execute(mask)
    mask_inte  = filter_int.Execute(mask)
    mask_circle  = mask_exte - mask_inte
    return mask_circle
    
    
def local_variance(image, sigma = 1.5):
    # mean_filter = sitk.MeanImageFilter()
    # mean_filter.SetRadius(radius)
    mean_filter = sitk.SmoothingRecursiveGaussianImageFilter()
    mean_filter.SetSigma(sigma)
    mean_image = mean_filter.Execute(image)
    
    square_filter  =sitk.SquareImageFilter()
    # Créer un filtre de convolution pour calculer la moyenne des carrés locaux
    square_image = square_filter.Execute(image)
    mean_square_image = mean_filter.Execute(square_image)
    
    # Calculer la variance locale
    variance_image = mean_square_image - sitk.Square(mean_image)
    
    return variance_image

def kill_to_high(image_to_high, zone_high,image_target, quantile = 0.3,show = False):
    mask_filter = sitk.MaskImageFilter()
    np_image = sitk.GetArrayFromImage(image_to_high)
    np_target = sitk.GetArrayFromImage(image_target)
    np_zone_high = sitk.GetArrayFromImage(zone_high)
    quantile_limit_sup = np.quantile(np_image[np_zone_high == 1], quantile)
    mask_low = image_to_high <= quantile_limit_sup
    masked_image =  mask_filter.Execute(image_target, mask_low) #on conserve seulement le masque
    masked_image = mask_filter.Execute(masked_image, zone_high == 1) #on conserve seulement le masque
    quantile_limit_low = 0.15* np.quantile(np_target.flatten(),0.9)+ 0.85 *np.min(np_target.flatten())
    mask_high = image_target > quantile_limit_low
    masked_image = mask_filter.Execute(masked_image, mask_high) #on conserve seulement le masque
    # plt.imshow(sitk.GetArrayFromImage(mask_high)[36,:,:])
    # plt.show()
    return masked_image, (mask_low + zone_high + mask_high) == 3


def get_connected_mask(mask):
    mask = sitk.Cast(mask != 0, sitk.sitkUInt8)
    connected_component_filter = sitk.ConnectedComponentImageFilter()
    connected_compos = connected_component_filter.Execute(mask)
    relabel_filter = sitk.RelabelComponentImageFilter()
    connected_compos_sorted = relabel_filter.Execute(connected_compos)
    new_mask = sitk.Cast(connected_compos_sorted == 1, sitk.sitkUInt8)
    return new_mask
    
    
def pipeline_3D(image, mask):
    variance = local_variance(image, sigma = 1)
    zone_high= extrac_circle(image, mask)
    _ , mask_lum = kill_to_high(variance, zone_high,image, quantile = 0.2)
    mask_tot = get_connected_mask(mask_lum)
    return mask_tot

    