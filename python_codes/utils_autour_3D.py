from __future__ import print_function
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
    image = sitk.Cast(image, sitk.sitkFloat32)
    mean_filter = sitk.SmoothingRecursiveGaussianImageFilter()
    mean_filter.SetSigma(sigma)
    mean_image = mean_filter.Execute(image)
    
    square_filter  =sitk.SquareImageFilter()
    # Créer un filtre de convolution pour calculer la moyenne des carrés locaux
    square_image = square_filter.Execute(image)
    mean_square_image = mean_filter.Execute(square_image)
    squared = square_filter.Execute(mean_image)
    
    # Calculer la variance locale
    variance_image = mean_square_image - squared
    # print(np.max(sitk.GetArrayFromImage(image)) - np.max(sitk.GetArrayFromImage(mean_image)))
    # print(np.mean(sitk.GetArrayFromImage(image)) - np.mean(sitk.GetArrayFromImage(mean_image)))
    #print(np.mean(sitk.GetArrayFromImage(square_image)) - np.mean(sitk.GetArrayFromImage(squared)))
    # print(np.max(sitk.GetArrayFromImage(image)))
    # print(np.max(sitk.GetArrayFromImage(square_image)))
    # plt.imshow(sitk.GetArrayFromImage(image)[110,:,:])
    # plt.colorbar()
    # plt.show()
    # plt.imshow(sitk.GetArrayFromImage(square_image)[110,:,:])
    # plt.colorbar()
    # plt.show()
    # plt.imshow(sitk.GetArrayFromImage(mean_square_image)[110,:,:])
    # plt.colorbar()
    # plt.show()
    # plt.imshow(sitk.GetArrayFromImage(mean_image)[110,:,:])
    # plt.colorbar()
    # plt.show()
    # plt.imshow(sitk.GetArrayFromImage(squared)[110,:,:])
    # plt.colorbar()
    # plt.show()
    return variance_image

def kill_to_high(image_to_high, zone_high,image_target, quantile = 0.3,show = False, slice= None):
    mask_filter = sitk.MaskImageFilter()
    np_image = sitk.GetArrayFromImage(image_to_high)
    np_target = sitk.GetArrayFromImage(image_target)
    np_zone_high = sitk.GetArrayFromImage(zone_high)
    quantile_limit_sup = np.quantile(np_image[np_zone_high == 1], quantile)
    mask_low = image_to_high <= quantile_limit_sup
    masked_image =  mask_filter.Execute(image_target, mask_low) #on conserve seulement le masque
    masked_image = mask_filter.Execute(masked_image, zone_high == 1) #on conserve seulement le masque
    quantile_limit_low = 0.15* np.quantile(np_target.flatten(),0.9)+ 0.85 *np.min(np_target.flatten())
    #quantile_limit_high = 0.2* np.quantile(np_target.flatten(),0.1)+ 0.8*np.max(np_target.flatten())
    mask_high = image_target > quantile_limit_low
    #mask_low_lum = image_target < quantile_limit_high
    masked_image = mask_filter.Execute(masked_image, mask_high) #on conserve seulement le masque
    #masked_image = mask_filter.Execute(masked_image, mask_low_lum) 
    # plt.imshow(sitk.GetArrayFromImage(mask_high)[36,:,:])
    # plt.show()
    #return masked_image, (mask_low + zone_high + mask_high + mask_low_lum) == 4
    if show:
        plt.imshow(np_image[slice,:,:])
        plt.colorbar()
        plt.show()
    return masked_image, (mask_low + zone_high + mask_high) == 3

def get_connected_mask(mask, erosion = 8):
    mask = sitk.Cast(mask != 0, sitk.sitkUInt8)
    filter_erode = sitk.BinaryErodeImageFilter()
    filter_erode.SetForegroundValue(1)
    filter_erode.SetKernelRadius(erosion)
    mask_eroded = filter_erode.Execute(mask)
    connected_component_filter = sitk.ConnectedComponentImageFilter()
    connected_compos = connected_component_filter.Execute(mask_eroded)
    relabel_filter = sitk.RelabelComponentImageFilter()
    connected_compos_sorted = relabel_filter.Execute(connected_compos)
    new_mask_eroded = sitk.Cast(connected_compos_sorted == 1, sitk.sitkUInt8)
    filter_dilate = sitk.BinaryDilateImageFilter()
    filter_dilate.SetForegroundValue(1)
    filter_dilate.SetKernelRadius(erosion)
    new_mask = sitk.Cast(filter_dilate.Execute(new_mask_eroded) & mask, sitk.sitkUInt8)
    return new_mask
    
    
def pipeline_3D(image, mask, sigma = 1, quantile = 0.2, interieur = 5, exterieur = 15, erosion = 8, show = False, slice = None):
    variance = local_variance(image, sigma = sigma)
    zone_high= extrac_circle(image, mask, interieur = interieur, exterieur = exterieur)
    _ , mask_lum = kill_to_high(variance, zone_high,image, quantile = quantile, show = show, slice = slice)
    mask_tot = get_connected_mask(mask_lum, erosion = erosion)
    return mask_tot

    