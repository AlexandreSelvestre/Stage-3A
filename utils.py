from __future__ import print_function
import scipy.ndimage as nd
import cv2
import numpy as np
import tqdm
import matplotlib.pyplot as plt
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
import matplotlib


def rescale_image(image):
    # Normaliser les valeurs de l'image entre 0 et 1
    image_min = np.min(image)
    image_max = np.max(image)
    image_normalized = (image - image_min) / (image_max - image_min)

    # Redimensionner les valeurs normalisées entre 0 et 255
    image_rescaled = (image_normalized * 255).astype(np.uint8)

    return image_rescaled


def rescale_image_float(image):
    min_val = np.min(image)
    max_val = np.max(image)
    return (image - min_val) / (max_val - min_val)


def get_good_slice(image_name):
    image = sitk.ReadImage(image_name)
    name_mask = image_name.replace('_NAT', '').replace('.nii', '_masked.nii')
    mask = sitk.ReadImage(name_mask)
    n_slices = image.GetSize()[2]
    li_good = []
    for slice_num in range(n_slices):
        slice_try = image[:, :, slice_num]
        mask_slice_try = mask[:, :, slice_num]
        lsif = sitk.LabelStatisticsImageFilter()
        lsif.Execute(slice_try, mask_slice_try != 0)
        boundingBox = np.array(lsif.GetBoundingBox(label=1))
        # UBound - LBound + 1 = Size
        ndims = np.sum((boundingBox[1::2] - boundingBox[0::2] + 1) > 3)
        if (sitk.GetArrayFromImage(mask_slice_try != 0).sum() > 0) & (ndims >= 2):
            li_good.append(slice_num)
    if len(li_good) == 0:
        print('No good slices')
        return
    else:
        return li_good[len(li_good)//2]
    
def get_good_slices_from_li(li_image,li_mask):
    li_good = []
    for slice_num, slice_image in enumerate(li_image):
        slice_mask = li_mask[slice_num]
        lsif = sitk.LabelStatisticsImageFilter()
        lsif.Execute(slice_image, slice_mask != 0)
        boundingBox = np.array(lsif.GetBoundingBox(label=1))
        ndims = np.sum((boundingBox[1::2] - boundingBox[0::2] + 1) > 3)
        if (sitk.GetArrayFromImage(slice_mask != 0).sum() > 0) & (ndims >= 2):
            li_good.append(slice_num)
    if len(li_good) == 0:
        print('No good slices')
        return
    else:
        return li_good


def get_all_good_slices(image_name, image=None, mask=None):
    if image_name == None:
        pass
    else:
        image = sitk.ReadImage(image_name)
        name_mask = image_name.replace(
            '_NAT', '').replace('.nii', '_masked.nii')
        mask = sitk.ReadImage(name_mask)
    n_slices = image.GetSize()[2]
    li_good = []
    for slice_num in range(n_slices):
        slice_try = image[:, :, slice_num]
        mask_slice_try = mask[:, :, slice_num]
        lsif = sitk.LabelStatisticsImageFilter()
        lsif.Execute(slice_try, mask_slice_try != 0)
        boundingBox = np.array(lsif.GetBoundingBox(label=1))
        # UBound - LBound + 1 = Size
        ndims = np.sum((boundingBox[1::2] - boundingBox[0::2] + 1) > 3)
        if (sitk.GetArrayFromImage(mask_slice_try != 0).sum() > 0) & (ndims >= 2):
            li_good.append(slice_num)
    if len(li_good) == 0:
        print('No good slices')
        return
    else:
        return li_good


def give_ls(ls_image):
    def mask_superpose(image_num, slice_num=None, other_mask=None):
        # Changer ls_image
        if type(image_num) == str:
            image_name = image_num
        else:
            image_name = ls_image[image_num]
        image = sitk.ReadImage(image_name)
        name_mask = image_name.replace(
            '_NAT', '').replace('.nii', '_masked.nii')
        mask = sitk.ReadImage(name_mask)
        if slice_num is None:
            slice_num = get_good_slice(image_name)
        slice_good = sitk.GetArrayFromImage(image[:, :, slice_num])
        mask_slice_good = sitk.GetArrayFromImage(mask[:, :, slice_num]) != 0
        slice_rgb = rescale_image(np.stack([slice_good]*3, axis=-1))
        slice_rgb[:, :, 0:3:2] = 0  # Seul canal 1 accepté
        slice_rgb[mask_slice_good == 1, 0] = 255
        if other_mask is not None:
            mask_superpose = other_mask != 0
            slice_rgb[mask_superpose == 1, 2] = 255
        plt.figure(figsize=(15, 5))
        plt.subplot(1, 3, 1)
        plt.imshow(slice_rgb)
        plt.subplot(1, 3, 2)
        plt.imshow(slice_good)
        plt.subplot(1, 3, 3)
        plt.imshow(sitk.GetArrayFromImage(mask[:, :, slice_num]))
    return mask_superpose


def mask_superpose_simple(slice_image, slice_mask, other_mask=None, num=1, max_num=1, legend="No legend"):
    # Show juste la superposition et attend une slice prédéfinie ainsi qu'une image name
    if isinstance(slice_image, sitk.Image):
        slice_good = sitk.GetArrayFromImage(slice_image)
    else:
        slice_good = slice_image
    if isinstance(slice_mask, sitk.Image):
        mask_slice_good = sitk.GetArrayFromImage(slice_mask) != 0
    else:
        mask_slice_good = slice_mask
    slice_rgb = rescale_image(np.stack([slice_good]*3, axis=-1))
    slice_rgb[:, :, 0:3:2] = 0  # Seul canal 1 accepté
    slice_rgb[mask_slice_good == 1, 0] = 255
    if other_mask is not None:
        if isinstance(other_mask, sitk.Image):
            other_mask = sitk.GetArrayFromImage(other_mask)
        mask_superpose = other_mask != 0
        slice_rgb[mask_superpose == 1, 2] = 255
    plt.subplot(1, max_num, num)
    plt.imshow(slice_rgb)
    plt.title(legend)


def check_key_num(num, key, ls, d):
    if key == None:
        names = ls[num]
    elif num == None:
        names = d[key]
    else:
        print("Please give a key or a num")
        raise Exception("arrêt")
    if (key != None and num != None):
        print("Please do not give two imputs")
        return Exception("arrêt")
    return names


def resample_image_to_reference(image, reference_image, force_spacing):
    # Obtenir l'espacement de l'image de référence
    if force_spacing is not None:
        reference_spacing = tuple(np.copy(force_spacing))
    else:
        reference_spacing = reference_image.GetSpacing()

    # Obtenir la taille de l'image de référence
    reference_size_initial = reference_image.GetSize()
    
    #nouvelle référence size
    reference_size = [int(np.ceil(reference_size_initial[i] * (reference_image.GetSpacing()[i] / reference_spacing[i]))) for i in range(len(reference_size_initial))]

    # Obtenir l'origine de l'image de référence
    reference_origin = reference_image.GetOrigin()

    # Obtenir la direction de l'image de référence
    reference_direction = reference_image.GetDirection()

    # new_spacing = calculate_new_spacing(original_size, original_spacing, reference_size)
    new_spacing = reference_spacing
    new_spacing = list(new_spacing)
    # print("ref", reference_spacing)
    # Définir le filtre de resegmentation
    resample = sitk.ResampleImageFilter()
    # new_spacing[2] = new_spacing[2]*1
    new_spacing = tuple(new_spacing)
    resample.SetOutputSpacing(new_spacing)
    resample.SetSize(reference_size)
    resample.SetOutputOrigin(reference_origin)
    resample.SetOutputDirection(reference_direction)
    resample.SetInterpolator(sitk.sitkLinear)

    # Appliquer la resegmentation
    resampled_image = resample.Execute(image)

    return resampled_image, resample
