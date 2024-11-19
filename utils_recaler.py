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
from natsort import natsorted
from utils import get_all_good_slices, give_ls, get_good_slice, rescale_image, rescale_image_float, mask_superpose_simple, check_key_num, resample_image_to_reference, get_good_slices_from_li


def eliminate_temps(name):
    useless_part_length = len(name.split('_')[-1]) + 1
    return (name[:-useless_part_length])


def equalize_slices(ls_image_full_time, dict_image_full_time, force_spacing = None, num=None, key=None, show=False, do_nothing=False):
    names = check_key_num(num, key, ls_image_full_time, dict_image_full_time)
    names_masks = [name.replace('_NAT', '').replace(
        '.nii', '_masked.nii') for name in names]
    li_masks = [(sitk.ReadImage(name_mask) != 0) for name_mask in names_masks]
    li_images = [sitk.ReadImage(name) for name in names]
    li_num_slice_tot = np.array([mask.GetSize()[2] for mask in li_masks])
    reference = max(li_num_slice_tot)
    differences = reference - li_num_slice_tot
    former_diff = np.copy(differences)
    if show:
        print("difference avant", differences)
    good_time = np.argmin(differences)
    if not do_nothing:
        for i in range(len(li_images)):
            li_images[i] = resample_image_to_reference(
                li_images[i], li_images[good_time], force_spacing=force_spacing)
            li_masks[i] = resample_image_to_reference(
                li_masks[i], li_images[good_time], force_spacing=force_spacing, interpolator = "nearest")


    li_num_slice_tot = np.array([mask.GetSize()[2] for mask in li_masks])
    li_slices_masks = [[sitk.Cast(mask[:, :, z], sitk.sitkUInt8) for z in range(li_num_slice_tot[i])] for i, mask in enumerate(li_masks)]
    li_slices_images = [[image[:, :, z] for z in range(
        li_num_slice_tot[i])] for i, image in enumerate(li_images)]
    reference = max(li_num_slice_tot)
    differences = reference - li_num_slice_tot
    if not do_nothing:
        for i, li_slices_time_t in enumerate(li_slices_masks):
            difference = differences[i]
            if difference > 0:
                print("TRES ETRANGE, new differences:",
                    differences, "avant on avait", former_diff)
                raise Exception("Le resampled n'a pas la bonne!")
    return li_slices_masks, li_slices_images, li_images, li_masks


def shift_vertical_with_zeros(mat, shift=1):
    shift = - shift
    if shift == 0:
        return mat
    elif shift > 0:
        if shift >= mat.shape[0]:
            return np.zeros_like(mat)
        return np.vstack([np.zeros((shift, mat.shape[1])), mat[:-shift]])
    else:
        shift = abs(shift)
        if shift >= mat.shape[0]:
            return np.zeros_like(mat)
        return np.vstack([mat[shift:], np.zeros((shift, mat.shape[1]))])


def shift_horizontal_with_zeros(mat, shift=1):
    if shift == 0:
        return mat
    elif shift > 0:
        if shift >= mat.shape[1]:
            return np.zeros_like(mat)
        return np.hstack([np.zeros((mat.shape[0], shift)), mat[:, :-shift]])
    else:
        shift = abs(shift)
        if shift >= mat.shape[1]:
            return np.zeros_like(mat)
        return np.hstack([mat[:, shift:], np.zeros((mat.shape[0], shift))])


def move(li_slices_time_t, delta_x, delta_y, delta_z):
    li_slices_time_t = [shift_horizontal_with_zeros(
        slice, delta_x) for slice in li_slices_time_t]
    li_slices_time_t = [shift_vertical_with_zeros(
        slice, delta_y) for slice in li_slices_time_t]
    if delta_z > 0:
        li_slices_time_t = li_slices_time_t[delta_z:]
        li_slices_time_t = li_slices_time_t + \
            [np.zeros_like(li_slices_time_t[0]) for i in range(delta_z)]
    if delta_z < 0:
        li_slices_time_t = li_slices_time_t[:delta_z]
        li_slices_time_t = [np.zeros_like(
            li_slices_time_t[0]) for i in range(-delta_z)] + li_slices_time_t

    return li_slices_time_t


def calc_ps(li_slices_time_t_1, li_slices_time_t_2, mode="simple"):
    ps = 0
    # print(len(li_slices_time_t_1), len(li_slices_time_t_2))
    if mode == "simple":
        for i in range(len(li_slices_time_t_1)):
            ps += np.sum(np.bitwise_and(li_slices_time_t_1[i].astype(
                bool), li_slices_time_t_2[i].astype(bool)))
    elif mode == "area":
        for i in range(len(li_slices_time_t_1)):
            ps += np.sum(li_slices_time_t_1[i]) * np.sum(li_slices_time_t_2[i])
    else:
        raise Exception("No mode of ps calculus")

    return ps


def calc_ps_decal(li_slices_time_t_ref, li_slices_time_t_2, delta_x, delta_y, delta_z, mode="simple"):
    li_slices_time_t_2 = move(li_slices_time_t_2, delta_x, delta_y, delta_z)
    # print(delta_x, delta_y, delta_z)
    return calc_ps(li_slices_time_t_ref, li_slices_time_t_2, mode=mode)


def find_best_decal(li_slices_time_t_ref, li_slices_time_t_2, fenetre_x, fenetre_y, fenetre_z):
    ps_max = 0
    best_decal = [0, 0, 0]
    for delta_x in fenetre_x:
        for delta_y in fenetre_y:
            for delta_z in fenetre_z:
                ps = calc_ps_decal(
                    li_slices_time_t_ref, li_slices_time_t_2, delta_x, delta_y, delta_z)
                # print(ps)
                if ps > ps_max:
                    ps_max = ps
                    best_decal = [delta_x, delta_y, delta_z]
    return best_decal


def find_best_z_decal(li_slices_time_t_ref, li_slices_time_t_2, fenetre_z):
    ps_max = 0
    best_decal = 0
    for delta_z in fenetre_z:
        ps = calc_ps_decal(li_slices_time_t_ref,
                           li_slices_time_t_2, 0, 0, delta_z, mode="area")
        # print(ps)
        if ps > ps_max:
            ps_max = ps
            best_decal = delta_z
    return best_decal


def find_best_decal_all_times(ls_image_full_time, dict_image_full_time, fenetre_z, fenetre_x=None, fenetre_y=None, num=None, key=None, mode="area"):
    names = check_key_num(num, key, ls_image_full_time, dict_image_full_time)
    # def resample_image_to_reference(image, reference_image):
    li_slices, _, _, _ = equalize_slices(
        ls_image_full_time, dict_image_full_time, num=num, key=key)
    li_slices_array  = [[sitk.GetArrayFromImage(slice) for slice in li_slices_time] for li_slices_time in li_slices]
    li_slices_time_t_ref = li_slices_array[0]
    if mode == "simple":
        li_best_decal = [[0, 0, 0]] + [find_best_decal(
            li_slices_time_t_ref, li_slices_time_t_other, fenetre_x, fenetre_y, fenetre_z) for li_slices_time_t_other in li_slices_array[1:]]
    if mode == "area":
        # vecteur des différentes décalages vs premioer temps
        li_best_decal = [0] + [find_best_z_decal(li_slices_time_t_ref, li_slices_time_t_other, fenetre_z)
                               for li_slices_time_t_other in li_slices_array[1:]]
    return li_best_decal


def find_best_decal_all_times_short(li_slices, fenetre_z, fenetre_x=None, fenetre_y=None, mode="area"):
    # def resample_image_to_reference(image, reference_image):
    li_slices_array  = [[sitk.GetArrayFromImage(slice) for slice in li_slices_time] for li_slices_time in li_slices]
    li_slices_time_t_ref = li_slices_array[0]
    if mode == "simple":
        li_best_decal = [[0, 0, 0]] + [find_best_decal(
            li_slices_time_t_ref, li_slices_time_t_other, fenetre_x, fenetre_y, fenetre_z) for li_slices_time_t_other in li_slices_array[1:]]
    if mode == "area":
        # vecteur des différentes décalages vs premioer temps
        li_best_decal = [0] + [find_best_z_decal(li_slices_time_t_ref, li_slices_time_t_other, fenetre_z)
                               for li_slices_time_t_other in li_slices_array[1:]]
    return li_best_decal

def find_best_all_area(li_masks, fenetre_z, fenetre_x, fenetre_y, ref_num = 0):
    li_slices_x = [[sitk.GetArrayFromImage(li_slices_time_t[x,:,:]) for x in range(li_slices_time_t.GetSize()[0])] for li_slices_time_t in li_masks]
    li_slices_y = [[sitk.GetArrayFromImage(li_slices_time_t[:,y,:]) for y in range(li_slices_time_t.GetSize()[1])] for li_slices_time_t in li_masks]
    li_slices_z = [[sitk.GetArrayFromImage(li_slices_time_t[:,:,z]) for z in range(li_slices_time_t.GetSize()[2])] for li_slices_time_t in li_masks]
    if fenetre_x is not None:
        li_best_decal_x = [find_best_z_decal(li_slices_x[ref_num], li_slices_time_t_other, fenetre_x) if i != ref_num else 0 for i,li_slices_time_t_other in enumerate(li_slices_x)]
    else:
        li_best_decal_x = [0 for i in range(len(li_masks))]
    if fenetre_y is not None:
        li_best_decal_y = [find_best_z_decal(li_slices_y[ref_num], li_slices_time_t_other, fenetre_y) if i != ref_num else 0 for i, li_slices_time_t_other in enumerate(li_slices_y)]
    else:
        li_best_decal_y = [0 for i in range(len(li_masks))]
    li_best_decal_z = [find_best_z_decal(li_slices_z[ref_num], li_slices_time_t_other, fenetre_z) if i != ref_num else 0 for i, li_slices_time_t_other in enumerate(li_slices_z)]
    spacing = li_masks[0].GetSpacing()
    li_best_shift_y = [li_best_shift_y * spacing[1] for li_best_shift_y in li_best_decal_y]
    li_best_shift_x = [li_best_shift_x * spacing[0] for li_best_shift_x in li_best_decal_x]
    return li_best_shift_x, li_best_shift_y, li_best_decal_z

def translate_mask(mask, x_shift, y_shift, z_shift):
    #rappel la translation sitk va de l'output vers l'input
    #Mais ici double inversion avec le fait que ps shift calcule aussi du nouveau temps vers la référence
    translation = sitk.TranslationTransform(3) #nb dimensions
    translation.SetOffset((-x_shift, -y_shift,-z_shift))
    default_value = 0
    resampled = sitk.Resample(mask,mask, translation,sitk.sitkNearestNeighbor , default_value)
    return resampled
    
if __name__ == "__main__":
    import matplotlib
    matplotlib.use('TkAgg')
    ls_image =sorted(glob.glob('./*_NAT*/*.nii'))
    name = ls_image[0]
    mask = sitk.ReadImage(name.replace('_NAT', '').replace('.nii', '_masked.nii')) !=0
    plt.subplot(1, 2, 1)
    plt.imshow(sitk.GetArrayFromImage(mask)[36,:,:])
    new_mask = translate_mask(mask, - 100, 10, 0)
    plt.subplot(1, 2, 2)
    plt.imshow(sitk.GetArrayFromImage(new_mask)[36,:,:])
    plt.show()
    
    
def find_new_scale(li_images, li_masks,n_slices, median_spacing):
    # retourne les spacings sur z qu'on voudrait pour que chaque tumeur soit de n_slices
    li_spacings_z = []
    li_new_num_slices = []
    li_heights_true_size = []
    for i, image in enumerate(li_images):
        spacing_small_z = [median_spacing[0], median_spacing[1], min(median_spacing[2]/2, image.GetSpacing()[2]/2)] #pour être fin sur z
        image_small_z , resampler = resample_image_to_reference(image, image, spacing_small_z)
        mask_small_z = resampler.Execute(li_masks[i]) != 0
        
        li_image_loc = [image_small_z[:,:,z] for z in range(image_small_z.GetSize()[2])]
        li_mask_loc = [mask_small_z[:,:,z] for z in range(mask_small_z.GetSize()[2])]
        li_true_slices_loc = get_good_slices_from_li(li_image_loc, li_mask_loc)
        height_pixels = li_true_slices_loc[len(li_true_slices_loc) - 1] - li_true_slices_loc[0] #pas de +1 car aucune cance d'être PILE au début de la tumeur sur le première slice
        spacing_z = image_small_z.GetSpacing()[2]
        height_true_size = height_pixels*spacing_z
        better_spacing_z = height_true_size/n_slices #doit garantir n_slices (car on a nécessairement un peu sous-estimé la hauteur de la tumeur à cause du sampling pas infiniment petit)
        # Il faut potentiellement augmenter le nombre de slices sinon la tumeur peut atterrir out of bounds
        new_num_slices = int(np.ceil(image.GetSize()[2]*image.GetSpacing()[2]/better_spacing_z))
        li_spacings_z.append(better_spacing_z)
        li_new_num_slices.append(new_num_slices)
        li_heights_true_size.append(height_true_size)
    return li_spacings_z, li_new_num_slices, li_heights_true_size
    
