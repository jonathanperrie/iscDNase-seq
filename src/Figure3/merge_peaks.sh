#!/bin/bash

# Merges peaks into groups based on assay and cell type 

# Merged single cell
cat sc_atac_b_cell_full_peaks.bed sc_atac_mono_full_peaks.bed sc_atac_t_cell_full_peaks.bed i_dnase_b_cell_full_peaks.bed i_dnase_mono_full_peaks.bed i_dnase_t_cell_full_peaks.bed > merged_sc_peaks.bed

# Merged assay
cat b_atac_*_sub_peaks.bed > b_atac_peaks.bed
cat b_dnase_*_sub_peaks.bed > b_dnase_peaks.bed
cat sc_atac_*_sub_peaks.bed > sc_atac_peaks.bed
cat i_dnase_*_sub_peaks.bed > i_dnase_peaks.bed

# Merged cell type
cat sc_atac_b_cell_sub_peaks.bed i_dnase_b_cell_sub_peaks.bed b_atac_b_cell_sub_peaks.bed b_dnase_b_cell_sub_peaks.bed > merged_b_cell_peaks.bed
cat sc_atac_t_cell_sub_peaks.bed i_dnase_t_cell_sub_peaks.bed b_atac_t_cell_sub_peaks.bed b_dnase_t_cell_sub_peaks.bed > merged_t_cell_peaks.bed
cat sc_atac_mono_sub_peaks.bed i_dnase_mono_sub_peaks.bed b_atac_mono_sub_peaks.bed b_dnase_mono_sub_peaks.bed > merged_mono_peaks.bed

