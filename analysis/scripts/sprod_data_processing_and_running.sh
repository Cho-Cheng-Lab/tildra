# Don't forget to activate the conda virtual environment
conda activate ./sprod_venv/

# input directory should contain:
- filtered_feature_bc_matrix.tar.gz:
	- "filtered_feature_bc_matrix" directory:
		- "barcodes.tsv.gz"
		- "features.tsv.gz"
		- "matrix.mtx.gz"
- spatial.tar.gz:
	- "spatial" directory:
		- aligned_fiducials.jpg
		- detected_tissue_image.jpg
		- scalefactors_json.json
		- spatial_enrichment.csv
		- tissue_hires_image.png
		- tissue_lowres_image.png
		- tissue_positions_list.csv
- sample.tif: tif image of the sample, should have the same name of the input directory: i.e. skin269, the tif image should be named "skin269.tif"

#######################################################################################################################

# run on all the samples
cd /media/kurdiabed/143F2A3651271F75/special_projects/8_single_cell_spatial_jeffrey/sprod_input_data
for f in *; do 
name=${f}; 
echo ${name}; 
python /media/kurdiabed/143F2A3651271F75/special_projects/8_single_cell_spatial_jeffrey/virtual_env/SPROD/sprod/data_preprocessing.py \
--input_type "visium" \
-v /media/kurdiabed/143F2A3651271F75/special_projects/8_single_cell_spatial_jeffrey/sprod_input_data/${f}/ \
-o /media/kurdiabed/143F2A3651271F75/special_projects/8_single_cell_spatial_jeffrey/sprod_outputs/sprod_counts/${name}; 
done

cd /media/pk3/143F2A3651271F75/special_projects/8_single_cell_spatial_jeffrey/sprod_outputs/sprod_counts
for f in *; do 
name=${f}; 
echo ${name}; 
python /media/pk3/143F2A3651271F75/special_projects/8_single_cell_spatial_jeffrey/virtual_env/SPROD/sprod.py \
--sprod_umap /media/pk3/143F2A3651271F75/special_projects/8_single_cell_spatial_jeffrey/sprod_outputs/sprod_counts/${f}/ \
/media/pk3/143F2A3651271F75/special_projects/8_single_cell_spatial_jeffrey/sprod_outputs/sprod_denoised_data/${name}
done
