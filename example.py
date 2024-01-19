import SPremiRNA as sm
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import pandas as pd


outdir = "./prostate_data/Prostate"
cancer = "Prostate"
results_folder = outdir+'_analysis'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

if __name__ == "__main__":

#Step 1 运行spatial_process
    # adata = sc.read_visium(
    #     path='./Spatial Cancer/' + cancer + "/",
    #     library_id=cancer,
    #     count_file='Visium_FFPE_Human_Prostate_IF_filtered_feature_bc_matrix.h5')
    # adata.obs_names_make_unique()
    # adata.var_names_make_unique()
    #
    # sm.QC(adata)
    # sm.QC_filter(adata,1500,30000)
    # # 数据标准化
    # sm.Standard(adata)
    # adata = sm.Select_highly(adata)
    # # 数据归一化
    # sc.pp.scale(adata, max_value=10)
    # tcga_mrna = pd.read_csv('./tcga_cancers/mrna/tcga_mrna_PRAD.csv', index_col=0).T
    # tcga_mirna = pd.read_csv('./tcga_cancers/mirna/tcga_mirna_PRAD.csv', index_col=0).T
    # ccle_mrna = pd.read_csv('./ccle_cancers/mrna/ccle_mrna_PROSTATE.csv', index_col=0).T
    # ccle_mirna = pd.read_csv('./ccle_cancers/mirna/ccle_mirna_PROSTATE.csv', index_col=0).T
    # sm.Save_ranked_mrna(tcga_mrna,ccle_mrna,adata,outdir)
    # sm.Save_ranked_mirna(tcga_mirna,ccle_mirna,outdir)
    # sm.Spatial_analyze(adata, outdir)
    # sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
    # sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")

# Step 2 运行spatial_xgboost
#     mrna = pd.read_csv(outdir + '_ranked_tcga_ccle_mrna.csv', index_col=0)
#     mirna = pd.read_csv(outdir + '_ranked_tcga_ccle_mirna.csv', index_col=0)
#     mrna, mirna = sm.standard(mrna, mirna)
#     x_train, x_test, y_train, y_test = sm.SplitData(mrna, mirna, split_size=0.1)
#     print(x_train.shape, x_test.shape, y_train.shape, y_test.shape)
#     model = sm.Build_Model(lr=0.05,n_estimators=300, max_depth=4)
#     sm.Train_Model(model, x_train, y_train,x_test,y_test,cancer,outdir)
#     st_mrna = pd.read_csv(outdir+"_st_mrna.csv", index_col=0)
#     st_mirna = model.predict(st_mrna)
#     st_mirna = pd.DataFrame(st_mirna, index=st_mrna.index, columns=y_train.columns)
#     st_mirna.to_csv(outdir + "_Pre" + "st_mirna.csv")

#Step 3 运行cell2loc
    # adata_ref = sc.read_loom('./prostate_data/PRAD_GSE141445_expression.loom')
    # meta = pd.read_table( "./prostate_data/PRAD_GSE141445_CellMetainfo_table.tsv")
    # print(adata_ref)
    # adata_ref = sm.Build_ref(adata_ref, meta, sample_counts=31000)
    # adata_ref, mod = sm.SC_estimation(adata_ref, "Celltype",ref_run_name=ref_run_name)
    # mod.plot_QC()
    # # 重建准确性，以评估推断是否存在任何问题。这个2D直方图图应该沿着有噪声的对角线具有大多数观测值。
    # # export estimated expression in each cluster
    # if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    #     inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
    #                                                           for i in adata_ref.uns['mod']['factor_names']]].copy()
    # else:
    #     inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
    #                               for i in adata_ref.uns['mod']['factor_names']]].copy()
    # inf_aver.columns = adata_ref.uns['mod']['factor_names']
    #
    # adata_vis = sc.read_visium(
    #     path='./Spatial Cancer/' + cancer + "/",
    #     library_id=cancer,
    #     count_file='Visium_FFPE_Human_Prostate_IF_filtered_feature_bc_matrix.h5')
    #
    # adata_vis, inf_aver = sm.PreprocessST(adata_vis, inf_aver)
    #
    # adata_vis, mod = sm.ST_estimation(adata_vis, inf_aver,run_name)
    # mod.plot_QC()
    # fig = mod.plot_spatial_QC_across_batches()
    # sm.Visual_abundance(adata_vis, adata_ref)

#Step 4 运行malig_mirna
    # adata_vis = sc.read_h5ad(outdir + "_analysis/cell2location_map/sp.h5ad")
    # adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    #
    # sp_mrna = sc.read_h5ad(outdir + '_mrna.h5ad')
    # pre_mirna = pd.read_csv(outdir + "_Pre" + "st_mirna.csv", index_col=0)
    # sp_mirna = sm.Build_STmirna(pre_mirna, sp_mrna)
    # sp_mirna = sm.Add_celltype(adata_vis, sp_mirna, sp_mrna,outdir)
    # sc.pl.umap(sp_mirna,
    #            color="celltype", legend_loc='on data',
    #            frameon=False, legend_fontsize=8, legend_fontoutline=2
    #            )
    #
    # sc.pl.spatial(sp_mirna, img_key="hires", color="celltype", size=1.5)
    # sc.pl.spatial(sp_mirna, img_key="hires", color="clusters", size=1.5)
    # top10 = sm.Find_top(sp_mirna, "Malignant",outdir)
    # print(top10)
    # sc.pl.rank_genes_groups_heatmap(sp_mirna, n_genes=3, groupby="celltype")


#Step 5 运行cor_analyze
